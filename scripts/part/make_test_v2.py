"""TEST tau taggers + HLT ParT: the TEST cache with the p2nano v2 ParT added, row-aligned.

The TEST tau taggers (hltJetTauTagTEST, ...) are scored on the v1 events only, in the row order of
the cache version V50nano_170pre6HLT_TESTtautag. p2nano v2 reprocessed the same events (v1 objects
bit-identical) and adds ParT for every jet. This script makes the cache version
V50nano_170pre6HLTv2_TESTtautag:
  * every file of V50nano_170pre6HLT_TESTtautag as a symlink (renamed to the new version);
  * <version>_<sample>_hltAK4PuppiJetParT.parquet and <version>_<sample>_HLTEMU.parquet taken from
    the v2 cache (V50nano_170pre6HLTv2, scripts/part/part_wps.py) and put in the TEST row order by
    the event ID (run, luminosityBlock, event).
Events of the TEST rows missing in v2 (minimum bias: v2 has 1,998,400 of the 1,999,040 v1 events)
get no ParT jets and HLTEMU bits False; the number is printed and written to the summary.

TEST samples and their v1 rows:
  HHTo2B2Tau, DYLL_M50: all v1 rows; MinBiasF1000: v1 MinBias files 0-999;
  MinBiasL1Pass: the v1 rows of the TOPO L1-pass event list (as make_hlt_jet_tautag_test.py).

Checks (stop on failure): v1 event IDs unique and matching the v1 row count; for the mapped rows
the HLT and L1 bits of the TEST cache (v1) and the v2 cache agree event by event, up to the
few events per bit of the non-reproducible L1 emulation (at most max(10, 0.1%) per bit).

v1 event IDs: <ids-dir>/ids_<HH|DY|MB>_v1.npy (sorted v1 file order, scratch scan); v2:
<v2 cache>/<v2>_<sample>_EventID.parquet (scripts/part/event_ids.py).

    python3 scripts/part/make_test_v2.py --ids-dir <dir with ids_*_v1.npy>
"""

import argparse
import glob
import json
import os
import sys

import awkward as ak
import numpy as np
import pandas as pd

from menu_tools.caching.make_hlt_jet_tautag_test import (
    HTAUPU,
    L1PASS_EVENTS,
    file_events,
    read_rows,
)
from menu_tools.utils import constants

V1 = "V50nano_170pre6HLT_TESTtautag"
V2 = "V50nano_170pre6HLTv2"
OUT = "V50nano_170pre6HLTv2_TESTtautag"
SAMPLES = {  # TEST sample: (v2 sample, v1 id file key)
    "HHTo2B2Tau": ("HHTo2B2Tau", "HH"),
    "DYLL_M50": ("DYLL_M50", "DY"),
    "MinBiasF1000": ("MinBias", "MB"),
    "MinBiasL1Pass": ("MinBias", "MB"),
}
NEW = ["hltAK4PuppiJetParT", "HLTEMU"]


def fail(msg):
    print(f"CHECK FAILED: {msg}")
    sys.exit(1)


def key(ids: np.ndarray) -> np.ndarray:
    ids = ids.astype(np.int64)
    return (ids[:, 0] << 52) ^ (ids[:, 1] << 32) ^ ids[:, 2]


def v1_rows(test_sample: str, n_v1: int) -> np.ndarray:
    if test_sample in ("HHTo2B2Tau", "DYLL_M50"):
        return np.arange(n_v1)
    n_ev = file_events(os.path.join(HTAUPU, "skim_jet"), "minbias")
    if (n_ev < 0).any() or n_ev.sum() != n_v1:
        fail(f"minbias skim file counts {n_ev.sum()} != {n_v1} v1 rows")
    if test_sample == "MinBiasF1000":
        return np.arange(int(n_ev[:1000].sum()))
    offset = np.concatenate([[0], np.cumsum(n_ev)[:-1]])
    req = pd.read_parquet(L1PASS_EVENTS)[["file", "entry"]]
    return np.sort(offset[req.file.to_numpy()] + req.entry.to_numpy())


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--ids-dir", required=True)
    p.add_argument("--samples", nargs="+", default=list(SAMPLES))
    args = p.parse_args()
    c1 = os.path.join(constants.CACHE_DIR, V1)
    c2 = os.path.join(constants.CACHE_DIR, V2)
    co = os.path.join(constants.CACHE_DIR, OUT)
    os.makedirs(co, exist_ok=True)
    # symlinks to all TEST v1 files (resolved)
    for f in sorted(glob.glob(os.path.join(c1, f"{V1}_*"))):
        dst = os.path.join(co, os.path.basename(f).replace(V1, OUT, 1))
        if not os.path.lexists(dst):
            os.symlink(os.path.realpath(f), dst)
    summary, summary_bits = {}, {}
    for ts in args.samples:
        s2, idk = SAMPLES[ts]
        ids1 = np.load(os.path.join(args.ids_dir, f"ids_{idk}_v1.npy"))
        k1 = key(ids1)
        if len(np.unique(k1)) != len(k1):
            fail(f"duplicate v1 event IDs ({idk})")
        rows1 = v1_rows(ts, len(ids1))
        n_test = ak.metadata_from_parquet(os.path.join(c1, f"{V1}_{ts}_L1.parquet"))[
            "num_rows"
        ]
        if n_test != len(rows1):
            fail(f"{ts}: {len(rows1)} v1 rows but {n_test} rows in the TEST cache")
        e2 = ak.from_parquet(os.path.join(c2, f"{V2}_{s2}_EventID.parquet"))
        k2 = key(np.stack([ak.to_numpy(e2[b]) for b in ("run", "luminosityBlock", "event")], 1))
        order = np.argsort(k2)
        pos = np.searchsorted(k2[order], k1[rows1])
        pos = np.clip(pos, 0, len(k2) - 1)
        found = k2[order][pos] == k1[rows1]
        r2 = np.where(found, order[pos], -1)
        print(f"{ts}: {len(rows1)} TEST rows, {found.sum()} found in v2 {s2}, {(~found).sum()} missing")
        # bit-identity check on the mapped rows
        uniq = np.unique(r2[found])
        for rec in ("HLT", "L1"):
            a1 = ak.from_parquet(os.path.join(c1, f"{V1}_{ts}_{rec}.parquet"))
            a2 = read_rows(os.path.join(c2, f"{V2}_{s2}_{rec}.parquet"), uniq)
            idx = np.searchsorted(uniq, r2[found])
            bad = {}
            for fld in a1.fields:
                if fld not in a2.fields:
                    continue
                d = int((ak.to_numpy(a1[fld])[found] != ak.to_numpy(a2[fld])[idx]).sum())
                if d:
                    bad[fld] = d
            print(f"  {rec}: {len(a1.fields)} bits compared, differing: {bad or 'none'}")
            # the L1 emulation is not bit-reproducible between runs (production README):
            # a few events per bit differ; a wrong mapping would differ in a large fraction
            tol = max(10, int(1e-3 * found.sum()))
            if any(d > tol for d in bad.values()):
                fail(f"{ts} {rec} bits differ between v1 and v2 beyond {tol} events: {bad}")
            summary_bits.setdefault(ts, {})[rec] = bad
        for obj in NEW:
            src = os.path.join(c2, f"{V2}_{s2}_{obj}.parquet")
            a2 = read_rows(src, uniq)
            idx = np.searchsorted(uniq, r2[found])
            if obj == "HLTEMU":
                cols = {}
                for fld in a2.fields:
                    v = np.zeros(len(rows1), dtype=bool)
                    v[found] = ak.to_numpy(a2[fld])[idx]
                    cols[fld] = v
                out = ak.zip(cols)
            else:
                # missing events -> empty list: append one empty event and point to it
                empty_event = ak.unflatten(ak.flatten(a2[:0]), [0])
                full = ak.concatenate([a2, empty_event])
                sel = np.where(found, np.searchsorted(uniq, np.where(found, r2, uniq[0])), len(a2))
                out = full[sel]
            tgt = os.path.join(co, f"{OUT}_{ts}_{obj}.parquet")
            ak.to_parquet(out, tgt + ".tmp")
            os.replace(tgt + ".tmp", tgt)
            print(f"  wrote {tgt}: {len(out)} events")
        summary[ts] = {
            "rows": int(len(rows1)),
            "found_in_v2": int(found.sum()),
            "bits_differing_v1_v2": summary_bits.get(ts, {}),
        }
    out = os.path.join(co, f"{OUT}_mapping.json")
    if os.path.exists(out):
        with open(out) as f:
            summary = json.load(f) | summary
    with open(out, "w") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main()
