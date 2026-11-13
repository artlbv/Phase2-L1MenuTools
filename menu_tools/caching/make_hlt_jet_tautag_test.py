"""TEST, not official: the experimental HLT-jet tau taggers as a cached object.

Writes the objects `hltJetTauTagTEST` and `hltJetTauTagRegTEST`, i.e. the HLT AK4 PUPPI
jets with the out-of-fold scores of the final tau taggers of the hlt-tau-pu study
(topo-hackathon-2026, results/hlt_tau_pu/README.md, sections 10.8 and 11): the
multi-task DeepSets `mcttreg_jv_s0` (3 classes tau / minimum bias / hard-scatter jet,
trained with inclusive tt, plus a tau pT regression), the BDT `bdt_full` and the BDT
pT regressor `bdt_reg`. The taggers are not HLT objects of any menu; this is a TEST.

Everything goes into a separate cache version (default V50nano_170pre6HLT_TESTtautag,
so the plots land in outputs/<that version>/ and never mix with the official ones):
  <version>_<sample>_hltJetTauTagTEST.parquet   the tagger object, one entry per
      cached hltAK4PuppiJet (raw nano order -> sorted by `pt`). Fields:
        pt, eta, phi   tau candidate: pT^sig (|sum p| of the PF candidates within
                       dR < 0.1 of the axis, the tau energy estimate the report's
                       trigger cuts on) and the axis (leading charged PF candidate)
                       for scored jets; pt = -1 and the jet eta / phi otherwise
        sigPt          pT^sig (= pt)
        ptReg          DeepSets regressed tau pT = pT^sig * exp(r), mcttreg_jv_s0
        ptRegBDT       BDT regressed tau pT (bdt_reg)
        jetPt, jetEta, jetPhi, rawPt   corrected / raw jet
        bdt            BDT score (bdt_full)
        ds             DeepSets trigger score p_tau / (p_tau + p_MB) (mcttreg_jv_s0)
        dsA01          p_tau / (p_tau + p_MB + 0.1 p_HS), the score that keeps the
                       b-jet rejection (report section 10.8)
        dsPtau, dsPmb, dsPhs           DeepSets class probabilities
        dsOld          the earlier DeepSets (mc_jv, no tt), p_tau / (p_tau + p_MB);
                       -1 where it was not scored (MinBiasL1Pass)
        scored         1 for scored jets (corrected pT > 20 GeV, |eta| < 2.5)
        label          skim label (0 tau, 1 hard-scatter jet, 2 pileup / MB, 3 e/mu)
      Unscored jets have all scores, pTs and the label at -1.
  <version>_<sample>_hltJetTauTagRegTEST.parquet   the same jets and fields with
      pt = ptReg (sorted by it): the tau candidate with the regressed pT.
  <version>_<sample>_<object>.parquet for the other objects the plots need: symlinks
      to the source version, or, with --max-file / --events, the row-aligned subset of
      the events (e.g. sample MinBiasF1000: the minimum-bias files in which every event
      was scored; MinBiasL1Pass: the events of --events).

Mapping (file, entry, jet_idx) of the skim to the cache: the skim's file index is the
position in its sorted file list, the cache is built from the sorted glob of the same
directory, so cache row = sum of the events of the earlier files + entry. Checked
before anything is written (the script stops on any failure):
  1. identical sorted file lists (cache ntuple_path vs the skim's files_<s>.txt);
  2. per-file event counts of the skim (uproot num_entries) add up to the cache rows;
  3. every scored jet: the cached raw jet at [row][jet_idx] has the skim's raw pT and
     eta, and its corrected pT (make_hlt_jets_corr) equals the skim's;
  4. completeness: the number of cached jets with corrected pT > 20, |eta| < 2.5 in the
     covered events equals the number of scored jets.
With --events (MinBiasL1Pass), also: the skim of those events (skim_jet_topo) holds
exactly the requested events, and every event of the full sample with the --seed-bit
fired is among them (so rates ANDed with that bit are complete).

MinBiasL1Pass normalisation: the sample is the 31,786 minimum-bias events of the L1
menu OR the three L1 bbtautau TOPOs at +20 kHz (topo-hackathon-2026
results/topo_lep/skim/tau_request), out of the 1,999,040 events of the full MinBias
cache. Rates ANDed with an L1 bit inside that set (L1_pDoublePuppiTau52_52) are complete;
normalise them to the full sample, 1 event = 31,038 kHz / 1,999,040 = 15.5 Hz
(`purity.norm_events: 1999040`). The normalisation is also written to
<version>_<out-sample>_norm.json.

Usage (1 core, a few minutes per sample; MinBiasL1Pass ~10 min of EOS reads):
    python3 menu_tools/caching/make_hlt_jet_tautag_test.py --sample HHTo2B2Tau \\
        --skim-sample hh
    python3 menu_tools/caching/make_hlt_jet_tautag_test.py --sample MinBias \\
        --skim-sample minbias --max-file 1000 --out-sample MinBiasF1000
    python3 menu_tools/caching/make_hlt_jet_tautag_test.py --sample MinBias \\
        --skim-sample minbias --events-l1pass --out-sample MinBiasL1Pass
"""

import argparse
import glob
import json
import os
import sys
from typing import Optional

import awkward as ak
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml

from menu_tools.caching.cache_objects import sorted_ntuple_files
from menu_tools.caching.make_hlt_jets_corr import (
    JEC_DIR,
    JEC_L2,
    correct_jets,
    l2relative,
    load_jec_table,
)
from menu_tools.utils import constants

TARGET = "hltJetTauTagTEST"
TARGET_REG = "hltJetTauTagRegTEST"
SRC_JETS = "hltAK4PuppiJet"
HTAUPU = (
    "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/alobanov/phase2/menu/"
    "ntuples/hlt_tau_pu"
)
L1PASS_EVENTS = (
    "/eos/home-a/alobanov/SWAN_projects/L1TML/topo-hackathon-2026/results/topo_lep/"
    "skim/tau_request/minbias_events.parquet"
)
CACHING_CFG = os.path.join(
    "configs", "V50nano", "cache_objects", "caching_170pre6_HLT.yaml"
)
DS = "ds_mcttreg_jv_s0"  # final DeepSets
DS_OLD = "ds_mc_jv"  # earlier DeepSets (no tt), kept as `dsOld`
ALPHA = 0.1
SCORES = [
    os.path.join(HTAUPU, "scored_jet.parquet"),  # bdt_full (and other BDTs)
    os.path.join(HTAUPU, "deepsets", f"scores_{DS[3:]}.parquet"),
    os.path.join(HTAUPU, "deepsets", "scores_bdtreg.parquet"),
    os.path.join(HTAUPU, "deepsets", f"scores_{DS_OLD[3:]}.parquet"),
]
SCORES_L1PASS = [os.path.join(HTAUPU, "deepsets", "topo_minbias.parquet")]
# Objects the TEST plots read besides the tagger objects
OBJECTS = {
    "HHTo2B2Tau": [
        "GenVisTau",
        "L1nnPuppiTau",
        "hltHpsPFTau",
        "hltAK4PuppiJetCorr",
        "Pileup",
        "L1",
        "HLT",
    ],
    "DYLL_M50": ["L1nnPuppiTau", "hltHpsPFTau", "hltAK4PuppiJetCorr", "Pileup"],
    "MinBias": [
        "L1nnPuppiTau",
        "hltHpsPFTau",
        "hltAK4PuppiJetCorr",
        "Pileup",
        "L1",
        "HLT",
    ],
}
CAND_COLS = [
    "file",
    "entry",
    "jet_idx",
    "jet_pt",
    "jet_eta",
    "jet_cpt",
    "sig_pt",
    "axis_eta",
    "axis_phi",
    "label",
]
KEY = ["file", "entry", "jet_idx"]
JET_PT, JET_ETA = 20.0, 2.5


def fail(msg: str) -> None:
    print(f"ALIGNMENT CHECK FAILED: {msg}")
    sys.exit(1)


def file_events(skim_dir: str, skim_sample: str) -> np.ndarray:
    """Events per file index from the skim's <sample>_files_*.parquet (failed reads,
    n_events = -1, were re-run: the valid row is kept)."""
    fs = sorted(glob.glob(os.path.join(skim_dir, f"{skim_sample}_files_*.parquet")))
    df = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
    df = df[df.n_events >= 0].drop_duplicates("file")
    n = np.full(int(df.file.max()) + 1, -1, dtype=np.int64)
    n[df.file.to_numpy()] = df.n_events.to_numpy()
    return n


def read_rows(fname: str, rows: Optional[np.ndarray]) -> ak.Array:
    """Rows of a cache parquet (all if None), sorted ascending; reads only the row
    groups that hold them (one group at a time)."""
    if rows is None:
        return ak.from_parquet(fname)
    pf = pq.ParquetFile(fname)
    md = pf.metadata
    start = np.cumsum([0] + [md.row_group(i).num_rows for i in range(md.num_row_groups)])
    parts = []
    for i in range(md.num_row_groups):
        sel = rows[(rows >= start[i]) & (rows < start[i + 1])]
        if len(sel):
            parts.append(ak.from_parquet(fname, row_groups=[i])[sel - start[i]])
    return ak.concatenate(parts) if len(parts) > 1 else parts[0]


def load_candidates(
    skim_dir: str, skim_sample: str, scores: list[str], max_file: Optional[int]
) -> pd.DataFrame:
    fs = sorted(glob.glob(os.path.join(skim_dir, f"{skim_sample}_cand_*.parquet")))
    filt = [("file", "<", max_file)] if max_file is not None else None
    cand = pd.concat(
        [pd.read_parquet(f, columns=CAND_COLS, filters=filt) for f in fs],
        ignore_index=True,
    )
    if cand.duplicated(KEY).any():
        fail("duplicate (file, entry, jet_idx) in the skim")
    for path in scores:
        s = pd.read_parquet(path)
        s = s[s["sample"] == skim_sample]
        s = s.drop(columns=[c for c in ("sample", "parity", "label") if c in s])
        s = s.drop(columns=[c for c in s if c in cand and c not in KEY])
        if max_file is not None:
            s = s[s.file < max_file]
        if s.duplicated(KEY).any():
            fail(f"duplicate keys in {path}")
        n0 = len(cand)
        cand = cand.merge(s, on=KEY, how="inner", validate="one_to_one")
        if len(cand) != n0 or len(s) != n0:
            fail(f"{path}: {len(s)} scores for {n0} candidates, {len(cand)} matched")
        print(f"scores {os.path.basename(path)}: {len(s)} jets matched")
    return cand


def check_and_index(
    cand: pd.DataFrame, raw: ak.Array, row_of: dict, table: np.ndarray
) -> np.ndarray:
    """Checks 3 and 4; returns the flat index of each candidate in the cached jets.

    row_of: offset (cache row of each file's first event) and, with --events, the
    sorted cache rows of the selected events (raw holds only those rows).
    """
    f = cand.file.to_numpy()
    e = cand.entry.to_numpy()
    n_ev = row_of["n_ev"]
    if (f >= len(n_ev)).any() or (e >= n_ev[f]).any() or (e < 0).any():
        fail("entry outside the file's event range")
    row = row_of["offset"][f] + e
    if row_of.get("rows") is not None:
        pos = np.searchsorted(row_of["rows"], row)
        if (pos >= len(row_of["rows"])).any() or (row_of["rows"][pos] != row).any():
            fail("scored jets in events outside the selected events")
        row = pos
    counts = ak.to_numpy(ak.num(raw[f"{SRC_JETS}_pt"]))
    j = cand.jet_idx.to_numpy().astype(np.int64)
    if (j >= counts[row]).any():
        fail(f"{(j >= counts[row]).sum()} jet_idx beyond the cached jets")
    start = np.concatenate([[0], np.cumsum(counts)])
    flat = start[row] + j
    pt = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_pt"])).astype(np.float64)
    eta = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_eta"])).astype(np.float64)
    dpt = np.abs(pt[flat] - cand.jet_pt.to_numpy()) / np.maximum(pt[flat], 1e-3)
    deta = np.abs(eta[flat] - cand.jet_eta.to_numpy())
    bad = (dpt > 1e-5) | (deta > 1e-5)
    print(
        f"check 3: {len(cand)} scored jets, raw pT / eta mismatches: {bad.sum()} "
        f"(max rel. dpT {dpt.max():.1e}, max deta {deta.max():.1e})"
    )
    if bad.any():
        fail(f"{bad.sum()} jets with different pT / eta")
    cpt = pt * l2relative(pt, eta, table)
    dc = np.abs(cpt[flat] - cand.jet_cpt.to_numpy()) / cpt[flat]
    print(f"check 3: corrected pT max rel. difference {dc.max():.1e}")
    if dc.max() > 1e-4:
        fail("corrected pT differs")
    n_events = len(raw)
    n_cand_cache = int(((cpt > JET_PT) & (np.abs(eta) < JET_ETA)).sum())
    print(
        f"check 4: {n_cand_cache} cached jets with corrected pT > {JET_PT:g}, "
        f"|eta| < {JET_ETA} in {n_events} events; {len(cand)} scored"
    )
    if n_cand_cache != len(cand):
        fail("scored jets != candidate jets in the cache")
    return flat


def build_objects(cand: pd.DataFrame, raw: ak.Array, flat: np.ndarray, table):
    """The two tagger objects: pt = pT^sig (TARGET) and pt = regressed pT (TARGET_REG)."""
    counts = ak.num(raw[f"{SRC_JETS}_pt"])
    pt = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_pt"])).astype(np.float64)
    eta = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_eta"]))
    phi = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_phi"]))
    n = len(pt)

    def full(values, default=-1):
        out = np.full(n, default, dtype=np.float32)
        out[flat] = values
        return out

    def vsmb(p_tau, p_mb, p_hs=None, alpha=0.0):
        den = p_tau + p_mb + (alpha * p_hs if p_hs is not None else 0)
        return np.where(den > 0, p_tau / np.maximum(den, 1e-30), 0.0)

    p_tau = cand[DS].to_numpy()
    p_mb = cand[f"{DS}_pmb"].to_numpy()
    p_hs = cand[f"{DS}_phs"].to_numpy()
    if DS_OLD in cand:
        ds_old = vsmb(cand[DS_OLD].to_numpy(), cand[f"{DS_OLD}_pmb"].to_numpy())
    else:
        ds_old = np.full(len(cand), -1.0)
    fields = {
        "pt": full(cand.sig_pt.to_numpy()),
        "eta": eta.astype(np.float32),
        "phi": phi.astype(np.float32),
        "sigPt": full(cand.sig_pt.to_numpy()),
        "ptReg": full(cand[f"{DS}_ptcorr"].to_numpy()),
        "ptRegBDT": full(cand.bdt_reg_ptcorr.to_numpy()),
        "jetPt": (pt * l2relative(pt, eta, table)).astype(np.float32),
        "jetEta": eta.astype(np.float32),
        "jetPhi": phi.astype(np.float32),
        "rawPt": pt.astype(np.float32),
        "bdt": full(cand.bdt_full.to_numpy()),
        "ds": full(vsmb(p_tau, p_mb)),
        "dsA01": full(vsmb(p_tau, p_mb, p_hs, ALPHA)),
        "dsPtau": full(p_tau),
        "dsPmb": full(p_mb),
        "dsPhs": full(p_hs),
        "dsOld": full(ds_old),
        "scored": full(np.ones(len(cand)), 0).astype(np.int32),
        "label": full(cand.label.to_numpy()).astype(np.int32),
    }
    fields["eta"][flat] = cand.axis_eta.to_numpy()
    fields["phi"][flat] = cand.axis_phi.to_numpy()
    out = {}
    for target, ptf in ((TARGET, "sigPt"), (TARGET_REG, "ptReg")):
        f = dict(fields, pt=fields[ptf])
        arr = ak.unflatten(ak.zip({f"{target}_{k}": v for k, v in f.items()}), counts)
        out[target] = arr[ak.argsort(arr[f"{target}_pt"], axis=-1, ascending=False)]
    return out


def write(arr: ak.Array, target: str) -> None:
    ak.to_parquet(arr, target + ".tmp")
    os.replace(target + ".tmp", target)
    print(f"Wrote {target}: {len(arr)} events")


def check_events(args, offset: np.ndarray, n_ev: np.ndarray, src: str):
    """--events: the selected cache rows; the skim covers exactly these events; every
    event of the full sample with the seed bit is among them."""
    req = pd.read_parquet(args.events)[["file", "entry"]]
    if req.duplicated().any():
        fail("duplicate events in --events")
    if (req.entry.to_numpy() >= n_ev[req.file.to_numpy()]).any():
        fail("--events entry outside its file")
    rows = np.sort(offset[req.file.to_numpy()] + req.entry.to_numpy())
    ev = pd.concat(
        [
            pd.read_parquet(f, columns=["file", "entry"])
            for f in sorted(
                glob.glob(os.path.join(args.skim_dir, f"{args.skim_sample}_ev_*.parquet"))
            )
        ],
        ignore_index=True,
    ).drop_duplicates()
    n_both = len(req.merge(ev, on=["file", "entry"]))
    print(
        f"check events: {len(req)} requested, {len(ev)} in the skim, {n_both} in both"
    )
    if not (n_both == len(req) == len(ev)):
        fail("the skim does not hold exactly the requested events")
    l1 = os.path.join(src, f"{args.src_version}_{args.sample}_L1.parquet")
    bit = ak.to_numpy(ak.from_parquet(l1, columns=[args.seed_bit])[args.seed_bit])
    fired = np.flatnonzero(bit.astype(bool))
    inside = np.isin(fired, rows)
    print(
        f"check seed: {args.seed_bit} fires in {len(fired)} of {len(bit)} events of the "
        f"full sample, {inside.sum()} of them in the selected events"
    )
    if not inside.all():
        fail(f"{(~inside).sum()} {args.seed_bit} events missing from --events")
    return rows, {"n_seed_full": int(len(fired)), "n_total": int(len(bit))}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--sample", required=True, choices=list(OBJECTS))
    p.add_argument("--skim-sample", required=True, help="hh, dy, minbias")
    p.add_argument("--out-sample", default=None, help="default: --sample")
    p.add_argument("--max-file", type=int, default=None)
    p.add_argument(
        "--events",
        default=None,
        help="parquet (file, entry) of the events to keep (row subset of the sample)",
    )
    p.add_argument(
        "--events-l1pass",
        action="store_true",
        help="--events = the L1-menu-pass minimum-bias events, with their skim "
        "(skim_jet_topo) and scores (deepsets/topo_minbias.parquet)",
    )
    p.add_argument("--seed-bit", default="L1_pDoublePuppiTau52_52_final")
    p.add_argument("--src-version", default="V50nano_170pre6HLT")
    p.add_argument("--version", default="V50nano_170pre6HLT_TESTtautag")
    p.add_argument("--skim-dir", default=None, help="default: HTAUPU/skim_jet")
    p.add_argument("--file-list", default=None, help="default: skim/files_<s>.txt")
    p.add_argument("--scores", nargs="+", default=None)
    p.add_argument("--caching-cfg", default=CACHING_CFG)
    p.add_argument("--jec", default=os.path.join(JEC_DIR, JEC_L2))
    args = p.parse_args()
    if args.events_l1pass:
        args.events = args.events or L1PASS_EVENTS
        args.skim_dir = args.skim_dir or os.path.join(HTAUPU, "skim_jet_topo")
        args.scores = args.scores or SCORES_L1PASS
    args.skim_dir = args.skim_dir or os.path.join(HTAUPU, "skim_jet")
    args.scores = args.scores or SCORES
    if args.events and args.max_file is not None:
        p.error("--events and --max-file exclude each other")
    out_sample = args.out_sample or args.sample
    table = load_jec_table(args.jec)

    # 1. file lists
    with open(args.caching_cfg) as fh:
        cfg = yaml.safe_load(fh)
    cache_files = list(
        sorted_ntuple_files(cfg[args.src_version][args.sample]["ntuple_path"])
    )
    flist = args.file_list or os.path.join(
        HTAUPU, "skim", f"files_{args.skim_sample}.txt"
    )
    skim_files = open(flist).read().split()
    print(f"check 1: cache {len(cache_files)} files, skim {len(skim_files)} files")
    if cache_files != skim_files:
        fail("sorted file lists differ")

    # 2. events per file (the main skim has every file of the sample)
    n_ev = file_events(os.path.join(HTAUPU, "skim_jet"), args.skim_sample)
    src = os.path.join(constants.CACHE_DIR, args.src_version)
    jets_file = os.path.join(
        src, f"{args.src_version}_{args.sample}_{SRC_JETS}.parquet"
    )
    n_cache = pq.ParquetFile(jets_file).metadata.num_rows
    print(
        f"check 2: {len(n_ev)} files in the skim tables ({(n_ev < 0).sum()} without "
        f"count), {n_ev.sum()} events; cache {n_cache}"
    )
    if len(n_ev) != len(cache_files) or (n_ev < 0).any() or n_ev.sum() != n_cache:
        fail("per-file event counts do not add up to the cache")
    offset = np.concatenate([[0], np.cumsum(n_ev)])
    rows, norm = None, None
    if args.max_file is not None:
        rows = np.arange(int(n_ev[: args.max_file].sum()))
        print(f"events used: {len(rows)} (files 0-{args.max_file - 1})")
    elif args.events:
        rows, norm = check_events(args, offset, n_ev, src)
        print(f"events used: {len(rows)} (--events {args.events})")

    # 3./4. jets
    raw = read_rows(jets_file, rows)
    max_file = args.max_file
    cand = load_candidates(args.skim_dir, args.skim_sample, args.scores, max_file)
    row_of = {"n_ev": n_ev, "offset": offset, "rows": rows if args.events else None}
    flat = check_and_index(cand, raw, row_of, table)

    dst = os.path.join(constants.CACHE_DIR, args.version)
    os.makedirs(dst, exist_ok=True)
    pattern = os.path.join(dst, f"{args.version}_{out_sample}_%s.parquet")
    for target, arr in build_objects(cand, raw, flat, table).items():
        write(arr, pattern % target)
    if norm is not None:
        norm.update(
            {
                "sample": out_sample,
                "source": f"{args.src_version}_{args.sample}",
                "events": args.events,
                "n_selected": int(len(rows)),
                "norm_events": int(n_cache),
                "rate_per_event_kHz": constants.RATE_NORM_FACTOR / n_cache,
                "complete_for": args.seed_bit,
                "note": "rates are complete only ANDed with a selection inside the "
                "selected events (e.g. the seed bit): purity.norm_events",
            }
        )
        with open(pattern.replace(".parquet", ".json") % "norm", "w") as fh:
            json.dump(norm, fh, indent=2)
        print(norm)

    # other objects: symlinks, or row-aligned subsets
    for obj in OBJECTS[args.sample]:
        source = os.path.join(src, f"{args.src_version}_{args.sample}_{obj}.parquet")
        target = pattern % obj
        if os.path.lexists(target):
            print(f"{target} exists, kept")
        elif rows is None:
            os.symlink(source, target)
            print(f"Linked {target} -> {source}")
        elif obj == "hltAK4PuppiJetCorr":
            # the source has one row group: correct the subset of raw jets instead
            write(correct_jets(raw, table), target)
        else:
            write(read_rows(source, rows), target)


if __name__ == "__main__":
    main()
