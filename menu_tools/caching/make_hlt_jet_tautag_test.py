"""TEST, not official: the experimental HLT-jet tau taggers as a cached object.

Writes the object `hltJetTauTagTEST`, i.e. the HLT AK4 PUPPI jets with the out-of-fold
scores of the tau taggers of the hlt-tau-pu study (topo-hackathon-2026,
results/hlt_tau_pu/README.md): the BDT (`bdt_full`) and the 3-class DeepSets with jet
variables (`mc_jv`). The taggers are not HLT objects of any menu; this is a TEST.

Everything goes into a separate cache version (default V50nano_170pre6HLT_TESTtautag,
so the plots land in outputs/<that version>/ and never mix with the official ones):
  <version>_<sample>_hltJetTauTagTEST.parquet   the tagger object, one entry per
      cached hltAK4PuppiJet (raw nano order -> sorted by `pt`). Fields:
        pt, eta, phi   tau candidate: pT^sig (|sum p| of the PF candidates within
                       dR < 0.1 of the axis, the tau energy estimate the report's
                       trigger cuts on) and the axis (leading charged PF candidate)
                       for scored jets; pt = -1 and the jet eta / phi otherwise
        jetPt, jetEta, jetPhi, rawPt   corrected / raw jet
        bdt            BDT score (bdt_full)
        ds             DeepSets trigger score p_tau / (p_tau + p_MB)
        dsPtau, dsPmb, dsPhs           DeepSets class probabilities
        scored         1 for scored jets (corrected pT > 20 GeV, |eta| < 2.5)
        label          skim label (0 tau, 1 hard-scatter jet, 2 pileup / MB, 3 e/mu)
      Unscored jets have all scores and the label at -1.
  <version>_<sample>_<object>.parquet for the other objects the plots need: symlinks
      to the source version, or, with --max-file, the row-aligned subset of the
      events of files [0, max-file) (e.g. sample MinBiasF1000: the minimum-bias files
      in which every event was scored).

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

Usage (1 core, a few minutes per sample):
    python3 menu_tools/caching/make_hlt_jet_tautag_test.py --sample HHTo2B2Tau \\
        --skim-sample hh
    python3 menu_tools/caching/make_hlt_jet_tautag_test.py --sample MinBias \\
        --skim-sample minbias --max-file 1000 --out-sample MinBiasF1000
"""

import argparse
import glob
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
SRC_JETS = "hltAK4PuppiJet"
HTAUPU = (
    "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/alobanov/phase2/menu/"
    "ntuples/hlt_tau_pu"
)
CACHING_CFG = os.path.join(
    "configs", "V50nano", "cache_objects", "caching_170pre6_HLT.yaml"
)
# Objects the TEST plots read besides the tagger object
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


def read_rows(fname: str, n_rows: Optional[int]) -> ak.Array:
    """The first n_rows events of a cache parquet (all if None), reading only the row
    groups needed."""
    if n_rows is None:
        return ak.from_parquet(fname)
    md = pq.ParquetFile(fname).metadata
    groups, tot = [], 0
    for i in range(md.num_row_groups):
        if tot >= n_rows:
            break
        groups.append(i)
        tot += md.row_group(i).num_rows
    return ak.from_parquet(fname, row_groups=groups)[:n_rows]


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
        s = s[s["sample"] == skim_sample].drop(
            columns=[c for c in ("sample", "parity", "label") if c in s]
        )
        if max_file is not None:
            s = s[s.file < max_file]
        if s.duplicated(KEY).any():
            fail(f"duplicate keys in {path}")
        n0 = len(cand)
        cand = cand.merge(s, on=KEY, how="inner", validate="one_to_one")
        if len(cand) != n0 or len(s) != n0:
            fail(f"{path}: {len(s)} scores for {n0} candidates, {len(cand)} matched")
    return cand


def check_and_index(
    cand: pd.DataFrame, raw: ak.Array, n_ev: np.ndarray, table: np.ndarray
) -> np.ndarray:
    """Checks 3 and 4; returns the flat index of each candidate in the cached jets."""
    offset = np.concatenate([[0], np.cumsum(n_ev)])
    f = cand.file.to_numpy()
    e = cand.entry.to_numpy()
    if (e >= n_ev[f]).any() or (e < 0).any():
        fail("entry outside the file's event range")
    row = offset[f] + e
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


def build_object(cand: pd.DataFrame, raw: ak.Array, flat: np.ndarray, table):
    counts = ak.num(raw[f"{SRC_JETS}_pt"])
    pt = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_pt"])).astype(np.float64)
    eta = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_eta"]))
    phi = ak.to_numpy(ak.flatten(raw[f"{SRC_JETS}_phi"]))
    n = len(pt)

    def full(values, default):
        out = np.full(n, default, dtype=np.float32)
        out[flat] = values
        return out

    p_tau = cand.ds_mc_jv.to_numpy()
    p_mb = cand.ds_mc_jv_pmb.to_numpy()
    ds = np.where(p_tau + p_mb > 0, p_tau / np.maximum(p_tau + p_mb, 1e-30), 0.0)
    fields = {
        "pt": full(cand.sig_pt.to_numpy(), -1),
        "eta": eta.astype(np.float32),
        "phi": phi.astype(np.float32),
        "jetPt": (pt * l2relative(pt, eta, table)).astype(np.float32),
        "jetEta": eta.astype(np.float32),
        "jetPhi": phi.astype(np.float32),
        "rawPt": pt.astype(np.float32),
        "bdt": full(cand.bdt_full.to_numpy(), -1),
        "ds": full(ds, -1),
        "dsPtau": full(p_tau, -1),
        "dsPmb": full(p_mb, -1),
        "dsPhs": full(cand.ds_mc_jv_phs.to_numpy(), -1),
        "scored": full(np.ones(len(cand)), 0).astype(np.int32),
        "label": full(cand.label.to_numpy(), -1).astype(np.int32),
    }
    fields["eta"][flat] = cand.axis_eta.to_numpy()
    fields["phi"][flat] = cand.axis_phi.to_numpy()
    arr = ak.unflatten(ak.zip({f"{TARGET}_{k}": v for k, v in fields.items()}), counts)
    return arr[ak.argsort(arr[f"{TARGET}_pt"], axis=-1, ascending=False)]


def write(arr: ak.Array, target: str) -> None:
    ak.to_parquet(arr, target + ".tmp")
    os.replace(target + ".tmp", target)
    print(f"Wrote {target}: {len(arr)} events")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--sample", required=True, choices=list(OBJECTS))
    p.add_argument("--skim-sample", required=True, help="hh, dy, minbias")
    p.add_argument("--out-sample", default=None, help="default: --sample")
    p.add_argument("--max-file", type=int, default=None)
    p.add_argument("--src-version", default="V50nano_170pre6HLT")
    p.add_argument("--version", default="V50nano_170pre6HLT_TESTtautag")
    p.add_argument("--skim-dir", default=os.path.join(HTAUPU, "skim_jet"))
    p.add_argument("--file-list", default=None, help="default: skim/files_<s>.txt")
    p.add_argument(
        "--scores",
        nargs="+",
        default=[
            os.path.join(HTAUPU, "scored_jet.parquet"),
            os.path.join(HTAUPU, "deepsets", "scores_mc_jv.parquet"),
        ],
    )
    p.add_argument("--caching-cfg", default=CACHING_CFG)
    p.add_argument("--jec", default=os.path.join(JEC_DIR, JEC_L2))
    args = p.parse_args()
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

    # 2. events per file
    n_ev = file_events(args.skim_dir, args.skim_sample)
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
    if args.max_file is not None:
        n_ev = n_ev[: args.max_file]
    n_rows = int(n_ev.sum()) if args.max_file is not None else None
    print(f"events used: {n_ev.sum()} (files 0-{len(n_ev) - 1})")

    # 3./4. jets
    raw = read_rows(jets_file, n_rows)
    cand = load_candidates(args.skim_dir, args.skim_sample, args.scores, args.max_file)
    flat = check_and_index(cand, raw, n_ev, table)

    dst = os.path.join(constants.CACHE_DIR, args.version)
    os.makedirs(dst, exist_ok=True)
    pattern = os.path.join(dst, f"{args.version}_{out_sample}_%s.parquet")
    write(build_object(cand, raw, flat, table), pattern % TARGET)

    # other objects: symlinks, or row-aligned subsets
    for obj in OBJECTS[args.sample]:
        source = os.path.join(src, f"{args.src_version}_{args.sample}_{obj}.parquet")
        target = pattern % obj
        if os.path.lexists(target):
            print(f"{target} exists, kept")
        elif n_rows is None:
            os.symlink(source, target)
            print(f"Linked {target} -> {source}")
        elif obj == "hltAK4PuppiJetCorr":
            # the source has one row group: correct the subset of raw jets instead
            write(correct_jets(raw, table), target)
        else:
            write(read_rows(source, n_rows), target)


if __name__ == "__main__":
    main()
