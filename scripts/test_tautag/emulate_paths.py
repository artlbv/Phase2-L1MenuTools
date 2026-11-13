"""TEST, not official: emulated double-tau "paths" of the TEST tau taggers, as per-event
bits for the path-efficiency plots (configs/V50nano/object_performance/test_tautag/
test_tautag_path_eff.yaml), next to the real HLT path bits.

Emulated path: L1_pDoublePuppiTau52_52_final AND >= 2 tagger candidates (hltJetTauTagTEST,
pt = pT^sig, or hltJetTauTagRegTEST, pt = the DeepSets regressed pT) with |eta| < 2.1,
pt > threshold and score > cut, two of them with dR > 0.5 (the double-tau emulation of the
hlt-tau-pu report and of scripts/test_tautag/wps.py; the HLT double-tau paths also require
dR > 0.5 between the legs).

Bits (HLT_TEST_Double<tagger>[Reg]_<WP>_<thr>, Reg: threshold on the regressed pT):
  *_ParT_30   the score cut chosen at the rate of HLT_DoubleMediumPFPuppiParTTauh30_eta2p1:
              at most as many minimum-bias events as that path (11 of the 1,999,040) on
              MinBiasL1Pass (which holds every event of the 2M with the L1 seed), at 30 GeV
              like the ParT path;
  *_1kHz_35   the 1 kHz WP of configs/V50nano/objects/test_tautag.yaml at 35 GeV
              (DeepTau path: 66 events);
  DoubleDSReg_1kHz_26   ds > 0.607 with regressed pT > 26 GeV (1 kHz, the HH efficiency of
              pT^sig > 35).
Writes <version>_<sample>_HLTTEST.parquet (row-aligned with the other caches; read by
path_efficiency.py for HLT_TEST_* bits) for HHTo2B2Tau, MinBiasF1000 and MinBiasL1Pass, and
configs/V50nano/object_performance/test_tautag/emulated_paths.csv (cuts, minimum-bias counts and
rates with 68% Garwood intervals, 1 event = 31,038 kHz / 1,999,040 = 15.5 Hz).

    python3 scripts/test_tautag/emulate_paths.py      # 1 core, < 1 min
"""

import os

import awkward as ak
import numpy as np
import pandas as pd
from scipy.stats import chi2

from menu_tools.utils import constants

VERSION = "V50nano_170pre6HLT_TESTtautag"
SEED = "L1_pDoublePuppiTau52_52_final"
PATHS = {
    "DeepTau": "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1",
    "ParT": "HLT_DoubleMediumPFPuppiParTTauh30_eta2p1",
}
N_TOTAL = 1999040  # MinBiasL1Pass normalisation (the full MinBias cache)
OBJ = {False: "hltJetTauTagTEST", True: "hltJetTauTagRegTEST"}
TAGGERS = {"BDT": "bdt", "DS": "ds", "DSA": "dsA01"}
# fixed cuts: name -> (tagger, regressed pT, cut, threshold)
FIXED = {
    "HLT_TEST_DoubleBDT_1kHz_35": ("BDT", False, 0.363, 35.0),
    "HLT_TEST_DoubleDS_1kHz_35": ("DS", False, 0.460, 35.0),
    "HLT_TEST_DoubleDSA_1kHz_35": ("DSA", False, 0.398, 35.0),
    "HLT_TEST_DoubleDSReg_1kHz_26": ("DS", True, 0.607, 26.0),
}
ODIR = os.path.join("configs", "V50nano", "object_performance", "test_tautag")


def garwood(n, cl=0.6827):
    a = (1 - cl) / 2
    lo = chi2.ppf(a, 2 * n) / 2 if n > 0 else 0.0
    return lo, chi2.ppf(1 - a, 2 * (n + 1)) / 2


def pair_score(arr, prefix, score, thr):
    """Per event: the largest s such that two candidates with |eta| < 2.1, pt > thr,
    dR > 0.5 both have score >= s (-inf without such a pair)."""
    o = arr[(abs(arr[f"{prefix}_eta"]) < 2.1) & (arr[f"{prefix}_pt"] > thr)
            & (arr[f"{prefix}_scored"] == 1)]
    c = ak.zip({"eta": o[f"{prefix}_eta"], "phi": o[f"{prefix}_phi"], "s": o[f"{prefix}_{score}"]})
    p = ak.combinations(c, 2)
    dphi = (p["0"].phi - p["1"].phi + np.pi) % (2 * np.pi) - np.pi
    ok = (p["0"].eta - p["1"].eta) ** 2 + dphi**2 > 0.25
    s = ak.where(p["0"].s < p["1"].s, p["0"].s, p["1"].s)[ok]
    return ak.to_numpy(ak.fill_none(ak.max(s, axis=1), -np.inf)).astype(float)


def load(sample):
    d = os.path.join(constants.CACHE_DIR, VERSION)
    f = lambda o: os.path.join(d, f"{VERSION}_{sample}_{o}.parquet")  # noqa: E731
    seed = ak.to_numpy(ak.from_parquet(f("L1"), columns=[SEED])[SEED]).astype(bool)
    hlt = ak.from_parquet(f("HLT"), columns=list(PATHS.values()))
    arrs = {r: ak.from_parquet(f(o)) for r, o in OBJ.items()}
    return f, seed, hlt, arrs


def main():
    # cuts at the ParT-path rate on MinBiasL1Pass
    _, seed, hlt, arrs = load("MinBiasL1Pass")
    n_part = int(ak.sum(hlt[PATHS["ParT"]]))
    defs = dict(FIXED)
    for t, score in TAGGERS.items():
        for reg in (False, True):
            ps = np.where(seed, pair_score(arrs[reg], OBJ[reg], score, 30.0), -np.inf)
            s = np.sort(ps[np.isfinite(ps)])[::-1]
            cut = float(np.nextafter(s[n_part], np.inf))  # at most n_part events >= cut
            name = f"HLT_TEST_Double{t}{'Reg' if reg else ''}_ParT_30"
            defs[name] = (t, reg, cut, 30.0)
    rows = []
    for sample in ("MinBiasL1Pass", "MinBiasF1000", "HHTo2B2Tau"):
        f, seed, hlt, arrs = load(sample)
        out, cache = {}, {}
        for name, (t, reg, cut, thr) in defs.items():
            key = (reg, TAGGERS[t], thr)
            if key not in cache:
                cache[key] = pair_score(arrs[reg], OBJ[reg], TAGGERS[t], thr)
            out[name] = seed & (cache[key] >= cut)
        ak.to_parquet(ak.zip(out), f("HLTTEST") + ".tmp")
        os.replace(f("HLTTEST") + ".tmp", f("HLTTEST"))
        n_tot = N_TOTAL if sample == "MinBiasL1Pass" else len(seed)
        scale = constants.RATE_NORM_FACTOR / n_tot
        bits = {**{v: ak.to_numpy(hlt[v]).astype(bool) for v in PATHS.values()},
                SEED: seed, **out}
        for name, m in bits.items():
            n = int(m.sum())
            lo, hi = garwood(n)
            t, reg, cut, thr = defs.get(name, (None, None, np.nan, np.nan))
            rows.append(dict(sample=sample, bit=name, tagger=t, pt="ptReg" if reg else
                             ("sigPt" if t else None), cut=cut, threshold=thr, n=n,
                             frac=n / len(m), rate_khz=n * scale, rate_lo=lo * scale,
                             rate_hi=hi * scale))
    df = pd.DataFrame(rows)
    df.loc[df["sample"] == "HHTo2B2Tau", ["rate_khz", "rate_lo", "rate_hi"]] = np.nan
    pd.set_option("display.width", 250)
    print(df.to_string(float_format="%.4f"))
    df.to_csv(os.path.join(ODIR, "emulated_paths.csv"), index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
