"""TEST, not official: paired comparison of the emulated tagger paths (HLT_TEST_*, from
emulate_paths.py) with the HLT ParT and DeepTau double-tau paths on HH -> bbtautau and DY, and
what the two tagged jets are in the events a tagger path takes and the HLT path does not.

Acceptance: >= 2 GenVisTau with pT > 25 GeV and |eta| < 2.1 ("all-acceptance"); also per bin of
the subleading such tau pT. Paired counts against the HLT path: gained (tagger and not HLT), lost
(HLT and not tagger), exact McNemar p-value (two-sided binomial on gained vs gained + lost).

Tagged pair: of the tagger candidates passing the path (|eta| < 2.1, pT > thr, score >= cut, two
with dR > 0.5), the pair with the largest smaller score (the one that decides the bit). Each jet:
  tau     dR < 0.3 to a GenVisTau (any pT);
  lep     else dR < 0.3 to a status-1 e / mu that is prompt or from a prompt tau decay
          (statusFlags bit 0 or 5), pT > 5 GeV;
  b, c, light   else the hadronFlavour of the nearest GenJet within dR < 0.4;
  unmatched     else (pileup or soft).
ID-only efficiency: the path with at least one passing pair whose two jets are both "tau".
ParT and DeepTau cannot be classified this way (no per-jet scores or trigger objects in the nano).

Writes configs/V50nano/object_performance/test_tautag/path_compare_{allacc,bins,jets}.csv.

    python3 scripts/test_tautag/path_diagnosis.py     # 1 core, ~1 min
"""

import os

import awkward as ak
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from menu_tools.utils import constants

VERSION = "V50nano_170pre6HLT_TESTtautag"
SEED = "L1_pDoublePuppiTau52_52_final"
PART = "HLT_DoubleMediumPFPuppiParTTauh30_eta2p1"
DEEP = "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1"
# pT mode of emulated_paths.csv -> (object, pT field)
MODES = {"sigPt": ("hltJetTauTagTEST", "pt"), "ptReg": ("hltJetTauTagRegTEST", "pt"),
         "jetPt": ("hltJetTauTagTEST", "jetPt")}
SCORE = {"BDT": "bdt", "DS": "ds", "DSA": "dsA01"}
ODIR = os.path.join("configs", "V50nano", "object_performance", "test_tautag")
SAMPLES = ["HHTo2B2Tau", "DYLL_M50"]
BINS = {
    "HHTo2B2Tau": [25, 30, 35, 40, 45, 50, 60, 80, 1000],
    "DYLL_M50": [25, 30, 35, 40, 45, 50, 60, 80, 1000],
}
CATS = ["tau", "lep", "b", "c", "light", "unmatched"]


def rec(sample, obj):
    d = os.path.join(constants.CACHE_DIR, VERSION)
    return ak.from_parquet(os.path.join(d, f"{VERSION}_{sample}_{obj}.parquet"))


def dr2(eta1, phi1, eta2, phi2):
    dphi = (phi1 - phi2 + np.pi) % (2 * np.pi) - np.pi
    return (eta1 - eta2) ** 2 + dphi**2


def classify(jets, gen):
    """Category index (CATS) per jet of a jagged array of jets."""
    def near(eta, phi, coll, r):
        pairs = ak.cartesian({"j": ak.zip({"eta": eta, "phi": phi}), "g": coll}, nested=True)
        d = dr2(pairs.j.eta, pairs.j.phi, pairs.g.eta, pairs.g.phi)
        return d, ak.any(d < r * r, axis=-1)

    _, is_tau = near(jets.eta, jets.phi, gen["tau"], 0.3)
    _, is_lep = near(jets.eta, jets.phi, gen["lep"], 0.3)
    d, has_jet = near(jets.eta, jets.phi, gen["jet"], 0.4)
    best = ak.argmin(d, axis=-1, keepdims=True)
    flav = ak.broadcast_arrays(gen["jet"].flav[:, np.newaxis], d)[0]
    flav = ak.fill_none(ak.firsts(flav[best], axis=-1), 0)
    cat = ak.where(is_tau, 0, ak.where(is_lep, 1, ak.where(
        has_jet, ak.where(flav == 5, 2, ak.where(flav == 4, 3, 4)), 5)))
    return cat


def load_gen(sample):
    v = rec(sample, "GenVisTau")
    j = rec(sample, "GenJet")
    src = os.path.join(constants.CACHE_DIR, "V50nano_170pre6HLT",
                       f"V50nano_170pre6HLT_{sample}_GenPart.parquet")
    g = ak.from_parquet(src, columns=["GenPart_pdgId", "GenPart_status", "GenPart_statusFlags",
                                      "GenPart_pt", "GenPart_eta", "GenPart_phi"])
    ap = abs(g.GenPart_pdgId)
    fl = g.GenPart_statusFlags
    lep = ((ap == 11) | (ap == 13)) & (g.GenPart_status == 1) & (g.GenPart_pt > 5) & (
        ((fl & 1) > 0) | ((fl & 32) > 0))
    return {
        "tau": ak.zip({"eta": v.GenVisTau_eta, "phi": v.GenVisTau_phi}),
        "taupt": v.GenVisTau_pt,
        "taustatus": v.GenVisTau_status,
        "lep": ak.zip({"eta": g.GenPart_eta[lep], "phi": g.GenPart_phi[lep]}),
        "jet": ak.zip({"eta": j.GenJet_eta, "phi": j.GenJet_phi, "flav": j.GenJet_hadronFlavour}),
    }


def passing_pairs(arr, prefix, score, cut, thr, ptfield="pt"):
    o = arr[(abs(arr[f"{prefix}_eta"]) < 2.1) & (arr[f"{prefix}_{ptfield}"] > thr)
            & (arr[f"{prefix}_scored"] == 1) & (arr[f"{prefix}_{score}"] >= cut)]
    c = ak.zip({"eta": o[f"{prefix}_eta"], "phi": o[f"{prefix}_phi"], "s": o[f"{prefix}_{score}"]})
    p = ak.combinations(c, 2)
    p = p[dr2(p["0"].eta, p["0"].phi, p["1"].eta, p["1"].phi) > 0.25]
    return p


def mcnemar(g, lost):
    return binomtest(g, g + lost).pvalue if g + lost else 1.0


def main():
    defs = pd.read_csv(os.path.join(ODIR, "emulated_paths.csv"))
    defs = defs[(defs["sample"] == "MinBiasL1Pass") & defs["tagger"].notna()]
    rows_all, rows_bins, rows_jets = [], [], []
    for sample in SAMPLES:
        seed = ak.to_numpy(rec(sample, "L1")[SEED]).astype(bool)
        hlt = rec(sample, "HLT")
        ref = {"ParT": ak.to_numpy(hlt[PART]).astype(bool),
               "DeepTau": ak.to_numpy(hlt[DEEP]).astype(bool)}
        test = rec(sample, "HLTTEST")
        gen = load_gen(sample)
        acc_t = (gen["taupt"] > 25) & (abs(gen["tau"].eta) < 2.1)
        nacc = ak.to_numpy(ak.sum(acc_t, axis=1))
        acc = nacc >= 2
        subl = ak.to_numpy(ak.fill_none(ak.pad_none(ak.sort(gen["taupt"][acc_t], ascending=False),
                                                    2)[:, 1], 0.0))
        objs = {o: rec(sample, o) for o in {m[0] for m in MODES.values()}}
        n = int(acc.sum())
        for r in ("ParT", "DeepTau"):
            rows_all.append(dict(sample=sample, path=r, ref=None, n_acc=n,
                                 n_pass=int((ref[r] & acc).sum()),
                                 eff=(ref[r] & acc).sum() / n))
        for _, d in defs.iterrows():
            bit = d["bit"]
            r = "ParT" if "ParT" in bit else "DeepTau"
            t = ak.to_numpy(test[bit]).astype(bool)
            # tagged pair and its classification
            obj, ptf = MODES[d["pt"]]
            pp = passing_pairs(objs[obj], obj, SCORE[d["tagger"]], d["cut"], d["threshold"], ptf)
            p0 = classify(pp["0"], gen)
            p1 = classify(pp["1"], gen)
            both_tau = ak.to_numpy(ak.any((p0 == 0) & (p1 == 0), axis=1)) & seed
            assert np.array_equal(ak.to_numpy(ak.num(pp) > 0) & seed, t), bit
            mins = ak.where(pp["0"].s < pp["1"].s, pp["0"].s, pp["1"].s)
            best = ak.argmax(mins, axis=1, keepdims=True)
            c0 = ak.to_numpy(ak.fill_none(ak.firsts(p0[best], axis=1), -1))
            c1 = ak.to_numpy(ak.fill_none(ak.firsts(p1[best], axis=1), -1))
            g = t & ~ref[r]
            lost = ref[r] & ~t
            rows_all.append(dict(
                sample=sample, path=bit, ref=r, n_acc=n, n_pass=int((t & acc).sum()),
                eff=(t & acc).sum() / n, n_ref=int((ref[r] & acc).sum()),
                gained=int((g & acc).sum()), lost=int((lost & acc).sum()),
                p_mcnemar=mcnemar(int((g & acc).sum()), int((lost & acc).sum())),
                eff_idonly=(both_tau & acc).sum() / n,
            ))
            e = BINS[sample]
            for lo, hi in zip(e[:-1], e[1:]):
                m = acc & (subl >= lo) & (subl < hi)
                rows_bins.append(dict(
                    sample=sample, path=bit, ref=r, lo=lo, hi=hi, n=int(m.sum()),
                    eff=(t & m).sum() / max(m.sum(), 1), eff_ref=(ref[r] & m).sum() / max(m.sum(), 1),
                    gained=int((g & m).sum()), lost=int((lost & m).sum()),
                    p_mcnemar=mcnemar(int((g & m).sum()), int((lost & m).sum())),
                    eff_idonly=(both_tau & m).sum() / max(m.sum(), 1)))
            # jet composition: all events of the tagger path, and gained events, in acceptance
            # and in all events
            for sel_name, sel in [("fired_acc", t & acc), ("gained_acc", g & acc),
                                  ("gained_all", g), ("fired_all", t)]:
                cats = np.concatenate([c0[sel], c1[sel]])
                row = dict(sample=sample, path=bit, ref=r, selection=sel_name,
                           n_events=int(sel.sum()))
                for i, cn in enumerate(CATS):
                    row[cn] = (cats == i).mean() if len(cats) else np.nan
                row["pair_tautau"] = ((c0[sel] == 0) & (c1[sel] == 0)).mean() if sel.any() else np.nan
                row["pair_with_b"] = ((c0[sel] == 2) | (c1[sel] == 2)).mean() if sel.any() else np.nan
                row["pair_with_lep"] = ((c0[sel] == 1) | (c1[sel] == 1)).mean() if sel.any() else np.nan
                rows_jets.append(row)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    for name, rows in [("allacc", rows_all), ("bins", rows_bins), ("jets", rows_jets)]:
        df = pd.DataFrame(rows)
        print(f"== {name}")
        print(df.to_string(float_format="%.3f"))
        df.to_csv(os.path.join(ODIR, f"path_compare_{name}.csv"), index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
