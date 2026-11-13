"""TEST, not official: working points of the TEST tau taggers (hltJetTauTagTEST,
configs/V50nano/objects/test_tautag.yaml) with the machinery of the hlt-tau-pu study
(topo-hackathon-2026 projects/hlt-tau-pu/scripts, results/hlt_tau_pu/README.md).

  *1kHz   the cut of the report's double-tau emulation at the rate of
          HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1: L1_pDoublePuppiTau52_52 AND two
          candidates with pT^sig > 35 GeV, |eta| < 2.1, dR > 0.5, score >= cut; at most
          int(1 kHz / 31,038 kHz * 1,999,040) = 64 minimum-bias events (all 2M events:
          the skim holds every seeded event of every file).
  *Eff90 / *Eff80   90% / 80% per-candidate efficiency on genuine HH tau_h (skim label
          0) with pT^sig > 20 GeV.
Regressed pT (DeepSets ptReg): with the DeepSets 1 kHz cut fixed, the ptReg threshold
that gives the same minimum-bias count as pT^sig > 35 ("matched rate"), and the 1 kHz cut
re-chosen with ptReg > 35. HH efficiency: events with 2 visible tau_h > 25 GeV (report).

    python3 scripts/test_tautag/wps.py      # 1 core, ~2 min; prints and writes
        configs/V50nano/object_performance/test_tautag/wps.csv
"""

import os
import sys

import numpy as np
import pandas as pd

HTP = "/eos/home-a/alobanov/SWAN_projects/L1TML/topo-hackathon-2026/projects/hlt-tau-pu/scripts"
sys.path.insert(0, HTP)
import evaluate as ev  # noqa: E402
import evaluate_deepsets as ED  # noqa: E402
from common import OUT, MINBIAS_RATE_KHZ, DITAU_L1, garwood, load  # noqa: E402

KEY = ["file", "entry", "jet_idx"]
DS = "ds_mcttreg_jv_s0"
ODIR = os.path.join("configs", "V50nano", "object_performance", "test_tautag")


def main():
    cols = KEY + ["label", "sig_pt", "axis_eta", "axis_phi"]
    parts = []
    for s in ("hh", "minbias"):
        c = load(s, "cand", "jet", columns=cols)
        c["sample"] = s
        parts.append(c)
    c = pd.concat(parts, ignore_index=True)
    for f in ("scored_jet.parquet", f"deepsets/scores_{DS[3:]}.parquet",
              "deepsets/scores_mc_jv.parquet"):
        sc = pd.read_parquet(os.path.join(OUT, f))
        sc = sc[sc["sample"].isin(["hh", "minbias"])]
        sc = sc.drop(columns=[k for k in ("parity", "label") if k in sc])
        c = c.merge(sc, on=["sample"] + KEY, how="left", validate="1:1")
    eps = 1e-30
    c["ds"] = c[DS] / np.maximum(c[DS] + c[f"{DS}_pmb"], eps)
    c["dsA01"] = c[DS] / np.maximum(c[DS] + c[f"{DS}_pmb"] + 0.1 * c[f"{DS}_phs"], eps)
    c["dsOld"] = c.ds_mc_jv / np.maximum(c.ds_mc_jv + c.ds_mc_jv_pmb, eps)
    c["bdt"] = c.bdt_full
    c["ptReg"] = c[f"{DS}_ptcorr"]
    tags = ["bdt", "ds", "dsA01", "dsOld"]
    assert c[tags + ["ptReg"]].notna().all().all()

    ev_hh, ev_mb = load("hh", "ev", "jet"), load("minbias", "ev", "jet")
    files_mb = load("minbias", "files", "jet")
    n_mb = int(files_mb.loc[files_mb.n_events > 0, "n_events"].sum())
    nn = int(1.0 / MINBIAS_RATE_KHZ * n_mb)
    acc = ev_hh[ev_hh.n_gvt20 >= 2]
    a25 = acc.n_gvt25.values >= 2
    print(f"minimum bias {n_mb} events, 1 kHz = {nn} events; seed fires in "
          f"{int(ev_mb[DITAU_L1].sum())}; HH accepted events {a25.sum()}")
    hh, mb = c[c["sample"] == "hh"], c[c["sample"] == "minbias"]
    sig = hh[(hh.label == 0) & (hh.sig_pt > 20)]
    rows = []

    def ev_scores(cand, evt, ptcol, ptmin, taggers):
        return __import__("evaluate_reg").event_scores_pt(cand, evt, taggers, ptcol, ptmin)

    # 1 kHz cuts with pT^sig > 35
    S_hh = ev_scores(hh, acc, "sig_pt", 35.0, tags)
    S_mb = ev_scores(mb, ev_mb, "sig_pt", 35.0, tags)
    for t in tags:
        thr = ED.thr_for(S_mb[t].values, nn)
        k = int((S_mb[t].values >= thr).sum())
        eff = (S_hh[t].values >= thr)[a25].mean()
        lo, hi = garwood(k)
        rows.append(dict(score=t, wp="1kHz", pt="sigPt", pt_thr=35.0, cut=thr, n_mb=k,
                         rate_khz=k * MINBIAS_RATE_KHZ / n_mb,
                         rate_lo=lo * MINBIAS_RATE_KHZ / n_mb, rate_hi=hi * MINBIAS_RATE_KHZ / n_mb,
                         hh_eff_acc25=eff, tau_eff20=(sig[t] >= thr).mean()))
        for e in (0.9, 0.8):
            cut = float(np.quantile(sig[t].values, 1 - e))
            k = int((S_mb[t].values >= cut).sum())
            rows.append(dict(score=t, wp=f"Eff{int(e * 100)}", pt="sigPt", pt_thr=35.0, cut=cut,
                             n_mb=k, rate_khz=k * MINBIAS_RATE_KHZ / n_mb,
                             hh_eff_acc25=(S_hh[t].values >= cut)[a25].mean(),
                             tau_eff20=(sig[t] > cut).mean(),
                             mb_fake20=(mb[(mb.file < 1000) & (mb.sig_pt > 20)][t] > cut).mean()))
    wp = pd.DataFrame(rows)
    cut1k = wp[(wp.score == "ds") & (wp.wp == "1kHz")].cut.iloc[0]
    n35 = wp[(wp.score == "ds") & (wp.wp == "1kHz")].n_mb.iloc[0]

    # regressed pT: threshold scan with the DeepSets 1 kHz cut fixed
    scan = []
    for x in np.arange(20.0, 40.01, 0.5):
        Sh = ev_scores(hh, acc, "ptReg", x, ["ds"])
        Sm = ev_scores(mb, ev_mb, "ptReg", x, ["ds"])
        k = int((Sm.ds.values >= cut1k).sum())
        c1 = ED.thr_for(Sm.ds.values, nn)
        scan.append(dict(pt_thr=x, n_mb=k, rate_khz=k * MINBIAS_RATE_KHZ / n_mb,
                         hh_eff_acc25=(Sh.ds.values >= cut1k)[a25].mean(), cut_1khz=c1,
                         n_mb_1khz=int((Sm.ds.values >= c1).sum()),
                         hh_eff_acc25_1khz=(Sh.ds.values >= c1)[a25].mean()))
    scan = pd.DataFrame(scan)
    ok = scan[scan.n_mb <= n35]
    xm = float(ok.pt_thr.min())
    r = scan[scan.pt_thr == xm].iloc[0]
    rows = [dict(score="ds", wp="1kHz", pt="ptReg", pt_thr=xm, cut=cut1k, n_mb=int(r.n_mb),
                 rate_khz=r.rate_khz, hh_eff_acc25=r.hh_eff_acc25, note="matched rate, DS1kHz cut")]
    r35 = scan[scan.pt_thr == 35.0].iloc[0]
    rows.append(dict(score="ds", wp="1kHz", pt="ptReg", pt_thr=35.0, cut=cut1k, n_mb=int(r35.n_mb),
                     rate_khz=r35.rate_khz, hh_eff_acc25=r35.hh_eff_acc25, note="DS1kHz cut"))
    # 1 kHz cut re-chosen per ptReg threshold: at 35 GeV, and at the highest threshold
    # whose HH efficiency is at least that of pT^sig > 35 (equal efficiency at 1 kHz)
    e35 = wp[(wp.score == "ds") & (wp.wp == "1kHz")].hh_eff_acc25.iloc[0]
    xe = float(scan[scan.hh_eff_acc25_1khz >= e35].pt_thr.max())
    for x, note in ((35.0, "cut re-chosen for 1 kHz"),
                    (xe, "cut re-chosen for 1 kHz, equal HH eff. to pT^sig > 35")):
        r = scan[scan.pt_thr == x].iloc[0]
        rows.append(dict(score="ds", wp="1kHz", pt="ptReg", pt_thr=x, cut=r.cut_1khz,
                         n_mb=int(r.n_mb_1khz), rate_khz=r.n_mb_1khz * MINBIAS_RATE_KHZ / n_mb,
                         hh_eff_acc25=r.hh_eff_acc25_1khz, note=note))
    wp = pd.concat([wp, pd.DataFrame(rows)], ignore_index=True)
    pd.set_option("display.width", 250)
    print(wp.to_string(float_format="%.4f"))
    print(scan.to_string(float_format="%.4f"))
    wp.to_csv(os.path.join(ODIR, "wps.csv"), index=False, float_format="%.4f")
    scan.to_csv(os.path.join(ODIR, "wps_ptreg_scan.csv"), index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
