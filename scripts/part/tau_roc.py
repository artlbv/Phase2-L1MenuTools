"""Per-jet tau-ID ROC curves on a common denominator of HLT AK4 PUPPI jets (p2nano v2 + TEST taggers).

Denominator: hltAK4PuppiJetCorr jets with JEC-corrected pT > 20 GeV and |eta| < 2.1, from the cache
version V50nano_170pre6HLTv2_TESTtautag (TEST taggers and v2 ParT row-aligned with the v1 events).
Each jet gets one score per tagger, taken from that tagger's own object:
  part   hltAK4PuppiJetParT TauvsAll (same jets; nearest within dR < 0.05)
  ds, dsA01, bdt   hltJetTauTagTEST (same jets; jetEta / jetPhi nearest within dR < 0.05)
  hps    hltHpsPFTau nearest within dR < 0.3: deepTauVSjet - Medium(pT_HPS), so 0 = HLT Medium WP
         (HPS taus are stored only in events with an L1 tau or e/gamma seed)
  l1     L1nnPuppiTau nearest within dR < 0.3 (any pT): 1 matched, 2 passLooseNN, 3 passTightNN
         (the nano has no continuous NN score, only the flags)
A missing score is the lowest value, so the ROC curves end at the matching efficiency.

Labels: signal = jet within dR < 0.3 of a GenVisTau with pT > 20 GeV (HH and DY separately).
Backgrounds: MinBiasF1000 jets (all; minimum-bias files 0-999, every event); HH b jets and HH
light/gluon jets (nearest GenJet with pT > 15 GeV within dR < 0.4, hadronFlavour 5 / 0), not within
dR < 0.3 of any GenVisTau or of a prompt e / mu (GenPart pT > 10, isPrompt or
isDirectPromptTauDecayProduct, isLastCopy; GenPart from the v2 cache, mapped to the TEST rows by the
GenJet content of the event).

    python3 scripts/part/tau_roc.py [--outdir <www dir>] [--csv <csv>] [--max-events N]
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402

import awkward as ak  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CACHE = os.environ.get(
    "MENU_TOOLS_CACHE",
    "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/alobanov/phase2/menu/ntuples/cache_ar26",
)
VERSION = "V50nano_170pre6HLTv2_TESTtautag"
V2 = "V50nano_170pre6HLTv2"
WWW = (
    "/eos/user/a/alobanov/www/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag_v2/"
    "object_performance/roc"
)
CSV = "configs/V50nano/object_performance/part/tau_roc.csv"

MISSING = -1e9
TAGGERS = {  # key: (label, color, linestyle)
    "part": ("HLT ParT TauvsAll", "#e42536", "-"),
    "ds": (r"TEST DeepSets $p_\tau/(p_\tau+p_{MB})$", "#5790fc", "-"),
    "dsA01": (r"TEST DeepSets $\alpha$ = 0.1", "#964a8b", "-"),
    "bdt": ("TEST BDT", "#f89c20", "-"),
    "hps": ("HLT HPS DeepTau (- Medium)", "#7a21dd", "--"),
    "l1": ("L1 nnPuppiTau (match / LooseNN / TightNN)", "#9c9ca1", ":"),
}
WPS = {  # tagger: [(name, cut)]
    "part": [("MediumParT", 0.11645), ("ParT1kHz", 0.01834)],
    "ds": [("DS1kHz", 0.460), ("DSParT", 0.919)],
    "dsA01": [("DSA1kHz", 0.398)],
    "bdt": [("BDT1kHz", 0.363)],
    "hps": [("Medium DeepTau", 0.0)],
    "l1": [("LooseNN", 2.0), ("TightNN", 3.0)],
}
WP_MARKERS = ["o", "s", "^"]
ETA = {"barrel": (0.0, 1.5), "endcap": (1.5, 2.1), "all": (0.0, 2.1)}
PT = {"pt20": (20.0, 1e9), "pt20to40": (20.0, 40.0), "pt40to100": (40.0, 100.0)}
PT_LABEL = {"pt20": r"$p_T$ > 20 GeV", "pt20to40": r"20 < $p_T$ < 40 GeV", "pt40to100": r"40 < $p_T$ < 100 GeV"}
ETA_LABEL = {"barrel": r"|$\eta$| < 1.5", "endcap": r"1.5 < |$\eta$| < 2.1", "all": r"|$\eta$| < 2.1"}
BKG = {
    "MB": ("MinBiasF1000", "mb", "minimum-bias jets"),
    "HHb": ("HHTo2B2Tau", "b", "HH b jets"),
    "HHlight": ("HHTo2B2Tau", "light", "HH light/g jets"),
}
SIG = {"HH": "HHTo2B2Tau", "DY": "DYLL_M50"}
TARGETS = [0.01, 0.05, 0.10]


def path(sample, obj, version=VERSION):
    return os.path.join(CACHE, version, f"{version}_{sample}_{obj}.parquet")


def read(sample, obj, cols, n=None, version=VERSION):
    a = ak.from_parquet(path(sample, obj, version), columns=[f"{obj}_{c}" for c in cols])
    if n is not None:
        a = a[:n]
    return {c: a[f"{obj}_{c}"] for c in cols}


def flat_event(jets_eta):
    """event index of every flattened jet."""
    n = ak.to_numpy(ak.num(jets_eta))
    return np.repeat(np.arange(len(n)), n), n


def nearest(jeta, jphi, oeta, ophi, dr_max):
    """For each jet (jagged), flat index into the flattened other collection of the nearest
    object within dr_max, else -1. Returns a flat numpy array aligned with ak.flatten(jeta)."""
    pairs_eta = ak.cartesian([jeta, oeta], nested=True)
    pairs_phi = ak.cartesian([jphi, ophi], nested=True)
    deta = pairs_eta["0"] - pairs_eta["1"]
    dphi = (pairs_phi["0"] - pairs_phi["1"] + np.pi) % (2 * np.pi) - np.pi
    dr2 = deta**2 + dphi**2
    idx = ak.argmin(dr2, axis=2)
    mind = ak.min(dr2, axis=2)
    idx = ak.to_numpy(ak.fill_none(ak.flatten(idx), -1)).astype(np.int64)
    mind = ak.to_numpy(ak.fill_none(ak.flatten(mind), np.inf))
    evt, _ = flat_event(jeta)
    off = np.concatenate([[0], np.cumsum(ak.to_numpy(ak.num(oeta)))[:-1]]).astype(np.int64)
    ok = (idx >= 0) & (mind < dr_max**2)
    return np.where(ok, off[evt] + idx, -1)


def take(flat_vals, idx, missing=MISSING):
    out = np.full(len(idx), missing, dtype=np.float64)
    ok = idx >= 0
    out[ok] = flat_vals[idx[ok]]
    return out


def hps_medium(pt):
    """HLT Medium DeepTau WP (hltHpsSelectedPFTausMediumDitauWPDeepTau), piecewise linear in pT."""
    return np.interp(pt, [35.0, 100.0, 300.0], [0.649, 0.441, 0.05])


def ev_key(genjet_pt):
    """hash of the GenJet pT list per event (gen level is identical in v1 and v2)."""
    n = ak.to_numpy(ak.num(genjet_pt))
    flat = ak.to_numpy(ak.flatten(genjet_pt)).astype(np.float32)
    off = np.concatenate([[0], np.cumsum(n)])
    return np.array(
        [hashlib.blake2b(flat[off[i] : off[i + 1]].tobytes(), digest_size=8).hexdigest() for i in range(len(n))]
    )


def lepton_veto(jeta, jphi, n_events):
    """True for jets within dR < 0.3 of a prompt e / mu (HH; GenPart from the v2 cache)."""
    gj1 = read("HHTo2B2Tau", "GenJet", ["pt"], n_events)["pt"]
    gj2 = read("HHTo2B2Tau", "GenJet", ["pt"], None, V2)["pt"]
    k1, k2 = ev_key(gj1), ev_key(gj2)
    s2 = pd.Series(np.arange(len(k2)), index=k2)
    s2 = s2[~s2.index.duplicated(keep=False)]
    row2 = s2.reindex(k1).to_numpy()
    found = ~np.isnan(row2)
    print(f"  HH GenPart mapping v1->v2 by GenJet content: {found.sum()} / {len(k1)} events")
    gp = read("HHTo2B2Tau", "GenPart", ["pt", "eta", "phi", "pdgId", "statusFlags"], None, V2)
    flags = gp["statusFlags"]
    sel = (
        ((abs(gp["pdgId"]) == 11) | (abs(gp["pdgId"]) == 13))
        & (gp["pt"] > 10)
        & (((flags & 1) > 0) | ((flags & (1 << 5)) > 0))
        & ((flags & (1 << 13)) > 0)
    )
    # events without a v2 match point to an appended empty event: no leptons (cannot veto)
    leta, lphi = gp["eta"][sel], gp["phi"][sel]
    empty = ak.unflatten(ak.flatten(leta)[:0], [0])
    rows = np.where(found, row2, len(leta)).astype(np.int64)
    leta = ak.concatenate([leta, empty])[rows]
    lphi = ak.concatenate([lphi, empty])[rows]
    return nearest(jeta, jphi, leta, lphi, 0.3) >= 0, found


def load_sample(sample, n_events=None):
    print(f"loading {sample}")
    j = read(sample, "hltAK4PuppiJetCorr", ["pt", "eta", "phi"], n_events)
    sel = (j["pt"] > 20) & (abs(j["eta"]) < 2.1)
    jpt, jeta, jphi = j["pt"][sel], j["eta"][sel], j["phi"][sel]
    nev = len(jpt)
    df = pd.DataFrame(
        {
            "pt": ak.to_numpy(ak.flatten(jpt)),
            "eta": ak.to_numpy(ak.flatten(jeta)),
        }
    )
    df["aeta"] = df.eta.abs()
    df["event"] = flat_event(jeta)[0]
    # ParT
    p = read(sample, "hltAK4PuppiJetParT", ["eta", "phi", "TauvsAll"], n_events)
    i = nearest(jeta, jphi, p["eta"], p["phi"], 0.05)
    df["part"] = take(ak.to_numpy(ak.flatten(p["TauvsAll"])), i)
    # TEST taggers
    t = read(sample, "hltJetTauTagTEST", ["jetEta", "jetPhi", "ds", "dsA01", "bdt", "scored"], n_events)
    i = nearest(jeta, jphi, t["jetEta"], t["jetPhi"], 0.05)
    scored = take(ak.to_numpy(ak.flatten(t["scored"])).astype(float), i, 0)
    for c in ("ds", "dsA01", "bdt"):
        v = take(ak.to_numpy(ak.flatten(t[c])), i)
        df[c] = np.where((scored > 0) & (v > -0.5), v, MISSING)
    # HPS DeepTau
    h = read(sample, "hltHpsPFTau", ["pt", "eta", "phi", "deepTauVSjet"], n_events)
    i = nearest(jeta, jphi, h["eta"], h["phi"], 0.3)
    hpt = take(ak.to_numpy(ak.flatten(h["pt"])), i)
    hdt = take(ak.to_numpy(ak.flatten(h["deepTauVSjet"])), i)
    df["hps"] = np.where(i >= 0, hdt - hps_medium(hpt), MISSING)
    # L1 nnPuppiTau
    l1 = read(sample, "L1nnPuppiTau", ["eta", "phi", "passLooseNN", "passTightNN"], n_events)
    i = nearest(jeta, jphi, l1["eta"], l1["phi"], 0.3)
    loose = take(ak.to_numpy(ak.flatten(l1["passLooseNN"])).astype(float), i, 0)
    tight = take(ak.to_numpy(ak.flatten(l1["passTightNN"])).astype(float), i, 0)
    df["l1"] = np.where(i >= 0, 1 + (loose > 0) + (tight > 0), MISSING)
    # labels
    df["label"] = "mb" if sample.startswith("MinBias") else "other"
    if not sample.startswith("MinBias"):
        g = read(sample, "GenVisTau", ["pt", "eta", "phi"], n_events)
        i_any = nearest(jeta, jphi, g["eta"], g["phi"], 0.3)
        gpt = take(ak.to_numpy(ak.flatten(g["pt"])), i_any, -1)
        hi = g["pt"] > 20
        i_sig = nearest(jeta, jphi, g["eta"][hi], g["phi"][hi], 0.3)
        df["gen_tau_pt"] = gpt
        lab = np.where(i_sig >= 0, "tau", "other").astype(object)
        if sample == "HHTo2B2Tau":
            gj = read(sample, "GenJet", ["pt", "eta", "phi", "hadronFlavour"], n_events)
            m = gj["pt"] > 15
            ig = nearest(jeta, jphi, gj["eta"][m], gj["phi"][m], 0.4)
            hf = take(ak.to_numpy(ak.flatten(gj["hadronFlavour"][m])).astype(float), ig, -1)
            lep, _ = lepton_veto(jeta, jphi, n_events)
            clean = (i_any < 0) & ~lep
            lab = np.where(clean & (hf == 5), "b", lab)
            lab = np.where(clean & (hf == 4), "c", lab)
            lab = np.where(clean & (hf == 0), "light", lab)
            lab = np.where(~clean & (i_sig < 0), np.where(lep, "lepton", "lowpt_tau"), lab)
            df["lep"] = lep
        df["label"] = lab
    df["sample"] = sample
    print(f"  {nev} events, {len(df)} jets; labels: {df.label.value_counts().to_dict()}")
    for c in TAGGERS:
        print(f"  {c}: scored fraction {np.mean(df[c] > MISSING / 2):.4f}")
    return df


def eff_curve(s, b):
    """signal / background efficiency for every threshold among the observed (non-missing) values."""
    s = np.sort(s)
    b = np.sort(b)
    th = np.unique(np.concatenate([s[s > MISSING / 2], b[b > MISSING / 2]]))
    es = 1 - np.searchsorted(s, th, side="left") / max(len(s), 1)
    eb = 1 - np.searchsorted(b, th, side="left") / max(len(b), 1)
    return th, es, eb


def eff_at(s, b, cut):
    return np.mean(s >= cut), np.mean(b >= cut)


def at_bkg(th, es, eb, target):
    ok = eb <= target
    if not ok.any():
        return np.nan, np.nan, np.nan
    k = np.argmax(np.where(ok, es, -1))
    return es[k], eb[k], th[k]


def select(df, eta, pt):
    lo, hi = ETA[eta]
    plo, phi = PT[pt]
    return df[(df.aeta >= lo) & (df.aeta < hi) & (df.pt >= plo) & (df.pt < phi)]


def draw(ax, sig, bkg, title):
    for c, (lab, col, ls) in TAGGERS.items():
        s, b = sig[c].to_numpy(), bkg[c].to_numpy()
        if len(s) == 0 or len(b) == 0:
            continue
        th, es, eb = eff_curve(s, b)
        if c == "l1":
            ax.plot(es, eb, color=col, ls=ls, marker="D", ms=5, label=lab)
        else:
            if len(th) > 3000:
                k = np.unique(np.linspace(0, len(th) - 1, 3000).astype(int))
                es, eb = es[k], eb[k]
            ax.plot(es, eb, color=col, ls=ls, lw=1.8, label=lab)
        if c == "l1":
            continue
        for (wn, cut), mk in zip(WPS[c], WP_MARKERS):
            e_s, e_b = eff_at(s, b, cut)
            ax.plot(e_s, max(e_b, 1e-5), marker=mk, color=col, ms=8, mfc="white", mew=2, ls="none")
    ax.set_yscale("log")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(1e-4, 1.2)
    ax.grid(True, which="major", alpha=0.3)
    ax.set_xlabel(r"$\tau_h$ efficiency")
    ax.text(0.04, 0.96, title, transform=ax.transAxes, va="top", fontsize=13)


def wp_legend_handles():
    hs = []
    for c, wps in WPS.items():
        if c == "l1":
            continue
        col = TAGGERS[c][1]
        for (wn, _), mk in zip(wps, WP_MARKERS):
            hs.append(
                plt.Line2D([], [], color=col, marker=mk, mfc="white", mew=2, ls="none", ms=8, label=wn)
            )
    return hs


def label(ax):
    hep.cms.label("Phase-2 Simulation Preliminary", data=True, rlabel="PU 200 (14 TeV)", ax=ax, fontsize=16)


def figure_grid(sig, bkg, sig_name, bkg_name, fname, outdir):
    fig, axes = plt.subplots(2, 3, figsize=(24, 15), sharey=True)
    for r, eta in enumerate(("barrel", "endcap")):
        for cc, pt in enumerate(PT):
            s, b = select(sig, eta, pt), select(bkg, eta, pt)
            draw(
                axes[r, cc],
                s,
                b,
                f"{ETA_LABEL[eta]}, {PT_LABEL[pt]}\n{len(s)} $\\tau_h$ ({sig_name}), {len(b)} {bkg_name}",
            )
            if cc == 0:
                axes[r, cc].set_ylabel(f"{bkg_name} efficiency")
            label(axes[r, cc])
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h + wp_legend_handles(), lab + [x.get_label() for x in wp_legend_handles()],
               loc="lower center", ncol=6, fontsize=13, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, outdir, fname)


def figure_summary(sig, bkgs, names, outdir, fname="roc_summary", sig_name="HH"):
    fig, axes = plt.subplots(1, len(bkgs), figsize=(8 * len(bkgs), 9), sharey=True)
    for ax, (bk, bdf) in zip(axes, bkgs.items()):
        s = select(sig, "all", "pt20")
        b = select(bdf, "all", "pt20")
        draw(ax, s, b, f"{ETA_LABEL['all']}, {PT_LABEL['pt20']}\n{len(s)} $\\tau_h$ ({sig_name}), {len(b)} {names[bk]}")
        ax.set_ylabel("background efficiency")
        label(ax)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h + wp_legend_handles(), lab + [x.get_label() for x in wp_legend_handles()],
               loc="lower center", ncol=6, fontsize=13, frameon=False)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    save(fig, outdir, fname)


def save(fig, outdir, fname):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{fname}.{ext}"), dpi=110 if ext == "png" else None)
    plt.close(fig)
    print(f"  wrote {outdir}/{fname}.png/.pdf")


def table(sigs, bkgs):
    rows = []
    for sn, sdf in sigs.items():
        for bk, bdf in bkgs.items():
            for eta in ETA:
                for pt in PT:
                    s, b = select(sdf, eta, pt), select(bdf, eta, pt)
                    for c in TAGGERS:
                        sv, bv = s[c].to_numpy(), b[c].to_numpy()
                        th, es, eb = eff_curve(sv, bv)
                        base = dict(signal=sn, background=bk, eta=eta, pt=pt, tagger=c, n_sig=len(sv), n_bkg=len(bv))
                        ms, mb = np.mean(sv > MISSING / 2), np.mean(bv > MISSING / 2)
                        rows.append(base | dict(point="matched", sig_eff=ms, bkg_eff=mb, cut=np.nan))
                        for t in TARGETS:
                            e_s, e_b, cut = at_bkg(th, es, eb, t)
                            rows.append(base | dict(point=f"bkg{int(round(t * 100))}%", sig_eff=e_s, bkg_eff=e_b, cut=cut))
                        for wn, cut in WPS[c]:
                            e_s, e_b = eff_at(sv, bv, cut)
                            rows.append(base | dict(point=wn, sig_eff=e_s, bkg_eff=e_b, cut=cut))
    out = pd.DataFrame(rows)
    for c in ("sig_eff", "bkg_eff"):
        out[c] = out[c].round(5)
    out["cut"] = out["cut"].round(5)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--outdir", default=WWW)
    p.add_argument("--csv", default=CSV)
    p.add_argument("--max-events", type=int, default=None, help="per sample, for tests")
    p.add_argument("--save-jets", default=None, help="optional parquet with the per-jet table")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    hep.style.use("CMS")
    hh = load_sample("HHTo2B2Tau", args.max_events)
    dy = load_sample("DYLL_M50", args.max_events)
    mb = load_sample("MinBiasF1000", args.max_events)
    if args.save_jets:
        pd.concat([hh, dy, mb]).to_parquet(args.save_jets)
    sigs = {"HH": hh[hh.label == "tau"], "DY": dy[dy.label == "tau"]}
    bkgs = {"MB": mb, "HHb": hh[hh.label == "b"], "HHlight": hh[hh.label == "light"]}
    for sn, sdf in sigs.items():
        for bk, bdf in bkgs.items():
            figure_grid(sdf, bdf, sn, BKG[bk][2], f"roc_tau{sn}_vs_{bk}", args.outdir)
    names = {k: v[2] for k, v in BKG.items()}
    figure_summary(sigs["HH"], bkgs, names, args.outdir)
    # HPS-matched subset (HPS exists only after an L1 tau / e-gamma seed): every tagger on the jets
    # with an HPS tau within 0.3, minimum bias from MinBiasL1Pass (all seeded events of the 2M)
    mbl = load_sample("MinBiasL1Pass", args.max_events)
    hm = lambda d: d[d.hps > MISSING / 2]  # noqa: E731
    sigs_h = {"HH_hpsmatched": hm(sigs["HH"]), "DY_hpsmatched": hm(sigs["DY"])}
    bkgs_h = {"MBL1Pass_hpsmatched": hm(mbl), "HHb_hpsmatched": hm(bkgs["HHb"]),
              "HHlight_hpsmatched": hm(bkgs["HHlight"])}
    names_h = {"MBL1Pass_hpsmatched": "HPS-matched MB (L1Pass) jets",
               "HHb_hpsmatched": "HPS-matched HH b jets", "HHlight_hpsmatched": "HPS-matched HH light/g jets"}
    figure_summary(sigs_h["HH_hpsmatched"], bkgs_h, names_h, args.outdir, "roc_summary_hpsmatched",
                   "HH, HPS-matched")
    figure_grid(sigs_h["HH_hpsmatched"], bkgs_h["MBL1Pass_hpsmatched"], "HH, HPS-matched",
                names_h["MBL1Pass_hpsmatched"], "roc_tauHH_vs_MBL1Pass_hpsmatched", args.outdir)
    t = pd.concat([table(sigs, bkgs), table(sigs_h, bkgs_h)])
    t.to_csv(args.csv, index=False)
    print(f"wrote {args.csv} ({len(t)} rows)")
    key = t[(t.signal == "HH") & (t.pt == "pt20") & (t.point.isin(["matched", "bkg1%", "bkg5%"]))]
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        print(key.pivot_table(index=["background", "eta", "tagger"], columns="point", values="sig_eff"))


if __name__ == "__main__":
    main()
