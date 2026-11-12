# TEST, not official: tau taggers on HLT jets in the L1 vs HLT tau plots

Experimental tau taggers on HLT AK4 PUPPI jets compared with the L1 NN PUPPI tau and the HLT HPS
taus. They are not HLT objects of any menu, and these plots are not part of the official configs.
Branch `ar26_tautag_test`. It is not part of the upstream PR stack.

- **Taggers:** the hlt-tau-pu study (topo-hackathon-2026, `results/hlt_tau_pu/README.md`).
  - `bdt`: XGBoost BDT on the nano jet variables, PF-candidate and vertex features (`bdt_full`).
  - `ds`: 3-class DeepSets (τ / minimum bias / hard-scatter jet) on the PF candidates and jet
    variables (`mc_jv`). The trigger score is p_τ / (p_τ + p_MB); `dsPhs` is the hard-scatter
    probability.
  - Both are trained on genuine τ from HH → bbττ and DY against minimum-bias jets (files 0–999).
    The scores are out-of-fold: two-fold cross-fit by file parity.
- **Object** `hltJetTauTagTEST`, made by `menu_tools/caching/make_hlt_jet_tautag_test.py`
  (definitions in `configs/V50nano/objects/test_tautag.yaml`).
  - Candidates are hltAK4PuppiJet with JEC-corrected pT > 20 GeV and |η| < 2.5.
  - `pt` is pT^sig, the τ energy estimate (the PF candidates within ΔR < 0.1 of the leading
    charged candidate). Its median response is 1.00. The report's trigger emulation cuts on it.
  - `eta` / `phi` are that axis. Unscored jets have `pt` = −1 and scores of −1.
- **Working points.** All of them also cut |η| < 2.1.

  | id | cut | meaning |
  |---|---|---|
  | `BDT1kHz` / `DS1kHz` | bdt > 0.363 / ds > 0.424 | the report's double-τ emulation at the rate of `HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1` (≈ 1 kHz, 64 / 61 minimum-bias events) |
  | `BDTEff90` / `DSEff90` | bdt > 0.682 / ds > 0.866 | 90% per-candidate efficiency on genuine HH τ_h with pT^sig > 20 GeV |
  | `BDTEff80` / `DSEff80` | bdt > 0.901 / ds > 0.970 | 80% |

## Cache and coverage

Everything is in a separate cache version, `V50nano_170pre6HLT_TESTtautag`, under
`$MENU_TOOLS_CACHE`. It holds the tagger object and symlinks to (or subsets of) the
`V50nano_170pre6HLT` objects the plots need. Running with `--version V50nano_170pre6HLT_TESTtautag`
sends all outputs to `outputs/V50nano_170pre6HLT_TESTtautag/`, apart from the official ones.

- **HHTo2B2Tau** (49,103 events) and **DYLL_M50** (200,000): every event has scores.
- **MinBiasF1000:** the 159,920 events of minimum-bias files 0–999, the only files where every
  event was scored. Files ≥ 1000 were scored only in `L1_pDoublePuppiTau52_52` events, and
  using them would bias the rates.
  - 1 event = 31,039 kHz / 159,920 = 0.194 kHz.
  - The other objects are row-aligned subsets of the full MinBias cache, so the L1 and HPS curves
    use the same events.
- **TT and HHTo4B** are not scored yet. New score files go in through `--scores` (keys
  `sample, file, entry, jet_idx`).

**Alignment checks.** `make_hlt_jet_tautag_test.py` stops on any failure. For HH, DY and minimum
bias:
- the sorted file lists of the cache and the skim are identical;
- the per-file event counts add up to the cache rows;
- every scored jet has the same raw pT and η in the cache at `[row][jet_idx]` (0 mismatches out of
  213,575 / 359,987 / 237,494), and the same corrected pT (to within 6e-8);
- the number of cached jets with corrected pT > 20 GeV and |η| < 2.5 equals the number of scored
  jets.

## Plots

| config | content |
|---|---|
| `test_tautag_matching.yaml` | matching efficiency vs gen visible τ pT, barrel / endcap, and all tagger WPs |
| `test_tautag_trigger.yaml`, `test_tautag_trigger_wMatching.yaml` | turn-ons at 27 / 35 GeV (pT^sig for the taggers) |
| `test_tautag_eta.yaml` | efficiency vs gen η |
| `test_tautag_l1seed.yaml` | L1 × HLT per leg with `L1_pDoublePuppiTau52_52` (as `hlt_tau_l1seed.yaml`) |
| `test_tautag_resolution.yaml` | pT^sig response / resolution against HPS and L1 |
| `test_tautag_purity.yaml` | HH purity / fake rate; single-object rate vs threshold in MinBiasF1000 |
| `test_tautag_rates.yaml` | single- and double-τ rate vs threshold in MinBiasF1000, each WP alone and ∧ L1 seed (below) |
| `test_tautag_matching_dm.yaml`, `test_tautag_trigger_wMatching_dm.yaml` | matching and turn-on efficiencies split by generator decay mode (below) |

The path-level efficiency (`path_eff_taus.yaml`) needs trigger bits, so there is no tagger
version; the report has the emulated double-τ path. Plot names start with `TEST`.

## Rates with and without the L1 seed

`test_tautag_rates.yaml` (`object_performance/purity/TEST{,Di}TauRates_{BDT,DS}_event_fraction_vs_threshold_*`).
Rate of events with ≥ 1 (`TESTTauRates`) or ≥ 2 (`TESTDiTauRates`) objects above the online pT
threshold (pT^sig for the taggers) in MinBiasF1000, all HLT objects |η| < 2.1. Each tagger WP and
HPS Medium is drawn alone (solid) and ANDed with the per-event bit `L1_pDoublePuppiTau52_52_final`
(dashed), the L1 seed of `HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1`. The ∧ seed curves are the
ones to compare with HPS. The plots mark 1 kHz and the seed alone. Options `min_objects`, `seed`
and `rate_lines` of the `purity:` block (docs/object-performance.md).

- The seed fires in 91 events (17.7 kHz, Garwood 68% [15.8, 19.7]): every ∧ seed rate is at
  most that. `HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1` fires in 3 events (0.58 kHz), the same
  3 events as HPS Medium ≥ 2 at 35 GeV.
- The double-object rate is ≥ 2 objects without the ΔR > 0.5 between the legs of the path and of
  the report's emulation. For HPS Medium it makes no difference here (3 = 3 events).
- L1 nnPuppiTau is `default` (|η| < 2.172, the seed acceptance), for reference only.

Rate at 35 GeV, events (kHz, Garwood 68%); 1 event = 0.194 kHz:

| object | ≥ 1 alone | ≥ 1 ∧ L1 seed | ≥ 2 alone | ≥ 2 ∧ L1 seed |
|---|---|---|---|---|
| HPS, no ID | 80 (15.5 [13.8, 17.5]) | | 25 (4.85 [3.89, 6.02]) | |
| HPS Medium DeepTau | 38 (7.38 [6.19, 8.77]) | 34 (6.60 [5.48, 7.93]) | 3 (0.58 [0.27, 1.15]) | 3 (0.58 [0.27, 1.15]) |
| BDT > 0.363 (`BDT1kHz`) | 924 (179 [173, 185]) | 30 (5.82 [4.77, 7.09]) | 19 (3.69 [2.85, 4.74]) | 4 (0.78 [0.41, 1.39]) |
| BDT > 0.682 (`BDTEff90`) | 464 (90.1 [85.9, 94.4]) | 16 (3.11 [2.34, 4.09]) | 7 (1.36 [0.86, 2.09]) | 2 (0.39 [0.14, 0.90]) |
| BDT > 0.901 (`BDTEff80`) | 181 (35.1 [32.5, 37.9]) | 7 (1.36 [0.86, 2.09]) | 0 (< 0.36) | 0 (< 0.36) |
| DeepSets > 0.424 (`DS1kHz`) | 1106 (215 [208, 221]) | 27 (5.24 [4.24, 6.45]) | 32 (6.21 [5.12, 7.51]) | 5 (0.97 [0.55, 1.62]) |
| DeepSets > 0.866 (`DSEff90`) | 425 (82.5 [78.5, 86.7]) | 9 (1.75 [1.18, 2.54]) | 3 (0.58 [0.27, 1.15]) | 2 (0.39 [0.14, 0.90]) |
| DeepSets > 0.970 (`DSEff80`) | 158 (30.7 [28.2, 33.3]) | 2 (0.39 [0.14, 0.90]) | 0 (< 0.36) | 0 (< 0.36) |
| L1 nnPuppiTau (online pT) | 2486 (483 [473, 492]) | | 104 (20.2 [18.2, 22.4]) | |

Lowest threshold (0.5 GeV grid from 10 GeV) with rate ≤ 1 kHz, i.e. ≤ 5 events (`at_rate` in
the JSONs); "≤ 10" means the rate is already below 1 kHz at 10 GeV:

| object | ≥ 1 ∧ L1 seed | ≥ 2 ∧ L1 seed | ≥ 2 alone |
|---|---|---|---|
| HPS Medium DeepTau | 70.5 | 19.5 | 19.5 |
| `BDT1kHz` | > 120 (7 events at 120 GeV) | 13.5 | 52.5 |
| `BDTEff90` | 65.5 | ≤ 10 (3 events) | 42.0 |
| `BDTEff80` | 42.0 | ≤ 10 (0 events) | ≤ 10 (2 events) |
| `DS1kHz` | 80.0 | 27.5 | 71.5 |
| `DSEff90` | 42.0 | ≤ 10 (3 events) | 31.5 |
| `DSEff80` | ≤ 10 (3 events) | ≤ 10 (0 events) | ≤ 10 (1 event) |

With the seed, the double-τ rate of all six WPs is ≤ 1 kHz at 35 GeV. The `*1kHz` WPs give
4–5 events at 35 GeV, HPS Medium 3, consistent with how the WPs were set (the report's
emulation with ΔR > 0.5, at the path's rate). These are counts of a few events: the WPs cannot
be ranked by rate on 160k minimum-bias events.

## Decay-mode split

`test_tautag_{matching,trigger_wMatching}_dm.yaml`: the barrel / endcap matching and
matching + trigger plots (HH → bbττ), one plot per tagger (`TESTTau..._{BDT,DS}[_wMatching]_DM`),
with all three WPs, L1 nnPuppiTau and HPS. `split_by` on `GenVisTau_status`
(docs/object-performance.md): per test object a plot with the inclusive curve and one curve per
generator decay mode (`..._DM_<object>_<thr>_*.png`), and `..._DM_fractions_<thr>_*.json`.
The decay-mode shares of the gen τ (|η| < 2.1, barrel): 1-prong 16%, 1-prong + π0s 54%,
3-prong 16%, 3-prong + π0 8%, other 6%.

Matching + trigger efficiency at 35 GeV (pT^sig for the taggers), gen pT 70–150 GeV, in %
(inclusive / 1-prong / 1-prong + π0s / 3-prong / 3-prong + π0):

| object | barrel | endcap |
|---|---|---|
| L1 nnPuppiTau | 95.0 / 92.8 / 96.3 / 94.8 / 93.9 | 83.1 / 71.3 / 85.5 / 86.2 / 86.7 |
| HPS Medium (L1 τ events) | 49.0 / 47.4 / 50.2 / 48.3 / 50.0 | 24.5 / 29.0 / 23.9 / 26.3 / 23.6 |
| `BDT1kHz` | 95.0 / 92.0 / 94.8 / 97.4 / 96.6 | 85.5 / 84.6 / 83.7 / 92.6 / 88.7 |
| `BDTEff90` | 93.3 / 91.7 / 93.7 / 94.9 / 92.8 | 83.5 / 83.5 / 82.6 / 90.2 / 80.8 |
| `BDTEff80` | 89.5 / 90.1 / 90.5 / 90.2 / 86.1 | 78.3 / 82.4 / 78.4 / 80.9 / 70.9 |
| `DS1kHz` | 95.3 / 92.4 / 95.1 / 97.3 / 96.7 | 85.3 / 84.2 / 83.6 / 91.5 / 87.7 |
| `DSEff90` | 93.0 / 91.5 / 93.8 / 93.4 / 91.8 | 83.2 / 83.5 / 82.4 / 88.6 / 80.8 |
| `DSEff80` | 88.5 / 89.1 / 90.9 / 85.8 / 82.0 | 79.2 / 82.4 / 80.0 / 82.0 / 70.0 |

The loose (`*1kHz`) WPs are flat in decay mode, slightly higher for 3-prong. The tight
(`*Eff80`) WPs lose most on 3-prong + π0 (82–86% barrel, ~70% endcap), DeepSets more than the
BDT on 3-prong. HPS is ~50% (barrel) because hltHpsPFTau is stored only in L1 τ events, flat in
decay mode. Matching efficiencies (gen pT 40–150 GeV) show the same pattern: barrel `BDTEff80`
89.5 / 85.1 / 84.4 / 77.9% and `DSEff80` 88.7 / 87.1 / 80.2 / 75.1% (1-prong / + π0s / 3-prong
/ + π0).

## Caveats

- These are not HLT objects and are not in the HLT menu. The scores are offline, out-of-fold
  predictions from the nano.
- hltHpsPFTau is stored only in events where an L1 τ seed fired, so its efficiencies are
  L1 × HLT-like and its single-object rate includes the L1 seed. The HLT-jet candidates exist in
  every event. Compare the taggers with HPS on the `l1seed` plots, or use the report's
  double-τ trigger emulation.
- MinBiasF1000 was the minimum-bias training sample, with each event scored by the model of the
  other fold. The report's held-out rates (files ≥ 1000) agree, but the rates here have 160k
  events of statistics.
- The WPs and the DeepSets variant were picked on the same HH / minimum-bias events
  (≈ 30 variants in the report).

## Reproduce

From the worktree, after `source ../../env.sh` (1 core, about 1 min in total):

```bash
S=menu_tools/caching/make_hlt_jet_tautag_test.py
nice python3 $S --sample HHTo2B2Tau --skim-sample hh
nice python3 $S --sample DYLL_M50 --skim-sample dy
nice python3 $S --sample MinBias --skim-sample minbias --max-file 1000 --out-sample MinBiasF1000
nice object_performance configs/V50nano/object_performance/test_tautag/*.yaml \
    --version V50nano_170pre6HLT_TESTtautag              # ~4 min (rates 20 s, *_dm 1.5 min)
rsync -a outputs/V50nano_170pre6HLT_TESTtautag/ \
    /eos/user/a/alobanov/www/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag/
```

Web: https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag/
