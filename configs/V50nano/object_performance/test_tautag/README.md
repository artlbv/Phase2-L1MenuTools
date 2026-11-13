# TEST, not official: tau taggers on HLT jets in the L1 vs HLT tau plots

Experimental tau taggers on HLT AK4 PUPPI jets compared with the L1 NN PUPPI tau, the HLT HPS
taus and the two HLT double-τ paths. They are not HLT objects of any menu, and these plots are not
part of the official configs. Branch `ar26_tautag_test`. It is not part of the upstream PR stack.

**What changed (2026-11-13).**
- The DeepSets is the final model of the hlt-tau-pu report (sections 10.8 and 11), `mcttreg_jv_s0`:
  trained with inclusive tt̄ in the hard-scatter class, with a τ pT regression head. New WPs.
- New score `dsA01` = p_τ / (p_τ + p_MB + 0.1 p_HS), which keeps the b-jet rejection, with its own WPs (`DSA*`).
- New object `hltJetTauTagRegTEST`: the same candidates with pT = the regressed τ pT.
- New sample `MinBiasL1Pass`: the ∧ L1-seed rates with all 2M minimum-bias events (14× the seeded
  statistics of MinBiasF1000).
- WPs at the rate of the new reference path `HLT_DoubleMediumPFPuppiParTTauh30_eta2p1` (`*ParT`).
  The two HLT paths are drawn as reference rates. New path-level comparison of the HLT paths and
  the emulated tagger paths.
- Path level on DY → ττ (no b jets) with paired gained / lost counts and the composition of the
  tagged jets: at the ParT rate the taggers are below ParT on genuine τ (DY 0.074–0.087 vs 0.104);
  their gain over ParT on HH comes mostly from b jets ("Conclusion on τ ID" below).
- Menu-tools options used here: `purity.norm_events`, `purity.seed.only`, `purity.reference_bits`
  (docs/object-performance.md), and in `path_efficiency.py` the bits `HLT_TEST_*` are read from the
  TEST-only record `HLTTEST` (the emulated tagger paths).

## Taggers and objects

- **Taggers:** the hlt-tau-pu study (topo-hackathon-2026, `results/hlt_tau_pu/README.md`). All
  scores are out-of-fold: two-fold cross-fit by file parity.
  - `bdt`: XGBoost BDT on the nano jet variables, PF-candidate and vertex features (`bdt_full`,
    unchanged).
  - `ds`: multi-task DeepSets `mcttreg_jv_s0`. It has 3 classes (τ / minimum bias / hard-scatter
    jet) and runs on the PF candidates and jet variables. It is trained on genuine τ of HH, DY and
    inclusive tt̄, against minimum-bias jets (files 0–999) and HH / DY / tt̄ hard-scatter jets. The
    trigger score is p_τ / (p_τ + p_MB). The class probabilities are `dsPtau`, `dsPmb` and `dsPhs`.
  - `dsA01` = p_τ / (p_τ + p_MB + 0.1 p_HS). In the report it has the same 1 kHz efficiency as `ds`.
    Its HH b-jet fake rate at 80% τ efficiency is 4.3%, against 32% for `ds`.
  - `ptReg` = pT^sig × exp(r̂), from the regression head of the same model. In the report its
    resolution is 10.5% instead of 14.7%, and its response is flat within ±1% for 20–100 GeV.
    `ptRegBDT` is the BDT regressor.
  - `dsOld`: the earlier DeepSets (`mc_jv`), kept for comparison. It is −1 in MinBiasL1Pass.
- **Objects** (made by `menu_tools/caching/make_hlt_jet_tautag_test.py`; definitions in
  `configs/V50nano/objects/test_tautag.yaml`):
  - `hltJetTauTagTEST`: hltAK4PuppiJet candidates with JEC-corrected pT > 20 GeV and |η| < 2.5.
    `pt` is pT^sig (the PF candidates within ΔR < 0.1 of the leading charged candidate);
    `eta` / `phi` are that axis. Unscored jets have `pt` = −1 and scores of −1.
  - `hltJetTauTagRegTEST`: the same jets with `pt` = `ptReg`.
- **Working points** (`scripts/test_tautag/wps.py` → `wps.csv`, `wps_ptreg_scan.csv`). All of them also cut |η| < 2.1.
  - `*1kHz` is the report's double-τ emulation at 1 kHz: `L1_pDoublePuppiTau52_52` ∧ ≥ 2
    candidates with pT^sig > 35 GeV, |η| < 2.1 and ΔR > 0.5. The score cut allows at most 64 of
    the 1,999,040 minimum-bias events.
  - `*ParT` is the same emulation at the rate of `HLT_DoubleMediumPFPuppiParTTauh30_eta2p1`
    (11 events, 0.17 kHz).
  - `*Eff90` / `*Eff80` give 90% / 80% per-candidate efficiency on genuine HH τ_h with
    pT^sig > 20 GeV.
  - "HH eff." is the report's event efficiency at that cut (events with 2 visible τ_h > 25 GeV).
    "MB fake" is the minimum-bias jet fake rate (files 0–999, pT^sig > 20 GeV).

  | score | `1kHz` | `ParT` | `Eff90` | `Eff80` |
  |---|---|---|---|---|
  | `bdt` | > 0.363 (HH eff. 0.597) | > 0.792 (0.513) | > 0.682 (MB fake 5.8%) | > 0.901 (1.6%) |
  | `ds` (new) | > 0.460 (0.614) | > 0.919 (0.515) | > 0.884 (5.1%) | > 0.968 (1.6%) |
  | `dsA01` (`DSA*`) | > 0.398 (0.610) | > 0.857 (0.426) | > 0.795 (6.1%) | > 0.928 (2.2%) |
  | `dsOld` (`DSOld1kHz`) | > 0.412 (0.610) | | | |

  The new DeepSets is +0.017 above the BDT at 1 kHz but equal at the ParT rate. The BDT WPs and the
  old DeepSets' Eff90 / Eff80 (0.866 / 0.970) are reproduced exactly. The previous `DS1kHz` cut
  of 0.424 (61 events) becomes 0.412 at 64 events with the same model.

## Cache and coverage

Everything is in a separate cache version, `V50nano_170pre6HLT_TESTtautag`, under
`$MENU_TOOLS_CACHE`. It holds the tagger objects and symlinks to (or subsets of) the
`V50nano_170pre6HLT` objects the plots need. Running with `--version V50nano_170pre6HLT_TESTtautag`
sends all outputs to `outputs/V50nano_170pre6HLT_TESTtautag/`, apart from the official ones.

- **HHTo2B2Tau** (49,103 events) and **DYLL_M50** (200,000): every event has scores.
- **MinBiasF1000:** the 159,920 events of minimum-bias files 0–999, where every event was scored.
  1 event = 31,038 kHz / 159,920 = 0.194 kHz. The seed fires in 91 events.
- **MinBiasL1Pass** (new): the 31,786 minimum-bias events of the L1 menu OR the three L1 bbττ
  TOPOs at +20 kHz (topo-hackathon-2026 `results/topo_lep/skim/tau_request`).
  - They come from all 1,999,040 events of the MinBias cache. They were scored with the final
    taggers for the TOPO study (`hlt_tau_pu/topo_scores`, `deepsets/topo_minbias.parquet`).
  - **They contain all 1,306 events of the 2M in which `L1_pDoublePuppiTau52_52_final` fired.**
    The cache script checks this against the bit and stops otherwise. So every rate ∧ that seed
    is complete.
  - The rate configs set `norm_events: 1999040`: 1 event = 31,038 kHz / 1,999,040 = 15.5 Hz
    (`<version>_MinBiasL1Pass_norm.json` in the cache). With `seed.only: true` only the ∧ seed
    curves are drawn: the "alone" rates are not defined on this subset.
- **HLTTEST records** (`<version>_<sample>_HLTTEST.parquet`, from
  `scripts/test_tautag/emulate_paths.py`): per-event bits `HLT_TEST_*` of the emulated tagger
  paths, row-aligned with the other records, for HH, MinBiasF1000 and MinBiasL1Pass.
- **TT** is not in the TEST cache. The menu-tools TT sample is tt̄ → 2ℓ2ν; the tagger scored
  inclusive tt̄, which is a different sample.

**Alignment checks.** `make_hlt_jet_tautag_test.py` stops on any failure. For HH, DY,
MinBiasF1000 and MinBiasL1Pass:
- the sorted file lists of the cache and the skim are identical;
- the per-file event counts add up to the 1,999,040 cache rows;
- every scored jet has the same raw pT and η in the cache at `[row][jet_idx]`: 0 mismatches out of
  213,575 / 359,987 / 237,494 / 76,654 jets. The corrected pT agrees to within 6e-8;
- the number of cached jets with corrected pT > 20 GeV and |η| < 2.5 equals the number of scored
  jets.

For MinBiasL1Pass, also:
- its skim (`skim_jet_topo`) holds exactly the 31,786 requested events;
- all 1,306 seed events are among them.

## Plots

| config | content |
|---|---|
| `test_tautag_matching.yaml` | matching efficiency vs gen visible τ pT, barrel / endcap, WPs (`TESTTauMatchingWPs`), α = 0.1 and earlier DeepSets (`TESTTauMatchingWPs_DSA`) |
| `test_tautag_trigger.yaml`, `test_tautag_trigger_wMatching.yaml` | turn-ons at 27 / 35 GeV (pT^sig for the taggers) |
| `test_tautag_eta.yaml` | efficiency vs gen η |
| `test_tautag_l1seed.yaml` | L1 × HLT per leg with `L1_pDoublePuppiTau52_52` (as `hlt_tau_l1seed.yaml`) |
| `test_tautag_resolution.yaml` | response / resolution of pT^sig and of the regressed pT against HPS and L1 |
| `test_tautag_purity.yaml` | HH purity / fake rate; single-object rate vs threshold in MinBiasF1000 |
| `test_tautag_rates.yaml` | single- and double-τ rate vs threshold in MinBiasF1000, each WP alone and ∧ L1 seed |
| `test_tautag_rates_l1pass.yaml` | the ∧ L1 seed rates on MinBiasL1Pass (2M events) |
| `test_tautag_matching_dm.yaml`, `test_tautag_trigger_wMatching_dm.yaml` | split by generator decay mode, one plot per tagger (`_BDT_`, `_DS_`, `_DSA_`) |
| `test_tautag_regpt.yaml` | per-τ turn-on: regressed pT against pT^sig at matched rate |
| `test_tautag_path_eff.yaml` | path level: HLT paths and emulated tagger paths vs the subleading gen τ pT |

Plot names start with `TEST`.

## Rates

`test_tautag_rates.yaml` (MinBiasF1000) and `test_tautag_rates_l1pass.yaml` (MinBiasL1Pass):
`object_performance/purity/TEST{,Di}TauRates{,L1Pass}_{BDT,DS,DSA,DSReg}_event_fraction_vs_threshold_*`.
- The plots show the rate of events with ≥ 1 (`TESTTauRates`) or ≥ 2 (`TESTDiTauRates`) objects
  above the online pT threshold (pT^sig for the taggers). All HLT objects have |η| < 2.1.
- Dashed curves: ∧ `L1_pDoublePuppiTau52_52_final`, the L1 seed of both HLT double-τ paths. Lines:
  the seed alone, 1 kHz and, in the double-τ plots, the two HLT paths.
  - DeepTau `HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1`: 66 events of the 2M (1.02 kHz), and
    3 events in MinBiasF1000.
  - ParT `HLT_DoubleMediumPFPuppiParTTauh30_eta2p1`: 11 events (0.17 kHz), and 0 events in
    MinBiasF1000.
- The double-object rate counts ≥ 2 objects without ΔR between them. The path emulation (below)
  requires ΔR > 0.5, which changes the counts by ≤ 1 event here.

Double-τ rate ∧ L1 seed at 35 GeV, events (kHz, Garwood 68%):

| object | MinBiasF1000 (0.194 kHz / event) | MinBiasL1Pass (15.5 Hz / event) |
|---|---|---|
| L1 seed alone | 91 (17.7) | 1306 (20.3 [19.7, 20.8]) |
| HPS Medium DeepTau (≥ 2 HPS) | 3 (0.58 [0.27, 1.15]) | 65 (1.01 [0.89, 1.15]) |
| `BDT1kHz` | 4 (0.78 [0.41, 1.39]) | 64 (0.99 [0.87, 1.13]) |
| `BDTParT` | 1 (0.19 [0.03, 0.64]) | 11 (0.17 [0.12, 0.24]) |
| `BDTEff90` / `BDTEff80` | 2 / 0 | 23 (0.36 [0.28, 0.45]) / 1 (0.016) |
| `DS1kHz` | 4 (0.78 [0.41, 1.39]) | 65 (1.01 [0.89, 1.15]) |
| `DSParT` | 0 (< 0.36) | 11 (0.17 [0.12, 0.24]) |
| `DSEff90` / `DSEff80` | 1 / 0 | 17 (0.26 [0.20, 0.34]) / 4 (0.062 [0.032, 0.111]) |
| `DSA1kHz` | 4 (0.78 [0.41, 1.39]) | 65 (1.01 [0.89, 1.15]) |
| `DSAParT` | 2 (0.39 [0.14, 0.90]) | 11 (0.17 [0.12, 0.24]) |
| `DSAEff90` / `DSAEff80` | 2 / 0 | 19 (0.30 [0.23, 0.38]) / 2 (0.031 [0.011, 0.072]) |

Single-τ rate at 35 GeV, events (kHz):

| object | MinBiasF1000 alone | MinBiasF1000 ∧ seed | MinBiasL1Pass ∧ seed |
|---|---|---|---|
| HPS Medium DeepTau | 38 (7.38 [6.19, 8.77]) | 34 (6.60 [5.48, 7.93]) | 511 (7.93 [7.59, 8.30]) |
| `BDT1kHz` | 924 (179 [173, 185]) | 30 (5.82 [4.77, 7.09]) | 497 (7.72 [7.37, 8.08]) |
| `DS1kHz` | 1100 (214 [207, 220]) | 27 (5.24 [4.24, 6.45]) | 455 (7.07 [6.74, 7.41]) |
| `DSA1kHz` | 1127 (219 [212, 225]) | 27 (5.24 [4.24, 6.45]) | 457 (7.10 [6.77, 7.44]) |
| `DSEff90` / `DSEff80` | 405 / 156 (78.6 / 30.3) | 14 / 4 (2.72 / 0.78) | 198 / 79 (3.07 / 1.23) |
| `DSParT` | 319 (61.9) | 11 (2.14) | 159 (2.47 [2.27, 2.68]) |

Without the seed, the single-τ rates of the taggers are 30–40× higher than with it.
Lowest double-τ ∧ seed threshold with ≤ 1 kHz on MinBiasL1Pass (`at_rate` in the JSONs):
- HPS Medium 35.5 GeV (64 events), `BDT1kHz` 35.0 GeV (64), `DS1kHz` / `DSA1kHz` 35.5 GeV (63);
- `hltJetTauTagRegTEST` with the `DS1kHz` cut 33.0 GeV (62), with the `DSReg1kHz` cut 26.0 GeV (64);
- the `*Eff90` / `*Eff80` / `*ParT` WPs are below 1 kHz at 10 GeV.

With 160k events, MinBiasF1000 has 3–4 events per WP and cannot rank the WPs. MinBiasL1Pass
reproduces the WP definitions within an event (64 / 65 / 11 events at 35 GeV).

## Path level (HLT paths and emulated tagger paths)

`test_tautag_path_eff.yaml` (stock `path_efficiency` with the `paths:` map, as
`path_eff_taus_compare.yaml`): `object_performance/path_efficiency/TESTDiTauPath_*` (HH) and
`TESTDiTauPathDY_*` (DY).
- Per-event efficiency, no matching, vs the subleading gen visible τ pT, the leading τ pT (events
  with two τ in |η| < 2.1) and the subleading |η| (subleading pT > 60 GeV on HH, > 35 GeV on DY).
  Each plot is also split by the decay-mode pair of the two leading τ (`_DMpair_*`).
- **Emulated paths** `HLT_TEST_Double<tagger>[Reg|Jet]_<WP>_<thr>` = L1 seed ∧ ≥ 2 tagger
  candidates with |η| < 2.1, pT > thr, score > cut, and ΔR > 0.5 between two of them (as the
  report's emulation and the HLT double-τ paths). Threshold on pT^sig; `Reg`: on the regressed pT;
  `Jet`: on the corrected AK4 jet pT (`jetPt`, as the ParT path). Per-event bits in
  `<version>_<sample>_HLTTEST.parquet` (HH, DY, MinBiasF1000, MinBiasL1Pass), written by
  `scripts/test_tautag/emulate_paths.py`; `path_efficiency.py` reads `HLT_TEST_*` bits from that
  record. Cuts, counts and rates: `emulated_paths.csv`.
  - `*_ParT_30`: the 4-decimal cut just above the 12th minimum-bias event, so that 11 of the 2M
    pass (the ParT path's count), at 30 GeV like ParT. `ParTRatePtDef` plots compare the pT
    definitions.
  - `*_1kHz_35`: the 1 kHz WPs at 35 GeV; `DSReg_1kHz_26`: ds > 0.607, regressed pT > 26 GeV.
- **Paired comparison and jet composition**: `scripts/test_tautag/path_diagnosis.py` →
  `path_compare_allacc.csv`, `path_compare_bins.csv`, `path_compare_jets.csv`.
  - Acceptance: ≥ 2 GenVisTau with pT > 25 GeV and |η| < 2.1 (HH 9,595 events, DY 2,764).
  - Gained / lost: events the tagger path takes and the HLT path does not, and the reverse. p-value:
    exact McNemar (binomial on gained vs gained + lost).
  - Tagged pair: of the passing candidates, the pair with the highest smaller score (the one that
    fires the bit). Each jet is classified in order: τ (ΔR < 0.3 to a GenVisTau); lepton (ΔR < 0.3
    to a status-1 e/μ, prompt or from a prompt τ decay, pT > 5); b / c / light (hadronFlavour of
    the nearest GenJet with ΔR < 0.4); unmatched.
  - "ID-only": the tagger path with a passing pair whose two jets are both τ-matched. The HLT
    paths cannot be classified (no per-jet ParT / DeepTau objects in the nano).
- DY GenVisTau and GenJet were cached for this (`configs/V50nano/cache_objects/
  caching_170pre6_HLT_DYgen.yaml`, new files in V50nano_170pre6HLT, same sorted file list and row
  order; 200,000 rows as the other DY records; 99.6% of the GenVisTau mothers in the cached GenPart
  are τ).

Rates: minimum-bias events of the 2M (MinBiasL1Pass), 1 event = 15.5 Hz, 68% Garwood intervals.
The `_ParT_30` paths: 11 events, 0.17 [0.12, 0.24] kHz each; the `_1kHz` paths: 64–65 events,
0.99–1.01 [0.87, 1.15] kHz. On MinBiasF1000 (0.19 kHz per event): DeepTau 3, ParT 0, `_1kHz`
4 each, `_ParT_30` 0 (DSA 1).

### HH → bbττ

Plateau: subleading gen τ pT ≥ 80 GeV (938 events). All (|η|): two gen τ in |η| < 2.1, any pT
(17,690 events, from the leading-pT plot). Acc.: the acceptance above (pT > 25 GeV, 9,595 events),
with gained / lost against the HLT path of the same rate.

| path | cut | min. bias ev. | plateau ≥ 80 | 40–45 GeV | all (\|η\|) | acc. | gained / lost | ID-only acc. |
|---|---|---|---|---|---|---|---|---|
| L1 seed `pDoublePuppiTau52_52` | | 1306 (20.3 kHz) | 0.855 | 0.759 | 0.484 | | | |
| HLT DeepTau HPS35 | | 66 (1.02 [0.90, 1.17] kHz) | 0.710 | 0.517 | 0.272 | 0.429 | | |
| HLT ParT Tauh30 | | 11 (0.17 [0.12, 0.24] kHz) | 0.771 | 0.611 | 0.304 | 0.515 | | |
| `DoubleBDT_ParT_30` | 0.808 | 11 | 0.799 | 0.608 | 0.353 | 0.533 | 724 / 552 | 0.436 |
| `DoubleBDTReg_ParT_30` | 0.8152 | 11 | 0.797 | 0.597 | 0.339 | 0.521 | 683 / 623 | 0.441 |
| `DoubleBDTJet_ParT_30` | 0.8367 | 11 | 0.809 | 0.603 | 0.362 | 0.544 | 816 / 537 | 0.462 |
| `DoubleDS_ParT_30` | 0.919 | 11 | 0.796 | 0.606 | 0.357 | 0.542 | 769 / 513 | 0.443 |
| `DoubleDSReg_ParT_30` | 0.9341 | 11 | 0.788 | 0.579 | 0.324 | 0.511 | 655 / 694 | 0.439 |
| `DoubleDSJet_ParT_30` | 0.9464 | 11 | 0.794 | 0.580 | 0.351 | 0.533 | 799 / 628 | 0.453 |
| `DoubleDSA_ParT_30` | 0.8598 | 11 | 0.713 | 0.537 | 0.272 | 0.464 | 348 / 832 | 0.440 |
| `DoubleBDT_1kHz_35` | 0.363 | 64 | 0.838 | 0.687 | 0.418 | 0.597 | 1741 / 130 (vs DeepTau) | 0.454 |
| `DoubleDS_1kHz_35` | 0.460 | 65 | 0.839 | 0.703 | 0.434 | 0.614 | 1866 / 92 | 0.463 |
| `DoubleDSA_1kHz_35` | 0.398 | 65 | 0.836 | 0.700 | 0.430 | 0.610 | 1841 / 109 | 0.463 |
| `DoubleDSReg_1kHz_26` | 0.607 | 64 | 0.837 | 0.690 | 0.422 | 0.616 | 1923 / 130 | 0.544 |

Per bin of the subleading τ pT, `DoubleDS_ParT_30` against ParT (gained / lost, McNemar p):

| subl. pT (GeV) | 25–30 | 30–35 | 35–40 | 40–45 | 45–50 | 50–60 | 60–80 | ≥ 80 |
|---|---|---|---|---|---|---|---|---|
| events | 1678 | 1411 | 1257 | 1054 | 800 | 1248 | 1209 | 938 |
| ParT | 0.107 | 0.256 | 0.484 | 0.611 | 0.682 | 0.734 | 0.796 | 0.771 |
| DS | 0.183 | 0.310 | 0.498 | 0.606 | 0.680 | 0.736 | 0.807 | 0.796 |
| gained / lost | 180 / 53 | 137 / 60 | 106 / 88 | 76 / 81 | 57 / 59 | 81 / 78 | 75 / 61 | 57 / 33 (p = 0.015) |
| DS ID-only | 0.062 | 0.218 | 0.385 | 0.519 | 0.594 | 0.635 | 0.728 | 0.700 |

**Jet composition** of the tagged pair, HH in acceptance (fractions of jets; "pair with b": events
with at least one b-matched leg). Leptons, c and unmatched are ≤ 1% everywhere.

| path | events | τ | b | light | ττ pair | pair with b |
|---|---|---|---|---|---|---|
| `DoubleDS_ParT_30`, all fired | 5196 | 0.85 | 0.13 | 0.01 | 0.71 | 0.26 |
| `DoubleDS_ParT_30`, gained vs ParT | 769 | 0.67 | 0.31 | 0.02 | 0.36 | 0.60 |
| `DoubleBDT_ParT_30`, gained | 724 | 0.67 | 0.29 | 0.02 | 0.37 | 0.57 |
| `DoubleDSJet_ParT_30`, gained | 799 | 0.74 | 0.23 | 0.02 | 0.50 | 0.46 |
| `DoubleDSA_ParT_30`, all fired | 4456 | 0.96 | 0.03 | 0.01 | 0.93 | 0.05 |
| `DoubleDSA_ParT_30`, gained | 348 | 0.88 | 0.08 | 0.02 | 0.77 | 0.17 |
| `DoubleDS_1kHz_35`, gained vs DeepTau | 1866 | 0.68 | 0.29 | 0.02 | 0.39 | 0.57 |

### DY (DYLL_M50, no b jets)

Acceptance as above (2,764 events; few have a subleading τ above 60 GeV). Plateau: subleading
τ pT ≥ 45 GeV (80 events).

| path | acc. | gained / lost vs ParT (DeepTau) | McNemar p | ID-only | 35–40 GeV | 40–45 GeV | ≥ 45 GeV |
|---|---|---|---|---|---|---|---|
| L1 seed | 0.158 | | | | | 0.581 | 0.800 |
| HLT ParT Tauh30 | **0.104** (287) | | | | 0.228 | 0.476 | 0.662 |
| HLT DeepTau HPS35 | 0.059 (163) | | | | | | 0.562 |
| `DoubleBDT_ParT_30` | 0.074 | 18 / 100 | < 1e-4 | 0.073 | 0.185 | 0.339 | 0.512 |
| `DoubleBDTReg_ParT_30` | 0.082 | 29 / 89 | < 1e-4 | 0.081 | | | |
| `DoubleBDTJet_ParT_30` | 0.087 | 39 / 85 | < 1e-4 | 0.086 | 0.204 | 0.339 | 0.512 |
| `DoubleDS_ParT_30` | 0.076 | 20 / 97 | < 1e-4 | 0.074 | 0.185 | 0.387 | 0.512 |
| `DoubleDSReg_ParT_30` | 0.081 | 34 / 96 | < 1e-4 | 0.080 | 0.201 | 0.347 | 0.512 |
| `DoubleDSJet_ParT_30` | 0.082 | 36 / 97 | < 1e-4 | 0.080 | 0.188 | 0.323 | 0.487 |
| `DoubleDSA_ParT_30` | 0.080 | 17 / 83 | < 1e-4 | 0.079 | 0.198 | 0.395 | 0.525 |
| `DoubleBDT_1kHz_35` | 0.079 | 79 / 25 (DeepTau) | < 1e-4 | 0.073 | | | 0.650 |
| `DoubleDS_1kHz_35` | 0.080 | 79 / 21 (DeepTau) | < 1e-4 | 0.076 | | | 0.675 |
| `DoubleDSReg_1kHz_26` | 0.116 | 183 / 24 (DeepTau) | < 1e-4 | 0.114 | | | 0.700 |

- In DY the tagged pairs in acceptance are τ-matched (98% ττ pairs), so ID-only ≈ fired.
- Over all DY events the tagger paths fire in 4.6–6% of the events (α = 0.1: 1.0–1.3%), 96% of
  the tagged jets being prompt e / μ (DY → ee, μμ); the HLT ParT and DeepTau paths fire in 5.4% and 5.1% of all DY events
  but only 287 and 163 times in the ττ acceptance, so they too fire mostly on leptons.

### Conclusion on τ ID

- **The taggers do not beat ParT on genuine τ ID at the ParT rate.** On DY (no b jets) every
  tagger path is below ParT: 0.074–0.087 against 0.104, losing 3–5 times more events than it gains
  (p < 1e-4), and the deficit is largest at 40–60 GeV (0.32–0.39 vs 0.48–0.68). The cut on the jet
  pT (as ParT) or on the regressed pT recovers only 0.005–0.013 of it.
- **On HH the advantage over ParT is mostly b jets.** In the acceptance the DS path gains 769 and
  loses 513 events, but 60% of the gained events have a b-matched leg. Below 35 GeV the gain is
  large (0.18 vs 0.11 at 25–30 GeV) while the ID-only efficiency is below ParT (0.06). At
  35–80 GeV the paths are equal (p > 0.2). Only the HH plateau ≥ 80 GeV keeps a +0.025 gain
  (57 / 33, p = 0.015) whose composition is not separated here.
- The α = 0.1 score removes the b jets (5% of the pairs with a b leg) and then sits below ParT on
  HH (0.464 vs 0.515, 348 / 832).
- **Against DeepTau at 1 kHz the taggers win on genuine τ**: on DY 0.079–0.116 vs 0.059
  (79–183 gained / 21–25 lost), and on HH ID-only 0.45–0.54 vs DeepTau's 0.43 (all its
  events, unclassified).
- Caveats: DY ττ statistics are small above 50 GeV (80 events above 45 GeV); ParT's own pairs
  cannot be classified, so ParT may also fire on b jets on HH; the ParT rate rests on 11
  minimum-bias events.

## Regressed pT

`test_tautag_regpt.yaml` (per-τ matching + threshold efficiency, the pT cut in the id):
`TESTTauRegPt_{Barrel,Endcap,Inclusive}_*`. Plus `hltJetTauTagRegTEST` in the resolution and
`TESTDiTauRates{,L1Pass}_DSReg` plots.

| emulated double-τ path | MB events (2M) | HH event eff. | per-τ eff., gen 33–36 GeV | 36–39 | 70–150 |
|---|---|---|---|---|---|
| pT^sig > 35, ds > 0.460 (`DS1kHz`) | 64 | 0.614 | 0.39 | 0.62 | 0.942 |
| ptReg > 35, ds > 0.460 | 56 | 0.555 | 0.35 | 0.64 | 0.949 |
| ptReg > 33, ds > 0.460 (matched rate) | 62 | 0.575 | 0.57 | 0.72 | 0.950 |
| ptReg > 35, ds > 0.407 (1 kHz) | 64 | 0.557 | 0.36 | 0.64 | 0.950 |
| **ptReg > 26, ds > 0.607 (1 kHz)** | 64 | **0.616** | 0.79 | 0.81 | 0.951 |

Per-τ efficiencies are inclusive in |η| < 2.1.

- **To match pT^sig > 35 GeV at 1 kHz, the regressed pT needs about 26 GeV.**
  - The cut is re-chosen per threshold (`wps_ptreg_scan.csv`). The HH efficiency at 1 kHz is 0.616
    at 26 GeV, 0.614 at 26.5 GeV and 0.609 at 28 GeV, against 0.614 for pT^sig > 35. This
    reproduces the report's "≈ 25–28 GeV".
  - With the score cut fixed, the matched-rate threshold is 33 GeV, with a lower efficiency (0.575).
- The regressed pT lowers the minimum-bias candidates more than the τ_h: 56 instead of 64 events
  at 35 GeV with the same score cut.
- In the endcap the regression raises the plateau (70–150 GeV, `DS1kHz` cut) from 0.868 to 0.902.
  pT^sig is 12% low at 1.5 < |η| < 2.1.

## Decay-mode split

`test_tautag_{matching,trigger_wMatching}_dm.yaml`: one plot per tagger
(`TESTTau..._{BDT,DS,DSA}[_wMatching]_DM`) with its three WPs, L1 nnPuppiTau and HPS.
`split_by` is on `GenVisTau_status`. The decay-mode shares of the gen τ (|η| < 2.1, barrel) are:
1-prong 16%, 1-prong + π0s 54%, 3-prong 16%, 3-prong + π0 8%, other 6%.

Matching + trigger efficiency at 35 GeV (pT^sig for the taggers), gen pT 70–150 GeV, in %
(inclusive / 1-prong / 1-prong + π0s / 3-prong / 3-prong + π0):

| object | barrel | endcap |
|---|---|---|
| L1 nnPuppiTau | 95.0 / 92.8 / 96.3 / 94.8 / 93.9 | 83.1 / 71.3 / 85.5 / 86.2 / 86.7 |
| HPS Medium (L1 τ events) | 49.0 / 47.4 / 50.2 / 48.3 / 50.0 | 24.5 / 29.0 / 23.9 / 26.3 / 23.6 |
| `BDT1kHz` | 95.0 / 92.0 / 94.8 / 97.4 / 96.6 | 85.5 / 84.6 / 83.7 / 92.6 / 88.7 |
| `BDTEff90` | 93.3 / 91.7 / 93.7 / 94.9 / 92.8 | 83.5 / 83.5 / 82.6 / 90.2 / 80.8 |
| `BDTEff80` | 89.5 / 90.1 / 90.5 / 90.2 / 86.1 | 78.3 / 82.4 / 78.4 / 80.9 / 70.9 |
| `DS1kHz` | 95.2 / 92.3 / 95.0 / 97.2 / 97.0 | 85.2 / 84.2 / 83.4 / 91.8 / 87.7 |
| `DSEff90` | 92.7 / 91.2 / 93.5 / 93.2 / 91.3 | 83.3 / 84.2 / 82.5 / 87.5 / 81.8 |
| `DSEff80` | 88.1 / 87.7 / 90.2 / 86.5 / 83.7 | 80.2 / 83.5 / 80.4 / 83.3 / 72.4 |
| `DSA1kHz` | 95.2 / 92.3 / 95.0 / 97.3 / 97.0 | 85.2 / 84.2 / 83.4 / 91.8 / 87.7 |
| `DSAEff90` | 92.7 / 91.1 / 93.4 / 93.4 / 91.9 | 83.1 / 83.8 / 82.1 / 89.4 / 80.3 |
| `DSAEff80` | 87.4 / 88.6 / 87.9 / 88.2 / 84.9 | 78.7 / 80.9 / 79.2 / 82.0 / 71.9 |

Matching efficiency, gen pT 40–150 GeV, barrel, in % (inclusive / 1-prong / + π0s / 3-prong / + π0):
- `DS1kHz`: 94.7 / 94.4 / 94.4 / 95.9 / 95.4
- `DSEff80`: 83.8 / 85.9 / 86.0 / 80.8 / 76.7
- `BDTEff80`: 84.5 / 89.5 / 85.1 / 84.4 / 77.9

**Plateaus.**
- The loose (`*1kHz`) WPs are flat in decay mode and ≥ 95% (barrel) / 85% (endcap) for every
  tagger.
- The new DeepSets has the same plateau as the earlier one (95.3 → 95.2% barrel at `DS1kHz`).
  Its Eff90 / Eff80 cuts are slightly tighter than before, by 0.3 / 0.4 points in the barrel.
- The tight WPs lose most on 3-prong + π0 (84–86% barrel, 71–72% endcap).

## Caveats

- These are not HLT objects and are not in the HLT menu. The scores are offline, out-of-fold
  predictions from the nano.
- hltHpsPFTau is stored only in events where an L1 τ seed fired, so its efficiencies are
  L1 × HLT-like and its single-object rate includes the L1 seed. The HLT-jet candidates exist in
  every event. Compare the taggers with HPS ∧ seed or at the path level.
- The ParT path can only be compared at the path level: the nano has no ParT scores per jet.
- The minimum bias of files 0–999 was the minimum-bias training sample (each event scored by the
  model of the other fold). In the report, the held-out files agree within 0.001.
- The WPs and the DeepSets variant were chosen on the same HH / minimum-bias events
  (≈ 30 variants in the report).
- MinBiasL1Pass: the DeepSets scores agree with the main skim for the jets in both. The BDT
  differs for 1.4% of the jets (by up to 0.11), from float rounding of the jet-mass feature
  (topo_scores README).

## Reproduce

From the worktree, after `source ../../env.sh` (1 core; the caches take about 1 min, the plots
about 5 min):

```bash
cache_objects configs/V50nano/cache_objects/caching_170pre6_HLT_DYgen.yaml   # DY GenVisTau, GenJet (3 min)
S=menu_tools/caching/make_hlt_jet_tautag_test.py
nice python3 $S --sample HHTo2B2Tau --skim-sample hh
nice python3 $S --sample DYLL_M50 --skim-sample dy
nice python3 $S --sample MinBias --skim-sample minbias --max-file 1000 --out-sample MinBiasF1000
nice python3 $S --sample MinBias --skim-sample minbias --events-l1pass --out-sample MinBiasL1Pass
nice python3 scripts/test_tautag/wps.py                 # WPs -> wps.csv, wps_ptreg_scan.csv (20 s)
nice python3 scripts/test_tautag/emulate_paths.py       # HLTTEST records, path rates (1 min)
nice python3 scripts/test_tautag/path_diagnosis.py      # paired comparison, jet composition (1 min)
nice object_performance configs/V50nano/object_performance/test_tautag/*.yaml \
    --version V50nano_170pre6HLT_TESTtautag
rsync -a outputs/V50nano_170pre6HLT_TESTtautag/ \
    /eos/user/a/alobanov/www/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag/
```

Web: https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag/
