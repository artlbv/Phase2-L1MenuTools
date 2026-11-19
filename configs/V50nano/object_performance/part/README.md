# HLT ParT (AK4 ParticleTransformer) for τ and b tagging, p2nano v2

This folder holds the HLT ParT tagger on HLT AK4 PUPPI jets, added to the menu-tools τ and b-tag plots.
It uses only the p2nano **v2** L1+HLT production (CMSSW_17_0_0_pre6, `P2Nano_L1HLT_noSel_v2`).
v1 stored ParT only in events where the ParT path ran, so v1 ParT is biased.
v2 stores ParT and DeepFlavour for every HLT jet with pT > 5 and |η| < 2.1 (`ParT_inputFilled`). Other jets get −10.

Branch: `ar26_part`, made from `ar26_integration` with `ar26_tautag_test` and `ar26_s12_sample_label` merged in.
Outputs:
- default plots: https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/V50nano_170pre6HLTv2/
- TEST plots: https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag_v2/

## Samples and cache

The cache version is `V50nano_170pre6HLTv2` (`configs/V50nano/cache_objects/caching_170pre6_HLTv2.yaml`). It has the v1 objects plus all ParT and DeepFlavour branches, hadron/parton flavour, GenVisTau, GenJet and GenPart.

| Sample | Events | Note |
|---|---|---|
| MinBias | 1,998,400 (2,094 files) | rate normalisation, 15.53 Hz/event |
| DYLL_M50 | 985,000 | 5× v1 (v1 had 200k) |
| TT | 300,000 | **inclusive tt̄**; v1 was tt̄ → 2ℓ2ν, so these are different events |
| HHTo2B2Tau | 49,103 | the same events as v1 |
| HHTo4B | 50,000 (50 files) | the same 50,000 events as v1 (task `261116_104350`) |

Post-processing steps:
- `make_hlt_jets_corr`, now streamed by row group, writes `hltAK4PuppiJetCorr`. It checks the JEC factor against the nano `jecFactor`: they agree to < 1.4e-4, except 74 of 227k jets at |η| bin edges.
- `make_hlt_jets_part` (new) writes `hltAK4PuppiJetParT`. It keeps jets with `ParT_inputFilled` and raw pT > 10. Its pT is the **raw** pT, as in the HLT path; `ptCorr` holds the corrected pT.
- `scripts/part/event_ids.py` writes `<version>_<sample>_EventID.parquet` with (run, luminosityBlock, event).
- `scripts/part/part_wps.py` sets the ParT working points (WPs) and writes the `HLTEMU` record of emulated ParT paths (bits `HLT_EMU_*`). `path_efficiency.py` reads `HLT_EMU_*` bits from this record.
- HHTo4B was added later (2026-11-19) with the same objects and steps. Its `HLTEMU` record (`scripts/part/hltemu_sample.py`) uses the WPs of `part_wps.csv` unchanged (the emulated Medium path reproduces the HLT ParT bit in 172 / 172 events). All 52 records have 50,000 rows. The JEC check flags 202 of 540k jets (HHTo2B2Tau: 209 of 520k), the same |η|-edge jets as in the other samples.
- Independent check against the HH→4b skim of `results/hlt_hh/skim_v2/` (topo-hackathon repo, same sorted file list): the 31 HLT bits agree in all events, the jet count (raw pT > 10, |η| < 2.5) agrees in all 50,000 events, and pT, η, DeepFlavour, ParT and hadronFlavour agree exactly in the 300 events compared.

## ParT τ ID

The reference path is `HLT_DoubleMediumPFPuppiParTTauh30_eta2p1`. It requires ≥ 2 **uncorrected** `hltAK4PFPuppiJets` with pT > 30 and |η| < 2.1 and `TauvsAll` > 0.11645. Its L1 seed is `L1_pDoublePuppiTau52_52`. It has no ΔR matching to the seed.
The emulation from the stored jets reproduces the HLT bit:
- HH: 8433/8433 events;
- minimum bias: 11/11;
- TT: 7852/7852;
- DY: 52,495 events in both; 5 only in HLT and 2 only in the emulation, from the 10-bit stored precision.

WPs (`objects/hlt.yaml`, object `hltAK4PuppiJetParT`, all with |η| < 2.1). The per-jet efficiency is for HH τh (GenVisTau, ΔR < 0.3); the fake rate is per minimum-bias jet.

| id | cut | definition | HH τ eff. | MB jet fake |
|---|---|---|---|---|
| MediumParT | 0.11645 | the HLT Medium WP | 0.952 | 6.3% |
| ParT1kHz | 0.01834 | ≤ 64 MB events (seed ∧ ≥ 2 jets pT > 30) | 0.976 | 19% |
| ParT1kHz35 | 0.01102 | same at jet pT > 35 | 0.981 | 25% |
| ParTEff90 | 0.4194 | 90% per-jet efficiency | 0.900 | 2.2% |
| ParTEff80 | 0.7651 | 80% per-jet efficiency | 0.800 | 0.75% |

Minimum-bias rates of the full paths (seed ∧ HLT): rate = 31,038 kHz × fraction, with Garwood 68% intervals.

| path | MB events | rate (kHz) |
|---|---|---|
| L1 seed | 1305 | 20.27 [19.71, 20.84] |
| HLT DeepTau (HPS 35) | 66 | 1.03 [0.90, 1.17] |
| HLT ParT Medium (jets 30) | 11 | 0.17 [0.12, 0.24] |
| ParT 1 kHz, jets 30 | 64 | 0.99 [0.87, 1.13] |
| ParT 1 kHz, jets 35 | 64 | 0.99 [0.87, 1.13] |
| ParT Eff90 / Eff80 | 2 / 0 | 0.03 / < 0.03 |

The rates rest on O(10–60) minimum-bias events. The MB/NuGun samples reuse pileup collisions, so the effective count can be lower than the event count. The 64 ParT 1 kHz events have 64 distinct leading tagged jets, so no repeated collision is visible at that level.

Path efficiency in the acceptance of ≥ 2 GenVisTau with pT > 25 and |η| < 2.1, at the same HLT rate of about 1 kHz (68% CL):

| path | HH (9,595 ev.) | DY v2 (13,583 ev.) | DY v1 subset (2,764 ev.) |
|---|---|---|---|
| HLT DeepTau | 0.429 (4120) | 0.057 (778) | 0.059 (163) |
| HLT ParT Medium (0.17 kHz) | 0.515 (4940) | 0.094 (1283) | 0.104 (287) |
| ParT 1 kHz, jets 30 | **0.556** (5335) [0.551, 0.561] | **0.099** (1343) | 0.109 (301) |
| ParT 1 kHz, jets 35 | 0.527 (5061) | 0.077 (1051) | 0.086 (237) |
| best TEST tagger (DSReg 1 kHz, 26) | 0.616 (5913) | — | 0.117 (322) |

Paired against ParT 1 kHz, as gained / lost events:
- HH DeepSets: +742 / −183.
- HH DeepTau: +303 / −1518.
- DY DSReg: +51 / −30.

The TEST taggers gain over ParT on HH mostly from b jets (TEST_tautag README, "Conclusion on τ ID"). On DY only the pT-regressed DeepSets is above ParT.

## ParT b tag (tt̄, v2)

`btag_tt.yaml` now shows two ParT discriminants next to DeepFlavour prob(b) and L1:
- `ParT b/(b+c+uds+g)`;
- `ParT prob(b)`.

New WPs on `hltAK4PuppiJetCorr` at the v2 light mis-id:
- `bTagParTLight10` / `bTagParTLight1` (b/(b+c+uds+g));
- `bTagParTbLight10` / `bTagParTbLight1` (prob(b)).

The block `BTagTT_eta2p1` is the same plot restricted to |η| < 2.1. In |η| < 2.4, ParT leaves jets with |η| 2.1–2.4 untagged.

b-jet efficiency at fixed light mis-id, jet pT > 30:

| region | light mis-id | DeepFlavour prob(b) | ParT b/(b+c+uds+g) | ParT prob(b) | L1 |
|---|---|---|---|---|---|
| \|η\| < 2.4 | 1% | 0.707 | 0.751 | 0.790 | 0.071 |
| \|η\| < 2.4 | 5% | 0.843 | 0.853 | 0.870 | 0.239 |
| \|η\| < 2.4 | 10% | 0.892 | 0.885 | 0.893 | |
| \|η\| < 2.1 | 1% | 0.746 | 0.792 | 0.832 | |
| \|η\| < 2.1 | 5% | 0.871 | 0.900 | 0.920 | |

ParT gains most at high pT: at 1% light mis-id with pT > 100 and |η| < 2.1, ParT reaches 0.864 against 0.722 for DeepFlavour.

## ParT b tag (HH → 4b, v2)

`btag_tt.yaml` has two more blocks on the HHTo4B sample, with the same taggers, WPs and truth definition as `BTagTT`:
- `BTagHH4b` (|η| < 2.4);
- `BTagHH4b_eta2p1` (|η| < 2.1, the ParT acceptance).

`btag_hh4b.yaml` (`BTagHHTo4B`, DeepFlavour and L1 only) is the v1 plot, also run on v2 for the v1/v2 comparison.

b-jet efficiency at fixed light mis-id, jet pT > 30 (`BTagHH4b*_roc_light_pt30toInf`, `target_misid`). HH→4b has 139k b jets but only 43k light jets (|η| < 2.4), so 1% light mis-id rests on about 430 jets.

| region | light mis-id | DeepFlavour prob(b) | ParT b/(b+c+uds+g) | ParT prob(b) | L1 |
|---|---|---|---|---|---|
| \|η\| < 2.4 | 1% | 0.667 | 0.799 | 0.790 | 0.088 (0.3%) |
| \|η\| < 2.4 | 5% | 0.824 | 0.890 | 0.883 | 0.276 (4.6%) |
| \|η\| < 2.4 | 10% | 0.875 | 0.913 | 0.908 | 0.338 (7.5%) |
| \|η\| < 2.1 | 1% | 0.692 | 0.825 | 0.818 | 0.091 (0.3%) |
| \|η\| < 2.1 | 5% | 0.850 | 0.927 | 0.920 | 0.284 (5.0%) |
| \|η\| < 2.1 | 10% | 0.903 | 0.952 | 0.947 | 0.348 (7.9%) |

The L1 score is discrete: its numbers are at the largest mis-id below the target (in brackets).
- The ParT gain over DeepFlavour is larger in HH→4b than in tt̄: at 1% light mis-id and |η| < 2.1, +0.13 (0.825 vs 0.692) against +0.05–0.09 in tt̄.
- In HH→4b, b/(b+c+uds+g) is slightly ahead of prob(b) (0.799 vs 0.790); in tt̄ prob(b) is ahead.
- With pT > 100 and |η| < 2.1, at 1%: ParT 0.839 / 0.831, DeepFlavour 0.636.

Outputs: https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/V50nano_170pre6HLTv2/object_performance/btag/ (`BTagHH4b*`, `BTagHHTo4B*`), and the path table `rate_tables/hlt_paths_V50nano_170pre6HLTv2_HHTo4B.csv`.

## Cross-checks v1 vs v2

- Events are matched by (run, luminosityBlock, event) with `scripts/part/compare_v1_v2.py`:
  - HH has the same 49,103 events.
  - DY v2 is a superset of v1.
  - MB: 1,998,240 events are common; 800 v1-only and 160 v2-only.
- Event by event, the v1 objects are identical except for a few events. The v2 production README says the L1 emulation is not bit-reproducible. Events that differ:

  | sample | HLT bits | L1 seed | HPS τ | jets | L1nnPuppiTau |
  |---|---|---|---|---|---|
  | HH (of 49,103) | 1 | 2 | 20 | 24 | |
  | DY (of 200k) | 2 | 3 | 24 | 141 | |
  | MB | 0 | 1 | 2 | | 776 |

- HH→4b (50,000 events, the same in v1 and v2, matched by event ID):
  - identical in every event: hltAK4PuppiJet (pT, η, φ, DeepFlavour b / c), hltAK4PuppiJetCorr, hltPFPuppiHT, GenJet, GenPart, all 31 HLT bits except `HLT_IsoMu24_FromL1TkMuon` (1 event);
  - L1 differs in a few events: L1puppiJetSC4 200, L1puppiJetSC4NG 99, L1 seeds 1–3 events per bit.
  - `BTagHHTo4B` outputs (152 curves): the HLT curves are identical; L1 curves differ by at most 5e-4 (max |pull| 0.09).
  - `hlt_paths` HHTo4B table: the HLT counts are equal except IsoMu24 (306 vs 307); L1 seed counts differ by at most 3 events.
- Outputs are compared with `scripts/part/compare_outputs.py` → `compare_v1_v2_outputs.csv`, curve by curve against the v1 www jsons:
  - HH, 370 curves: max |pull| 0.28.
  - DY, 221 curves: v2 has 5× the events; max |pull| 3.8 (L1 barrel turn-on), mean χ²/ndf 0.53.
  - The DeepTau and DeepFlavour curves agree within statistics.
- DeepFlavour WPs defined on v1 (tt̄→2ℓ2ν):
  - b efficiency is unchanged on v2: 0.888 / 0.682 at bTagLight10 / bTagLight1.
  - Light mis-id is lower on v2 inclusive tt̄: 0.093 vs 0.100 and 0.0078 vs 0.0100. The sample is different; the b-tag curves are on v2 TT and not comparable 1:1 with v1.

## TEST v2 (`TEST_tautag_v2`)

`scripts/part/make_test_v2.py` builds `V50nano_170pre6HLTv2_TESTtautag`:
- the v1 TEST cache files, as symlinks;
- `hltAK4PuppiJetParT` and `HLTEMU` from v2, put in the TEST row order by event ID.

The 12 MinBiasL1Pass events missing in v2 get no ParT jets. The HLT/L1 bits of the mapped rows are checked against v1 (at most a few events per bit differ). All 13 TEST configs gained ParT MediumParT and ParT1kHz curves, the ParT rates in `rates_l1pass`, and the ParT 1 kHz paths (30 and 35) in `path_eff`.

## τ ROC curves

Per-jet τ ID of all taggers on one common set of jets: `scripts/part/tau_roc.py` → `tau_roc.csv` and the plots in
https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag_v2/object_performance/roc/ (PNG + PDF).
Cache: `V50nano_170pre6HLTv2_TESTtautag` (v2 ParT row-aligned with the TEST taggers; HH 49,103, DY 200,000, MinBiasF1000
159,920 events).
The same plots with only the official objects (HLT ParT, HLT HPS DeepTau, L1 nnPuppiTau; no TEST taggers) are in the
default folder https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/V50nano_170pre6HLTv2/object_performance/roc/
(`--taggers part,hps,l1 --csv none`; same jets and numbers).

**HLT HPS DeepTau is only stored after an L1 seed.** hltHpsPFTau exists in every event where `L1_pDoublePuppiTau52_52`
fired (HH 16,398 / 16,398, DY 17,670 / 17,670, minimum bias 1,306 / 1,306 of the 2M) and in a few events of the
τ-cross seeds, never without one of them. In the unseeded plots the DeepTau curve therefore contains the seed
efficiency (it ends at 0.50 of the HH τh), while ParT and L1 are counted on all jets; the legend says so.

**Seeded version (`--seeded`, `roc_*_seeded`, default v2 folder only).** ParT, DeepTau and the L1 NN score all
∧ `L1_pDoublePuppiTau52_52`, the seed of both HLT double-τ paths (ParT Tauh30 and DeepTau HPS35); all three curves
end at the seed efficiency (0.45–0.47 of the HH τh jets). Minimum bias: the seeded jets are those of MinBiasL1Pass
(all 1,306 seed events of the 2M), divided by all minimum-bias jets of the 2M (MinBiasF1000 jets per η / pT bin
× 1,999,040 / 159,920; check: the HPS-matched MB fraction is 0.00166 in MinBiasF1000 and 0.00174 from MinBiasL1Pass).

HH τh, |η| < 2.1, pT > 20 GeV, all ∧ seed:

| | τh eff. | MB jet eff. | HH b eff. | HH light eff. |
|---|---|---|---|---|
| seed only (ParT-scored jets) | 0.466 | 1.8e-3 | 0.358 | 0.390 |
| ParT MediumParT | 0.449 | 1.9e-4 | 0.031 | 0.019 |
| ParT ParT1kHz | 0.458 | 4.0e-4 | 0.080 | 0.046 |
| HPS Medium DeepTau | 0.403 | 3.9e-4 | 0.065 | 0.036 |
| L1 TightNN | 0.437 | 1.1e-3 | 0.102 | 0.078 |

At 1% HH b (light) efficiency ∧ seed: ParT 0.415 (0.434), DeepTau 0.259 (0.352), L1 0.255 (0.265).
With the seed applied to both, ParT Medium keeps more τh than DeepTau Medium (0.449 vs 0.403) at half the
minimum-bias jet efficiency, and its curve is below DeepTau's everywhere except at the lowest τh efficiencies of the
20–40 GeV endcap bin (few jets).

![HH τh vs minimum bias, b and light jets](https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag_v2/object_performance/roc/roc_summary.png)

*HH τh efficiency against the background efficiency per jet, |η| < 2.1, corrected pT > 20 GeV. Left: MinBiasF1000 jets; middle: HH b jets;
right: HH light/g jets. Open markers: the WPs of the reports applied per jet. ParT (red) is far ahead against b and light jets;
against minimum-bias jets all jet taggers are close, the TEST taggers slightly ahead.*

1. **Denominator.** HLT AK4 PUPPI jets (`hltAK4PuppiJetCorr`) with JEC-corrected pT > 20 GeV and |η| < 2.1. Bins: barrel
   |η| < 1.5, endcap 1.5–2.1; corrected jet pT 20–40 and 40–100 GeV. The pT is the corrected jet pT for every tagger
   (not the raw pT of the ParT path nor the pT^sig of the TEST taggers).
2. **Scores per jet**, each from its own object; a missing score is the lowest value, so a curve ends at its matching
   efficiency (`point = matched` in the CSV):
   - ParT: `hltAK4PuppiJetParT` TauvsAll, nearest jet within ΔR < 0.05. Present for 100% of the jets (0.04% missing in
     MinBiasL1Pass: the 12 events absent in v2).
   - TEST `ds`, `dsA01`, `bdt`: `hltJetTauTagTEST`, matched on `jetEta`/`jetPhi` within ΔR < 0.05. 100% scored.
   - HPS DeepTau: nearest `hltHpsPFTau` within ΔR < 0.3; score = deepTauVSjet − Medium(pT_HPS), the HLT Medium WP
     (0.649 / 0.441 / 0.05 at 35 / 100 / 300 GeV, linear between), so 0 is the Medium WP. HPS-matched fraction:
     HH τh 50%, DY τh 9%, HH b / light jets 44% / 46%, MinBiasF1000 jets 0.17%.
   - L1 nnPuppiTau: nearest within ΔR < 0.3 (any pT); the score is its NN ID score, stored in the nano branch
     `L1nnPuppiTau_chargedIso` (despite the name). Check on all L1 taus of the four samples (508,307): `passLooseNN` ⇔
     score > 0.22 and `passTightNN` ⇔ score > 0.40, with no exception (largest failing / smallest passing values
     0.2197 / 0.2227 and 0.3994 / 0.4033, flat in pT), so the two flags lie exactly on the curve (grey markers).
     L1-matched fraction: HH τh 95%, DY τh 74%, HH b / light 82% / 62%, MinBiasF1000 jets 21%.
3. **Labels.** Signal: jet within ΔR < 0.3 of a GenVisTau with pT > 20 GeV, HH (43,766 jets) and DY (16,474) separately.
   Backgrounds: (a) all MinBiasF1000 jets (147,358; pileup fakes); (b) HH b jets (69,976) and (c) HH light/g jets (47,081):
   nearest GenJet with pT > 15 GeV within ΔR < 0.4 with hadronFlavour 5 / 0, and no GenVisTau (any pT) or prompt e/μ
   (GenPart pT > 10, isPrompt or isDirectPromptTauDecayProduct, last copy) within 0.3. The GenPart comes from the v2
   cache, mapped to the TEST rows by the GenJet content of the event (49,103 / 49,103 mapped). HH c jets (3,488) are not shown.
4. **Signal efficiency at 1% background efficiency** (HH τh, pT > 20 GeV; `tau_roc.csv`, `point = bkg1%`, the highest
   efficiency with background ≤ 1%):

   | tagger | MB all | MB barrel | MB endcap | HH b all | HH b barrel | HH b endcap |
   |---|---|---|---|---|---|---|
   | HLT ParT | 0.837 | 0.840 | 0.790 | **0.750** | **0.778** | **0.575** |
   | TEST DeepSets `ds` | 0.852 | 0.863 | 0.843 | 0.261 | 0.289 | 0.203 |
   | TEST DeepSets `dsA01` | 0.826 | 0.834 | 0.815 | 0.629 | 0.640 | 0.592 |
   | TEST BDT | **0.856** | **0.869** | 0.831 | 0.332 | 0.344 | 0.278 |
   | HPS DeepTau | (0.502) | (0.512) | (0.435) | 0.256 | 0.273 | 0.190 |
   | L1 nnPuppiTau NN score | 0.709 | 0.717 | 0.627 | 0.405 | 0.429 | 0.259 |

   The HPS entry against MB is its matching efficiency: only 0.17% of the minimum-bias jets have an HPS τ, so its
   whole curve is below 1% (item 7). Part of the L1 rejection is likewise the L1 τ reconstruction itself (21% of the
   MB jets have an L1 τ).
   At 5% / 10% background (CSV) the MB numbers converge (ParT 0.944 / 0.965, BDT 0.952 / 0.977) while the b-jet gap
   stays (ParT 0.927 / 0.959, `dsA01` 0.868 / 0.903, `ds` 0.555 / 0.685). HH light/g at 1%: ParT 0.844, `dsA01` 0.766,
   BDT 0.587, L1 0.413, `ds` 0.503, HPS 0.375 (L1 barrel / endcap 0.457 / 0.226).
5. **Against b and light jets ParT is the best tagger by far, against minimum bias the TEST taggers are equal or slightly
   better.** ParT is trained with b / c / light classes; `ds` and the BDT only against minimum bias (+ hard-scatter class
   for `ds`, read as p_τ/(p_τ+p_MB)), so b jets look τ-like to them. `dsA01` recovers most of the b rejection (0.63 at 1%)
   but stays below ParT. This matches the path-level finding that the TEST taggers' HH gain over ParT comes mostly from
   b-jet legs.
6. **In pT bins** (HH, |η| < 2.1, 1% background): against MB ParT / BDT / `ds` reach 0.716 / 0.762 / 0.759 at 20–40 GeV
   and 0.823 / 0.803 / 0.810 at 40–100 GeV; against HH b jets 0.399 / 0.244 / 0.101 (`dsA01` 0.372) at 20–40 and
   0.735 / 0.468 / 0.314 (`dsA01` 0.575) at 40–100 GeV. The TEST advantage against MB is at low pT; ParT is ahead above 40 GeV.
   Endcap is worse for every tagger, most for ParT against b jets (0.575 vs 0.778 barrel).
7. **HPS-matched subset (fair to DeepTau).** HPS taus are stored only in events with an L1 τ or e/γ seed, so on all jets
   the HPS curve is mostly the seed. `roc_summary_hpsmatched.png` and `roc_tauHH_vs_MBL1Pass_hpsmatched.png` use only jets
   with an HPS τ within 0.3, with minimum bias from `MinBiasL1Pass` (31,786 L1-menu/TOPO-pass events of the 2M, containing
   all seeded events; 3,258 HPS-matched jets, so 1% ≈ 33 jets). HH τh at 1% background:
   - MB: ParT 0.796, BDT 0.658, `ds` 0.647, `dsA01` 0.638, DeepTau 0.426;
   - HH b: ParT 0.763, `dsA01` 0.605, BDT 0.343, DeepTau 0.316, `ds` 0.242;
   - On these seeded, harder minimum-bias jets ParT beats the TEST taggers also against minimum bias; DeepTau is the
     weakest against MB in barrel and endcap, and against b jets is comparable to `ds` / the BDT.
8. **L1 nnPuppiTau NN score.** Against MB it is the weakest of the jet-level taggers apart from HPS (0.709 at 1%,
   ParT 0.837, `ds` 0.852), but against HH b jets it is better than `ds` and the BDT (0.405 vs 0.261 / 0.332; ParT
   0.750, `dsA01` 0.629). On the HPS-matched subset against MinBiasL1Pass it drops to 0.185 at 1% (ParT 0.796): on
   seeded, τ-like minimum-bias jets the L1 NN separates little. In pT bins at 1% MB: 0.488 (20–40) / 0.513 (40–100),
   lower than inclusive, where the softer pT spectrum of the MB jets helps.
9. **WPs per jet** (HH τh eff. / MB jet eff. / HH b eff.; all |η| < 2.1, pT > 20): MediumParT 0.951 / 6.0% / 7.9%,
   ParT1kHz 0.976 / 19% / 22%, DS1kHz 0.978 / 9.8% / 90%, DSParT 0.887 / 1.6% / 39%, DSA1kHz 0.979 / 10.5% / 88%,
   BDT1kHz 0.961 / 6.3% / 62%, HPS Medium DeepTau 0.443 / 0.03% / 7.8% (×HPS matching), L1 LooseNN 0.904 / 9.9% / 25%,
   L1 TightNN 0.869 / 5.5% / 17%.
   The TEST WPs were set on the double-τ rate with pT^sig > 35 GeV, so per-jet at pT > 20 they are loose.
10. **DY τh** are softer: at 1% MB background ParT 0.751, `ds` 0.662, BDT 0.653, `dsA01` 0.667; at 1% HH b ParT 0.623,
   `dsA01` 0.390, `ds` / BDT ≈ 0.04 (`roc_tauDY_vs_*.png`).
11. **Caveats.**
    - MinBiasF1000 (files 0–999) is the minimum-bias training sample of the TEST taggers; the scores are out-of-fold
      (two-fold by file parity), but the population is the training one. HH and DY τh were also in their training (out-of-fold).
      ParT was trained on other samples.
    - The TEST WPs were chosen on the same minimum-bias events (rate on the 2M, MinBiasF1000 ⊂ 2M).
    - HPS (and therefore DeepTau) exists only in L1-τ/e-γ-seeded events: on all jets its curve mixes the L1 seed with the
      ID; use the HPS-matched plots for the ID comparison. MinBiasL1Pass jets are selected by the L1 menu.
    - The L1 NN score is read from `chargedIso`; the identity with the pass flags is exact on these samples.
    - Per-jet efficiencies; the trigger-level comparison at equal rate is in the sections above and in the TEST README.

## Reproduce

```bash
export _myschedd_POOL=eossubmit _condor_CONDOR_HOST="tweetybird04.cern.ch, tweetybird03.cern.ch"
condor_submit scripts/condor/caching_170pre6HLTv2_minbias.sub     # MB chunks, 140 jobs
condor_submit scripts/condor/caching_170pre6HLTv2_signal.sub      # DY / TT / HH chunks
condor_submit scripts/condor/caching_170pre6HLTv2_post.sub        # merge, sums, jets Corr, ParT
condor_submit scripts/condor/caching_170pre6HLTv2_signal.sub \
    joblist=scripts/condor/caching_170pre6HLTv2_hh4b_jobs.txt logtag=_hh4b   # HHTo4B only, 2 jobs
condor_submit scripts/condor/caching_170pre6HLTv2_post.sub postlist=scripts/condor/caching_170pre6HLTv2_hh4b_post.txt
for s in HHTo2B2Tau HHTo4B DYLL_M50 MinBias; do                      # EventID records
  python3 scripts/part/event_ids.py configs/V50nano/cache_objects/caching_170pre6_HLTv2.yaml --sample $s; done
python3 scripts/part/part_wps.py                                 # WPs, HLTEMU, part_*.csv
python3 scripts/part/hltemu_sample.py --sample HHTo4B            # HLTEMU of HHTo4B with the same WPs
python3 scripts/part/make_test_v2.py --ids-dir <dir with ids_*_v1.npy>
condor_submit scripts/condor/menu_runs_jobsh.sub joblist=scripts/condor/part_170pre6HLTv2.txt
condor_submit scripts/condor/menu_runs_jobsh.sub joblist=scripts/condor/part_test_170pre6HLTv2.txt
condor_submit scripts/condor/menu_runs_jobsh.sub joblist=scripts/condor/hh4b_170pre6HLTv2.txt   # HH->4b b tag, paths
python3 scripts/part/compare_outputs.py <www>/V50nano_170pre6HLT/object_performance \
    outputs/V50nano_170pre6HLTv2/object_performance --subdirs btag --match BTagHHTo4B
python3 scripts/part/compare_outputs.py <www>/V50nano_170pre6HLT/object_performance \
    outputs/V50nano_170pre6HLTv2/object_performance --subdirs turnons path_efficiency resolution \
    --csv configs/V50nano/object_performance/part/compare_v1_v2_outputs.csv
nice python3 scripts/part/tau_roc.py                             # tau ROCs -> www .../roc, tau_roc.csv (1 core)
nice python3 scripts/part/tau_roc.py --taggers part,hps,l1 --seeded --csv none \
    --outdir <www>/V50nano_170pre6HLTv2/object_performance/roc     # official objects only, + seeded
```
