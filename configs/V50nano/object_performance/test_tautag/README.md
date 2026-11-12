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

The path-level efficiency (`path_eff_taus.yaml`) needs trigger bits, so there is no tagger
version; the report has the emulated double-τ path. Plot names start with `TEST`.

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
    --version V50nano_170pre6HLT_TESTtautag              # ~1.5 min
rsync -a outputs/V50nano_170pre6HLT_TESTtautag/ \
    /eos/user/a/alobanov/www/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag/
```

Web: https://alobanov.web.cern.ch/L1T/Phase2/menu/Validation/MenuTools_AR26/TEST_tautag/
