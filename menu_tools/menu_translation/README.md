# MenuTools <-> CMSSW P2GT menu translation

Translates seeds between the MenuTools rate-table menus
(`configs/<version>/rate_table/*menu*.yml`) and the CMSSW P2GT emulator menus
(`L1Trigger/Configuration/python/Phase2GTMenus/SeedDefinitions/step1_2024`).

```bash
# rate-table yaml -> CMSSW seed cff (no CMSSW needed)
python -m menu_tools.menu_translation yaml-to-cmssw \
    configs/V43nano/rate_table/v38_menu_Step1and2.yml -o l1tGTMenu_generated_cff.py

# CMSSW menu -> rate-table yaml
#   input: structured config dump, `dumpL1TGTMenu.py` in CMSSW or a ConfDB export
python -m menu_tools.menu_translation cmssw-to-yaml l1tGTMenu_cff_fullConfigDump.py \
    -o menu.yml --cmssw-src /path/to/CMSSW/src   # --cmssw-src not needed inside cmsenv
```

Seeds that cannot be expressed in the target format are listed as `SKIPPED` with
the reason and the exit code is 1; all other seeds are written.

## How it works

Both directions go through the format-neutral `Seed` description (`seed_model.py`):
legs (object, threshold, |eta| and primary-vertex dz cuts) and correlations between
legs, named after the emulator cuts (`minDR`, `maxInvMass`, `maxDz`, `os`, ...).

| file | purpose |
|---|---|
| `seed_model.py` | `Seed`, `Leg`, `Threshold` |
| `cross_mask_parser.py` | parses `cross_masks` (restricted grammar, python `ast`) |
| `rate_table_yaml_io.py` | yaml <-> `Seed` |
| `cmssw_menu_reader.py` | dumped CMSSW menu -> `Seed` |
| `cmssw_menu_writer.py` | `Seed` -> CMSSW cff (structured form: `object = gt_ref(...)`, `offlineMinPt`) |
| `object_mapping.py` + `configs/<version>/cmssw_object_map.yaml` | MenuTools object key <-> CMSSW object |

Object cuts are **not** translated (the cut languages differ, e.g. `{hwQual} >= 3`
vs the `qualityFlags` bit mask): `cmssw_object_map.yaml` states which objects
correspond, and both definitions must be kept in sync. Objects without an
entry (e.g. `L1gmtTkMuon:Medium`) make the seed be skipped.

## Supported / not supported

* thresholds: `offline_pt >= X` / `offline_pt > X` (offline, scaled) and `pt > X`
  (online). The emulator applies `>` to the online value, so `>=` vs `>` on
  offline thresholds is not preserved when going through CMSSW.
* `cross_masks`: `legI.deltaR(legJ)`, `(legI+legJ).mass`, `abs(legI.eta-legJ.eta)`,
  `abs(legI.z0-legJ.z0)`, `legI.charge*legJ.charge < 0` (`> 0`), `abs(legI.eta) < X`,
  and `abs(legI.z0-legPV.z0) < X` with a `L1PV:default` leg; combined with `&`.
* not supported: other expressions (e.g. b-tag score sums), more than 4 legs
  (plus the vertex), eta-dependent online thresholds, CMSSW algorithms that
  combine several paths (`pA and pB`), flat (non-structured) CMSSW menus.
* masks like `abs(leg2.z0-leg1.z0) < 1 & (leg2.deltaR(leg3) > 0)` are evaluated by
  python as `abs(...) < (1 & (...))` since `&` binds tighter than `<`; they are
  read as intended with a warning - please add parentheses.

Seed names: `L1_X` in yaml <-> condition `X`, path `pX` in CMSSW.
