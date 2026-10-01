"""Read seeds from a dumped CMSSW P2GT menu in the structured form.

Input is a full config dump as produced by `dumpL1TGTMenu.py` in CMSSW (or
exported from ConfDB): conditions reference the object definitions by name
(`object = cms.PSet(refToPSet_ = cms.string("..."))`) and use offline thresholds.
Loading the dump needs `FWCore.ParameterSet.Config`: either run inside `cmsenv`
or pass the path of a CMSSW source checkout (only its pure-python
FWCore/ParameterSet package is used).
"""

import importlib
import os
import sys
import tempfile
from typing import Any, Optional

from menu_tools.menu_translation.object_mapping import CmsswObject, ObjectMapping
from menu_tools.menu_translation.seed_model import (
    CORRELATION_CUTS,
    Leg,
    Seed,
    Threshold,
    TranslationError,
    cmssw_path_to_seed_name,
)

CONDITION_LEGS = {
    "L1GTSingleObjectCond": 1,
    "L1GTDoubleObjectCond": 2,
    "L1GTTripleObjectCond": 3,
    "L1GTQuadObjectCond": 4,
}

# parameters of a condition that are not part of the menu logic
_TECHNICAL = {
    "scales",
    "cosh_eta_lut",
    "cosh_eta_lut2",
    "cos_phi_lut",
    "primVertTag",
    "sanity_checks",
    "inv_mass_checks",
}
_LEG_PARAMETERS = {
    "object",
    "ptScaling",
    "offlineMinPt",
    "offlineMinScalarSumPt",
    "minPt",
    "regionsMinPt",
    "minScalarSumPt",
    "minEta",
    "maxEta",
    "maxPrimVertDz",
    "primVertex",
}


def import_cmssw_config(cmssw_src: Optional[str] = None):
    """Return the FWCore.ParameterSet.Config module."""
    try:
        return importlib.import_module("FWCore.ParameterSet.Config")
    except ImportError:
        if cmssw_src is None:
            raise TranslationError(
                "FWCore.ParameterSet.Config is not importable: run inside cmsenv "
                "or pass the path of a CMSSW source checkout"
            )
    source = os.path.join(
        os.path.abspath(cmssw_src), "FWCore", "ParameterSet", "python"
    )
    if not os.path.isdir(source):
        raise TranslationError(f"{source} does not exist")
    farm = tempfile.mkdtemp(prefix="menu_translation_")
    package = os.path.join(farm, "FWCore", "ParameterSet")
    os.makedirs(package)
    for name in ("FWCore/__init__.py", "FWCore/ParameterSet/__init__.py"):
        open(os.path.join(farm, name), "a").close()
    for entry in os.listdir(source):
        if entry.endswith(".py") and entry != "__init__.py":
            os.symlink(os.path.join(source, entry), os.path.join(package, entry))
    sys.path.insert(0, farm)
    return importlib.import_module("FWCore.ParameterSet.Config")


def load_process(path: str, cmssw_src: Optional[str] = None):
    import_cmssw_config(cmssw_src)
    namespace: dict[str, Any] = {}
    with open(path) as f:
        exec(compile(f.read(), path, "exec"), namespace)
    return namespace["process"]


def _ref_name(pset, what: str) -> str:
    if not pset.isRef_():
        raise TranslationError(f"'{what}' is not a reference to a top-level PSet")
    return pset.refToPSet_.value()


def _value(params: dict, name: str):
    return params[name].value() if name in params else None


def _leg_from_parameters(
    params: dict, process, mapping: ObjectMapping, where: str
) -> Leg:
    if "object" not in params:
        raise TranslationError(
            f"{where}: no 'object' reference (flat menus are not supported, "
            "dump the structured menu)"
        )
    unsupported = set(params) - _LEG_PARAMETERS - _TECHNICAL
    if unsupported:
        raise TranslationError(f"{where}: unsupported parameters {sorted(unsupported)}")

    object_name = _ref_name(params["object"], f"{where}.object")
    scaling = (
        _ref_name(params["ptScaling"], f"{where}.ptScaling")
        if "ptScaling" in params
        else None
    )
    default_scaling = None
    obj_pset = getattr(process, object_name, None)
    if obj_pset is not None and hasattr(obj_pset, "ptScaling"):
        default_scaling = _ref_name(obj_pset.ptScaling, f"{object_name}.ptScaling")
    if scaling == default_scaling:
        scaling = None

    threshold, variable = None, "pt"
    if "offlineMinPt" in params:
        threshold = Threshold(_value(params, "offlineMinPt"), offline=True)
    if "offlineMinScalarSumPt" in params:
        threshold = Threshold(_value(params, "offlineMinScalarSumPt"), offline=True)
        variable = "scalarSumPt"
    online = [n for n in ("minPt", "regionsMinPt", "minScalarSumPt") if n in params]
    if online and threshold is not None or len(online) > 1:
        raise TranslationError(f"{where}: more than one threshold")
    if online:
        values = (
            set(params[online[0]].value())
            if online[0] == "regionsMinPt"
            else {_value(params, online[0])}
        )
        if len(values) != 1:
            raise TranslationError(
                f"{where}: eta-dependent online thresholds {sorted(values)}"
            )
        threshold = Threshold(values.pop(), offline=False, inclusive=False)
        if online[0] == "minScalarSumPt":
            variable = "scalarSumPt"

    key = mapping.from_cmssw(CmsswObject(object_name, scaling, variable))
    leg = Leg(key, threshold)

    min_eta, max_eta = _value(params, "minEta"), _value(params, "maxEta")
    if min_eta is not None or max_eta is not None:
        if min_eta is None or max_eta is None or min_eta != -max_eta:
            raise TranslationError(f"{where}: only symmetric eta ranges are supported")
        leg.abs_eta_max = max_eta
    if "maxPrimVertDz" in params:
        if _value(params, "primVertex") != 0:
            raise TranslationError(f"{where}: only primVertex = 0 is supported")
        leg.pv_dz_max = _value(params, "maxPrimVertDz")
    elif "primVertex" in params:
        raise TranslationError(f"{where}: primVertex without maxPrimVertDz")
    return leg


def _add_correlations(seed: Seed, i: int, j: int, params: dict, where: str):
    unsupported = set(params) - set(CORRELATION_CUTS)
    if unsupported:
        raise TranslationError(
            f"{where}: unsupported correlation cuts {sorted(unsupported)}"
        )
    for cut, param in params.items():
        value = param.value()
        if cut in ("os", "ss") and not value:
            continue
        seed.add_correlation(i, j, cut, value)


def seed_from_condition(name: str, module, process, mapping: ObjectMapping) -> Seed:
    n_legs = CONDITION_LEGS[module.type_()]
    params = module.parameters_()
    seed = Seed(name, [])
    if n_legs == 1:
        leg_params = {k: v for k, v in params.items() if k not in _TECHNICAL}
        seed.legs.append(_leg_from_parameters(leg_params, process, mapping, name))
        return seed

    pairs = [(i, j) for i in range(1, n_legs + 1) for j in range(i + 1, n_legs + 1)]
    known = {f"collection{k}" for k in range(1, n_legs + 1)} | _TECHNICAL
    for k in range(1, n_legs + 1):
        collection = params[f"collection{k}"].parameters_()
        seed.legs.append(
            _leg_from_parameters(collection, process, mapping, f"{name}.collection{k}")
        )
    if n_legs == 2:
        correlations = {k: v for k, v in params.items() if k not in known}
        _add_correlations(seed, 1, 2, correlations, name)
    else:
        for i, j in pairs:
            label = f"correl{i}{j}"
            known.add(label)
            if label in params:
                _add_correlations(
                    seed, i, j, params[label].parameters_(), f"{name}.{label}"
                )
        leftover = set(params) - known
        if leftover:
            raise TranslationError(f"{name}: unsupported parameters {sorted(leftover)}")
    return seed


def read_cmssw_menu(
    path: str, mapping: ObjectMapping, cmssw_src: Optional[str] = None
) -> tuple[list[Seed], dict[str, str]]:
    """(seeds, {skipped algorithm: reason}) of all algorithms of a dumped menu."""
    return read_cmssw_process(load_process(path, cmssw_src), mapping)


def read_cmssw_process(
    process, mapping: ObjectMapping
) -> tuple[list[Seed], dict[str, str]]:
    """(seeds, {skipped algorithm: reason}) of all algorithms of a cms.Process.

    Each algorithm must consist of one path with one P2GT condition; logical
    combinations of paths cannot be expressed in the rate-table menu.
    """
    seeds, skipped = [], {}
    for algo in process.l1tGTAlgoBlockProducer.algorithms:
        expression = algo.expression.value().strip()
        name = (
            algo.name.value()
            if hasattr(algo, "name") and algo.name.value()
            else expression
        )
        try:
            if not hasattr(process, expression) or not hasattr(
                getattr(process, expression), "moduleNames"
            ):
                raise TranslationError(
                    f"logical expression '{expression}' is not supported"
                )
            path_ = getattr(process, expression)
            conditions = [
                m
                for m in path_.moduleNames()
                if getattr(process, m).type_() in CONDITION_LEGS
            ]
            if len(conditions) != 1:
                raise TranslationError(
                    f"path {expression} has {len(conditions)} conditions"
                )
            seeds.append(
                seed_from_condition(
                    cmssw_path_to_seed_name(expression),
                    getattr(process, conditions[0]),
                    process,
                    mapping,
                )
            )
        except TranslationError as e:
            skipped[name] = str(e)
    return seeds, skipped
