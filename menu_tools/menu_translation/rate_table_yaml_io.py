"""Read and write MenuTools rate-table menus (configs/<version>/rate_table/*.yml)."""

import re
from typing import Any, Optional

import yaml

from menu_tools.menu_translation.cross_mask_parser import LegCut, parse_cross_mask
from menu_tools.menu_translation.seed_model import (
    Leg,
    Seed,
    Threshold,
    TranslationError,
)

_THRESHOLD_RE = re.compile(
    r"^(?:leg\d+\.)?(offline_pt|pt)\s*(>=|>)\s*([-+]?(?:\d+\.?\d*|\.\d+))$"
)
_LEG_RE = re.compile(r"^leg(\d+)$")

# (quantity, op) -> correlation cut
_PAIR_CUTS = {
    ("dR", ">"): "minDR",
    ("dR", "<"): "maxDR",
    ("mass", ">"): "minInvMass",
    ("mass", "<"): "maxInvMass",
    ("dEta", ">"): "minDEta",
    ("dEta", "<"): "maxDEta",
    ("dz", ">"): "minDz",
    ("dz", "<"): "maxDz",
}


def _fmt(value: float) -> str:
    return str(float(value))


def parse_threshold(cut: Optional[str]) -> Optional[Threshold]:
    if cut is None:
        return None
    match = _THRESHOLD_RE.match(cut.strip())
    if not match:
        raise TranslationError(f"unsupported threshold_cut '{cut}'")
    variable, op, value = match.groups()
    return Threshold(
        float(value), offline=variable == "offline_pt", inclusive=op == ">="
    )


def format_threshold(threshold: Optional[Threshold]) -> Optional[str]:
    if threshold is None:
        return None
    variable = "offline_pt" if threshold.offline else "pt"
    op = ">=" if threshold.inclusive else ">"
    return f"{variable} {op} {_fmt(threshold.value)}"


def _set_once(leg: Leg, attribute: str, value: float, seed: str):
    current = getattr(leg, attribute)
    if current is not None and current != value:
        raise TranslationError(f"{seed}: conflicting {attribute} {current} and {value}")
    setattr(leg, attribute, value)


def seed_from_yaml(
    name: str, config: dict[str, Any], primary_vertex: str
) -> tuple[Seed, list[str]]:
    """Seed from one rate-table menu entry, plus warnings."""
    yaml_legs = sorted(
        (int(m.group(1)), value)
        for key, value in config.items()
        if (m := _LEG_RE.match(key))
    )
    unknown = [k for k in config if not _LEG_RE.match(k) and k != "cross_masks"]
    if unknown:
        raise TranslationError(f"{name}: unsupported keys {unknown}")

    pv_legs = [n for n, leg in yaml_legs if leg["obj"] == primary_vertex]
    if len(pv_legs) > 1:
        raise TranslationError(f"{name}: more than one primary vertex leg")
    pv_leg = pv_legs[0] if pv_legs else None

    seed = Seed(name, [])
    index: dict[int, int] = {}  # yaml leg number -> seed leg index (1-based)
    for n, leg in yaml_legs:
        if n == pv_leg:
            if leg.get("threshold_cut") is not None:
                raise TranslationError(f"{name}: threshold on the primary vertex leg")
            continue
        seed.legs.append(Leg(leg["obj"], parse_threshold(leg.get("threshold_cut"))))
        index[n] = len(seed.legs)

    warnings: list[str] = []
    for mask in config.get("cross_masks") or []:
        cuts, mask_warnings = parse_cross_mask(mask)
        warnings += [f"{name}: {w}" for w in mask_warnings]
        for cut in cuts:
            if isinstance(cut, LegCut):
                if cut.quantity != "absEta" or cut.op != "<" or cut.i == pv_leg:
                    raise TranslationError(f"{name}: unsupported leg cut {cut}")
                _set_once(seed.legs[index[cut.i] - 1], "abs_eta_max", cut.value, name)
                continue
            if (cut.i not in index and cut.i != pv_leg) or (
                cut.j not in index and cut.j != pv_leg
            ):
                raise TranslationError(f"{name}: '{mask}' refers to an undefined leg")
            if pv_leg in (cut.i, cut.j):
                other = cut.j if cut.i == pv_leg else cut.i
                if cut.quantity != "dz" or cut.op != "<":
                    raise TranslationError(
                        f"{name}: only |dz| < X is supported with the primary vertex"
                    )
                _set_once(seed.legs[index[other] - 1], "pv_dz_max", cut.value, name)
            elif cut.quantity == "chargeProduct":
                if cut.value != 0:
                    raise TranslationError(
                        f"{name}: charge product must be compared to 0"
                    )
                cut_name = "os" if cut.op == "<" else "ss"
                seed.add_correlation(index[cut.i], index[cut.j], cut_name, True)
            else:
                seed.add_correlation(
                    index[cut.i],
                    index[cut.j],
                    _PAIR_CUTS[(cut.quantity, cut.op)],
                    cut.value,
                )
    return seed, warnings


def read_rate_table_menu(
    path: str, primary_vertex: str
) -> tuple[list[Seed], dict[str, str], list[str]]:
    """(seeds, {skipped seed: reason}, warnings)"""
    with open(path) as f:
        menu = yaml.safe_load(f)
    seeds, skipped, warnings = [], {}, []
    for name, config in menu.items():
        try:
            seed, seed_warnings = seed_from_yaml(name, config, primary_vertex)
        except TranslationError as e:
            skipped[name] = str(e)
            continue
        seeds.append(seed)
        warnings += seed_warnings
    return seeds, skipped, warnings


def _correlation_masks(i: int, j: int, cut: str, value) -> str:
    li, lj = f"leg{i}", f"leg{j}"
    if cut in ("os", "ss"):
        if not value:
            return ""
        return f"{li}.charge*{lj}.charge {'<' if cut == 'os' else '>'} 0.0"
    op = ">" if cut.startswith("min") else "<"
    quantity = {
        "DR": f"{li}.deltaR({lj})",
        "InvMass": f"({li}+{lj}).mass",
        "DEta": f"abs({li}.eta-{lj}.eta)",
        "Dz": f"abs({li}.z0-{lj}.z0)",
    }[cut[3:]]
    return f"{quantity} {op} {_fmt(value)}"


def seed_to_yaml(seed: Seed, primary_vertex: str) -> dict[str, Any]:
    has_pv = any(leg.pv_dz_max is not None for leg in seed.legs)
    offset = 1 if has_pv else 0
    masks = []
    legs: dict[str, Any] = {}
    if has_pv:
        legs["leg1"] = {"threshold_cut": None, "obj": primary_vertex}
    for k, leg in enumerate(seed.legs, start=1):
        n = k + offset
        legs[f"leg{n}"] = {
            "threshold_cut": format_threshold(leg.threshold),
            "obj": leg.object,
        }
        if leg.pv_dz_max is not None:
            masks.append(f"abs(leg{n}.z0-leg1.z0) < {_fmt(leg.pv_dz_max)}")
        if leg.abs_eta_max is not None:
            masks.append(f"abs(leg{n}.eta) < {_fmt(leg.abs_eta_max)}")
    for (i, j), cuts in sorted(seed.correlations.items()):
        for cut, value in cuts.items():
            mask = _correlation_masks(i + offset, j + offset, cut, value)
            if mask:
                masks.append(mask)
    return {"cross_masks": masks, **legs}


def write_rate_table_menu(seeds: list[Seed], primary_vertex: str) -> str:
    menu = {seed.name: seed_to_yaml(seed, primary_vertex) for seed in seeds}
    return yaml.safe_dump(menu, sort_keys=False, default_flow_style=False)
