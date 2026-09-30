"""Format-neutral description of L1 seeds, shared by all readers and writers."""

from dataclasses import dataclass, field
from typing import Optional, Union

# Correlation cuts between two legs. The names follow the P2GT emulator
# parameters (L1GTCorrelationalCut), all cuts are strict in the emulator.
CORRELATION_CUTS = (
    "minDR",
    "maxDR",
    "minInvMass",
    "maxInvMass",
    "minDEta",
    "maxDEta",
    "minDz",
    "maxDz",
    "os",
    "ss",
)

CorrelationValue = Union[float, bool]


class TranslationError(Exception):
    """A seed cannot be expressed in the target format."""


@dataclass(frozen=True)
class Threshold:
    """pT (or scalar sum pT) threshold of a leg.

    offline: value is an offline threshold, converted with the object scaling;
             otherwise it is applied directly to the online (L1) value.
    inclusive: `>=` (MenuTools convention for offline thresholds) instead of `>`.
               The emulator always applies `>` to the online value.
    """

    value: float
    offline: bool = True
    inclusive: bool = True


@dataclass
class Leg:
    """One object of a seed.

    object: MenuTools object key, e.g. "L1gmtTkMuon:VLoose"
    abs_eta_max: additional |eta| < X requirement
    pv_dz_max: |z0 - z0(primary vertex)| < X requirement
    """

    object: str
    threshold: Optional[Threshold] = None
    abs_eta_max: Optional[float] = None
    pv_dz_max: Optional[float] = None


@dataclass
class Seed:
    """A seed (one condition): its legs and the correlations between them.

    correlations: {(i, j): {cut: value}} with 1-based leg indices i < j
    """

    name: str
    legs: list[Leg]
    correlations: dict[tuple[int, int], dict[str, CorrelationValue]] = field(
        default_factory=dict
    )

    def add_correlation(self, i: int, j: int, cut: str, value: CorrelationValue):
        if cut not in CORRELATION_CUTS:
            raise TranslationError(f"{self.name}: unknown correlation cut '{cut}'")
        if i == j:
            raise TranslationError(f"{self.name}: correlation of leg{i} with itself")
        key = (min(i, j), max(i, j))
        cuts = self.correlations.setdefault(key, {})
        if cut in cuts and cuts[cut] != value:
            raise TranslationError(
                f"{self.name}: conflicting values for {cut} of legs {key}: "
                f"{cuts[cut]} and {value}"
            )
        cuts[cut] = value


def seed_name_to_cmssw(name: str) -> str:
    """MenuTools seed name -> CMSSW condition (module) name: `L1_X` -> `X`"""
    return name[3:] if name.startswith("L1_") else name


def cmssw_path_to_seed_name(path: str) -> str:
    """CMSSW path name -> MenuTools seed name: `pX` -> `L1_X`"""
    return "L1_" + (path[1:] if path.startswith("p") else path)
