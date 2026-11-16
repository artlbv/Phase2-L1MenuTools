"""Purity, fake rate and multiplicity of trigger (L1 or HLT) objects.

Unlike the matching efficiencies (gen object -> trigger object), the loop runs over the
trigger objects: *all* test objects of each event that pass the object ID / eta cuts
(optionally only the leading one), and each is flagged `matched` if any selected gen
object lies within `match_dR` (the test object's `match_dR`, or the `purity.match_dR`
override). The matching is not exclusive: two trigger objects can match the same gen
object.

Modes (`purity.mode`), quantity and binning:
  purity        object-level purity  N(matched) / N(all)                vs pt/eta/pu
  fake_rate     fake fraction        N(unmatched) / N(all) = 1 - purity vs pt/eta/pu
                unmatched objects per event (mean) above threshold
                                                                    vs threshold/pu
                events with >= 1 unmatched object above threshold (fraction,
                in kHz if `rate`)                                   vs threshold/pu
                unmatched objects per event per bin (dN/dx)         vs pt/eta
  multiplicity  as fake_rate, but with all objects and no gen matching (for minimum
                bias, whose gen record is pruned): mean multiplicity, rate of events
                with >= 1 object above threshold (kHz), dN/dx per event.

Thresholds are on the online pT field `pt` of the test object, `pt >= threshold` as in
the rate plots (no offline scalings).

Event-level options: `min_objects` (events with >= N objects, e.g. 2 for a
double-object rate), `seed` (also each listed object ANDed with a per-event trigger
bit, e.g. the L1 seed of an HLT path; `only`: just the ANDed curves), `rate_lines`
(rates in kHz to mark, with the lowest threshold that gives at most that rate), `reference_bits` (per-event trigger
bits, e.g. HLT paths, drawn as reference rates) and
`norm_events` (the sample is a subset of that many events holding every event that can
pass, e.g. the L1-accepted minimum-bias events: fractions and rates are per
`norm_events`).
"""

import os
import re
from typing import Any, Optional, Sequence

import awkward as ak
import numpy as np
import numpy.typing as npt

from menu_tools.object_performance.resolution_collection import make_edges
from menu_tools.rate_table.menu_table import garwood_interval
from menu_tools.object_performance.turnon_collection import (
    ArrayLoader,
    TurnOnCollection,
)
from menu_tools.utils import constants, utils
from menu_tools.utils.objects import (
    CUT_FIELD_RE,
    Object,
    compute_selection_mask_for_object_cuts,
)

MODES = ("purity", "fake_rate", "multiplicity")
# Object-level binning variables (of the trigger object) and event-level ones.
OBJECT_VS = ("pt", "eta", "abseta", "phi")
EVENT_VS = ("pu", "threshold")
VS_VARIABLES = OBJECT_VS + EVENT_VS

DEFAULT_VS = {
    "purity": {"pt": {"min": 0, "max": 100, "step": 5}},
    "fake_rate": {
        "pt": {"min": 0, "max": 100, "step": 5},
        "threshold": {"min": 0, "max": 100, "step": 5},
    },
    "multiplicity": {"threshold": {"min": 0, "max": 100, "step": 5}},
}

# Events are processed in slices of this size for the trigger x gen pairing.
CHUNK = 200_000


def delta_phi(a, b):
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def match_any(trig: ak.Array, gen: ak.Array, dr: float, chunk: int = CHUNK) -> ak.Array:
    """Per trigger object: True if any gen object of the same event is within `dr`.

    trig, gen: jagged records with `eta` and `phi`, one row per event (same length).
    Returns a jagged boolean array with the structure of `trig`.
    """
    out = []
    for start in range(0, len(trig), chunk):
        t = trig[start : start + chunk]
        g = gen[start : start + chunk]
        pairs = ak.cartesian(
            {"t": t[["eta", "phi"]], "g": g[["eta", "phi"]]}, nested=True
        )
        deta = pairs.t.eta - pairs.g.eta
        dphi = delta_phi(pairs.t.phi, pairs.g.phi)
        out.append(ak.any(deta * deta + dphi * dphi < dr * dr, axis=2))
    if not out:
        return ak.Array([])
    return ak.concatenate(out, axis=0) if len(out) > 1 else out[0]


def ratio_with_cp(
    num: npt.ArrayLike, den: npt.ArrayLike
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """num/den with Clopper-Pearson 68% errors (NaN where den == 0)."""
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = np.where(den > 0, num / den, np.nan)
    err = utils.clopper_pearson_err(num, den)
    err = np.where(den > 0, err, np.nan)
    return y, err[0], err[1]


def mean_with_error(values: np.ndarray) -> tuple[float, float]:
    """Mean and its standard error (NaN for an empty sample)."""
    n = len(values)
    if n == 0:
        return np.nan, np.nan
    if n == 1:
        return float(values[0]), np.nan
    return float(np.mean(values)), float(np.std(values, ddof=1) / np.sqrt(n))


def linear_fit(x, y, yerr) -> Optional[dict[str, float]]:
    """Weighted straight-line fit y = a + b x; None with fewer than 2 usable points."""
    x, y, yerr = (np.asarray(v, dtype=float) for v in (x, y, yerr))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    if ok.sum() < 2:
        return None
    w = 1 / yerr[ok] ** 2
    (b, a), cov = np.polyfit(x[ok], y[ok], 1, w=np.sqrt(w), cov="unscaled")
    return {
        "intercept": float(a),
        "slope": float(b),
        "slope_err": float(np.sqrt(cov[0, 0])),
    }


class PurityConfig:
    """The `purity:` block of a plot config.

    purity:
      mode: purity                 # purity | fake_rate | multiplicity
      # gen objects to match to (not for multiplicity); a list is an OR
      reference_gen:
        # parquet name, or an object key with id (e.g. "GenPart:electron")
        object: GenPart
        cuts: ["abs({pdgId}) == 11", "(({statusFlags}>>0)&1) == 1"]
        label: "gen e"
      match_dR: 0.15               # optional, default: the test object's match_dR
      vs:
        pt: {min: 0, max: 100, step: 5}       # trigger-object pT
        eta: {min: -3, max: 3, step: 0.2}     # trigger-object eta (abseta: |eta|, phi)
        pu: {edges: [140, 180, 200, 220, 260]} # pileup of the event (Pileup parquet)
        threshold: {min: 0, max: 100, step: 5} # threshold scan (event-level quantities)
      thresholds: [0, 20]          # pT thresholds for the vs eta/phi/pu plots
                                   # (default [0])
      leading_only: false          # only the leading (highest-pT) test object of
                                   # each event
      rate: false                  # event fractions in kHz (default: true for
                                   # MinBias samples)
      pu_variable: nPU
      pileup_version: null         # cache version holding the Pileup parquet
                                   # (default: the version, then <version>PU)
      min_objects: 1               # event fractions: events with >= N objects above
                                   # threshold (2: double-object rate)
      seed:                        # also each object ANDed with a per-event bit
        bit: L1_pDoublePuppiTau52_52_final  # field of <v>_<sample>_<record>.parquet
        record: L1                 # default: the bit's prefix (L1 / HLT)
        label: "L1 seed"
        objects: [...]             # test object keys (default: all)
        only: false                # true: drop the plain curves of these objects
                                   # (e.g. a sample of seed-passing events only)
      reference_bits:              # rates of per-event bits (e.g. HLT paths) as lines
        - {bit: HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1, label: "DeepTau path"}
                                   # record: default the bit's prefix (L1 / HLT)
      norm_events: null            # event-level fractions / rates per this many events
                                   # instead of the sample's (a subset that contains
                                   # every event that can pass); not with `pu`
      rate_lines: [1000]           # kHz (needs `rate`): lines, and per curve the
                                   # lowest threshold with rate <= line
      rate_line_step: 0.5          # GeV, threshold step of that search
    """

    def __init__(self, cfg: Optional[dict], sample: str) -> None:
        cfg = cfg or {}
        self.mode: str = cfg.get("mode", "purity")
        if self.mode not in MODES:
            raise ValueError(f"purity: unknown mode `{self.mode}`, use one of {MODES}")
        gen = cfg.get("reference_gen")
        if self.mode != "multiplicity":
            if not gen:
                raise KeyError(f"purity: mode `{self.mode}` needs `reference_gen`")
            self.reference_gen: list[dict] = gen if isinstance(gen, list) else [gen]
            for spec in self.reference_gen:
                if "object" not in spec:
                    raise KeyError(
                        "purity: each `reference_gen` entry needs an `object`"
                    )
        else:
            self.reference_gen = []
        self.match_dR: Optional[float] = cfg.get("match_dR")
        vs = cfg.get("vs") or DEFAULT_VS[self.mode]
        for name in vs:
            if name not in VS_VARIABLES:
                raise ValueError(
                    f"purity: unknown `vs` variable `{name}`, use one of {VS_VARIABLES}"
                )
        self.vs: dict[str, Any] = dict(vs)
        self.thresholds: list[float] = [float(t) for t in cfg.get("thresholds", [0])]
        self.leading_only: bool = bool(cfg.get("leading_only", False))
        self.rate: bool = bool(cfg.get("rate", "minbias" in sample.lower()))
        self.pu_variable: str = cfg.get("pu_variable", "nPU").removeprefix("Pileup_")
        self.pileup_version: Optional[str] = cfg.get("pileup_version")
        self.min_objects: int = int(cfg.get("min_objects", 1))
        if self.min_objects < 1:
            raise ValueError("purity: `min_objects` must be >= 1")
        self.seed: Optional[dict] = cfg.get("seed")
        if self.seed is not None and "bit" not in self.seed:
            raise KeyError("purity: `seed` needs a `bit`")
        self.rate_lines: list[float] = [float(r) for r in cfg.get("rate_lines", [])]
        if self.rate_lines and not self.rate:
            raise ValueError("purity: `rate_lines` needs `rate`")
        self.rate_line_step: float = float(cfg.get("rate_line_step", 0.5))
        self.reference_bits: list[dict] = list(cfg.get("reference_bits") or [])
        for ref in self.reference_bits:
            if "bit" not in ref:
                raise KeyError("purity: each `reference_bits` entry needs a `bit`")
        norm = cfg.get("norm_events")
        self.norm_events: Optional[int] = int(norm) if norm is not None else None
        if self.norm_events is not None and "pu" in self.vs:
            raise ValueError("purity: `norm_events` does not work with the `pu` binning")

    @property
    def seed_label(self) -> str:
        if not self.seed:
            return ""
        return self.seed.get("label", self.seed["bit"])

    @property
    def gen_label(self) -> str:
        labels = [s.get("label", s["object"]) for s in self.reference_gen]
        return " or ".join(labels)

    def vs_edges(self, name: str) -> np.ndarray:
        return make_edges(self.vs[name])

    def vs_label(self, name: str) -> Optional[str]:
        spec = self.vs[name]
        return spec.get("label") if isinstance(spec, dict) else None


class PurityCollection(TurnOnCollection):
    """Trigger objects of each test object, their gen-match flag, and the derived
    quantities.

    After `create()`:
      objects[obj_key]: flat numpy arrays pt, eta, phi, matched (None in multiplicity
                        mode), event (event index), dr (matching cone)
      n_events, pileup (per event, or None)
      plots: list of dicts, one per output plot (see `_add_plot`)
    """

    def __init__(
        self,
        cfg_plot: dict,
        plot_name: Optional[str] = None,
        config_version: Optional[str] = None,
        override_version: Optional[str] = None,
    ):
        super().__init__(
            cfg_plot,
            threshold=-999,
            plot_name=plot_name,
            config_version=config_version,
            override_version=override_version,
        )
        self.pur_cfg = PurityConfig(
            self.cfg_plot.config_dict.get("purity"), self.cfg_plot.sample
        )
        self.objects: dict[str, dict[str, Any]] = {}
        self.n_events: int = 0
        self.pileup: Optional[np.ndarray] = None
        self.pileup_source: Optional[str] = None
        self.gen: Optional[ak.Array] = None
        self.plots: list[dict] = []
        # per-event seed bit and its rate (count, value, Garwood interval)
        self.seed_mask: Optional[np.ndarray] = None
        self.seed_rate: Optional[dict] = None
        # rates of the `reference_bits` (count, value, Garwood interval, label)
        self.reference_rates: list[dict] = []

    @property
    def matching(self) -> bool:
        return self.pur_cfg.mode != "multiplicity"

    @property
    def n_norm(self) -> int:
        """Events the per-event fractions and rates are normalised to: the sample's,
        or `norm_events` (the sample is a subset holding every event that can pass)."""
        return self.pur_cfg.norm_events or self.n_events

    # ---- loading -------------------------------------------------------------------

    @staticmethod
    def _cut_fields(cuts) -> set[str]:
        """Field names used in a list of cut strings (or a dict of such lists)."""
        if isinstance(cuts, dict):
            cuts = [c for v in cuts.values() for c in v]
        return {f for c in cuts or [] for f in re.findall(CUT_FIELD_RE, c)}

    def _load_parquet(self, name: str, fields: set[str]) -> ak.Array:
        """Only the needed columns of <version>_<sample>_<name>.parquet.

        As a record array with the field names of `ArrayLoader` (prefix stripped).
        """
        v = self.cfg_plot.version
        fname = f"{constants.CACHE_DIR}/{v}/{v}_{self.cfg_plot.sample}_{name}.parquet"
        if not os.path.exists(fname):
            raise FileNotFoundError(fname)
        loader = ArrayLoader(self, self.cfg_plot)
        cols = [
            c
            for c in ak.metadata_from_parquet(fname)["form"].fields
            if loader._transform_key(c, name) in fields
        ]
        try:
            arr = ak.from_parquet(fname, columns=cols)
        except Exception:
            # Caches written as one awkward record column cannot be pruned by pyarrow
            # ("only part of an extension's storage type"): read everything, keep the
            # fields.
            arr = ak.from_parquet(fname)[cols]
        return ak.zip({loader._transform_key(k, name): arr[k] for k in cols})

    def _load_arrays(self) -> None:
        for test_obj, _ in self.test_objects:
            fields = {"pt", "eta", "phi"} | self._cut_fields(test_obj.cuts)
            self.ak_arrays[str(test_obj)] = self._load_parquet(
                test_obj.nano_obj_name, fields
            )

    def _load_gen(self) -> None:
        """Selected gen objects of all `reference_gen` entries (eta, phi, pt),
        concatenated per event.
        """
        parts = []
        for spec in self.pur_cfg.reference_gen:
            key = spec["object"]
            name = key.split(":")[0]
            fields = {"pt", "eta", "phi"} | self._cut_fields(spec.get("cuts", []))
            if ":" in key:
                fields |= self._cut_fields(
                    Object(key, self.cfg_plot.version_for_objects).cuts
                )
            arr = self._load_parquet(name, fields)
            if ":" in key:
                # Object key with id: the cuts of the object definition
                obj = Object(key, self.cfg_plot.version_for_objects)
                arr = arr[compute_selection_mask_for_object_cuts(obj, arr)]
            for cut in spec.get("cuts", []):
                arr = arr[eval(re.sub(CUT_FIELD_RE, r"arr['\1']", cut))]
            parts.append(ak.zip({f: arr[f] for f in ("pt", "eta", "phi")}))
        self.gen = parts[0] if len(parts) == 1 else ak.concatenate(parts, axis=1)

    def _pileup_file(self) -> Optional[str]:
        version, sample = self.cfg_plot.version, self.cfg_plot.sample
        candidates = (
            [self.pur_cfg.pileup_version]
            if self.pur_cfg.pileup_version
            else [version, version + "PU"]
        )
        for v in candidates:
            fname = f"{constants.CACHE_DIR}/{v}/{v}_{sample}_Pileup.parquet"
            if os.path.exists(fname):
                return fname
        return None

    def _load_pileup(self, n_events: int) -> None:
        fname = self._pileup_file()
        if fname is None:
            print(
                f"WARNING: no Pileup parquet for {self.cfg_plot.version} / "
                f"{self.cfg_plot.sample}: skipping the pileup binning"
            )
            return
        arr = ak.from_parquet(fname)
        field = next(
            (
                f
                for f in arr.fields
                if f.removeprefix("Pileup_") == self.pur_cfg.pu_variable
            ),
            None,
        )
        if field is None:
            raise KeyError(
                f"{fname} has no field `{self.pur_cfg.pu_variable}` (fields: "
                f"{arr.fields})"
            )
        pileup = ak.to_numpy(arr[field]).astype(float)
        if len(pileup) != n_events:
            raise ValueError(
                f"{fname} has {len(pileup)} events, the test objects {n_events}: not "
                "row-aligned"
            )
        self.pileup = pileup
        self.pileup_source = fname

    def _load_bit(self, bit: str, record: Optional[str] = None) -> np.ndarray:
        """A per-event bit of <version>_<sample>_<record>.parquet, row-aligned with the
        test objects (record: default the bit's prefix)."""
        record = record or bit.split("_")[0]
        v, sample = self.cfg_plot.version, self.cfg_plot.sample
        fname = f"{constants.CACHE_DIR}/{v}/{v}_{sample}_{record}.parquet"
        if not os.path.exists(fname):
            raise FileNotFoundError(fname)
        arr = ak.from_parquet(fname)
        if bit not in arr.fields:
            raise KeyError(f"{fname} has no field `{bit}`")
        mask = ak.to_numpy(arr[bit]).astype(bool)
        if len(mask) != self.n_events:
            raise ValueError(
                f"{fname} has {len(mask)} events, the test objects {self.n_events}: "
                "not row-aligned"
            )
        return mask

    def _bit_rate(self, mask: np.ndarray) -> dict:
        """Count of a per-event bit, its rate (or fraction) and 68% Garwood interval."""
        n = int(mask.sum())
        scale = constants.RATE_NORM_FACTOR if self.pur_cfg.rate else 1
        scale /= self.n_norm
        lo, hi = garwood_interval(n)
        return {
            "n": n,
            "y": n * scale,
            "garwood_lo": lo * scale,
            "garwood_hi": hi * scale,
        }

    def _load_seed(self) -> None:
        """The per-event `seed` bit, row-aligned with the test objects."""
        seed = self.pur_cfg.seed
        if not seed:
            return
        self.seed_mask = self._load_bit(seed["bit"], seed.get("record"))
        self.seed_rate = self._bit_rate(self.seed_mask)

    def _load_reference_bits(self) -> None:
        """Rates of the `reference_bits` (e.g. the HLT paths an emulation compares to)."""
        for ref in self.pur_cfg.reference_bits:
            mask = self._load_bit(ref["bit"], ref.get("record"))
            self.reference_rates.append(
                dict(
                    self._bit_rate(mask),
                    bit=ref["bit"],
                    label=ref.get("label", ref["bit"]),
                )
            )

    def _add_seeded_objects(self) -> None:
        """`<key>+seed` copies of the objects listed in `seed.objects`, keeping only
        the objects of events in which the seed bit fired. Used for the event-level
        quantities only.
        """
        if self.seed_mask is None:
            return
        only = bool(self.pur_cfg.seed.get("only", False))
        keys = [
            str(Object(k, self.cfg_plot.version_for_objects))
            for k in self.pur_cfg.seed.get("objects") or []
        ] or list(self.objects)
        for key in keys:
            if key not in self.objects:
                raise KeyError(f"purity: seed object `{key}` is not a test object")
            obj = self.objects[key]
            sel = self.seed_mask[obj["event"]]
            seeded = {
                k: (v[sel] if isinstance(v, np.ndarray) else v) for k, v in obj.items()
            }
            seeded["base"] = key
            self.objects[f"{key}+seed"] = seeded
            if only:
                del self.objects[key]

    def _plain_objects(self) -> dict[str, dict[str, Any]]:
        """The test objects without the seeded copies."""
        return {k: o for k, o in self.objects.items() if "base" not in o}

    # ---- per-object flags ----------------------------------------------------------

    def _build_objects(self) -> None:
        lengths = {len(self.ak_arrays[str(o)]) for o, _ in self.test_objects}
        if len(lengths) != 1:
            raise ValueError(
                f"test objects of {self.cfg_plot.plot_name} have different event "
                f"counts {lengths}"
            )
        self.n_events = lengths.pop()
        if self.gen is not None and len(self.gen) != self.n_events:
            raise ValueError(
                f"gen objects have {len(self.gen)} events, the test objects "
                f"{self.n_events}"
            )
        for test_obj, _ in self.test_objects:
            key = str(test_obj)
            arr = self.ak_arrays[key]
            if arr.ndim < 2:
                raise ValueError(
                    f"{key}: purity needs a jagged object collection, not a per-event "
                    "quantity"
                )
            arr = ak.zip({f: arr[f] for f in ("pt", "eta", "phi")})
            if self.pur_cfg.leading_only:
                arr = arr[ak.argmax(arr.pt, axis=1, keepdims=True)]
                arr = arr[~ak.is_none(arr, axis=1)]
            dr = float(
                self.pur_cfg.match_dR
                if self.pur_cfg.match_dR is not None
                else test_obj.match_dR
            )
            matched = match_any(arr, self.gen, dr) if self.matching else None
            event = ak.broadcast_arrays(ak.local_index(arr, axis=0), arr.pt)[0]
            flat = lambda a: ak.to_numpy(ak.flatten(a, axis=None))  # noqa: E731
            self.objects[key] = {
                "pt": flat(arr.pt).astype(float),
                "eta": flat(arr.eta).astype(float),
                "phi": flat(arr.phi).astype(float),
                "event": flat(event).astype(np.int64),
                "matched": flat(matched).astype(bool) if matched is not None else None,
                "dr": dr,
            }
            del self.ak_arrays[key]

    # ---- quantities ----------------------------------------------------------------

    def _x(self, obj: dict, vs: str) -> np.ndarray:
        if vs == "abseta":
            return np.abs(obj["eta"])
        if vs == "pu":
            return self.pileup[obj["event"]]
        return obj[vs]

    def _add_plot(
        self,
        quantity: str,
        vs: str,
        threshold: Optional[float],
        series: dict,
        extra=None,
    ) -> None:
        edges = self.pur_cfg.vs_edges(vs)
        plot = {
            "quantity": quantity,
            "vs": vs,
            "threshold": threshold,
            "edges": edges,
            "series": series,
        }
        plot.update(extra or {})
        self.plots.append(plot)

    def _object_ratio(self, quantity: str, vs: str, threshold: Optional[float]) -> None:
        """purity / fake fraction of the trigger objects in bins of `vs`."""
        edges = self.pur_cfg.vs_edges(vs)
        series = {}
        for key, obj in self._plain_objects().items():
            sel = (
                obj["pt"] >= threshold
                if threshold is not None
                else np.ones(len(obj["pt"]), bool)
            )
            x = self._x(obj, vs)[sel]
            m = obj["matched"][sel]
            if quantity == "fake_fraction":
                m = ~m
            den = np.histogram(x, bins=edges)[0]
            num = np.histogram(x[m], bins=edges)[0]
            y, lo, hi = ratio_with_cp(num, den)
            total = ratio_with_cp([m.sum()], [len(m)])
            series[key] = {
                "y": y,
                "err_lo": lo,
                "err_hi": hi,
                "num": num,
                "den": den,
                "integral": float(total[0][0]) if len(m) else None,
                "integral_n": int(len(m)),
            }
        self._add_plot(quantity, vs, threshold, series)

    def _object_density(
        self, quantity: str, vs: str, threshold: Optional[float]
    ) -> None:
        """(Unmatched) objects per event in bins of an object variable,
        Poisson errors.
        """
        edges = self.pur_cfg.vs_edges(vs)
        series = {}
        for key, obj in self._plain_objects().items():
            sel = (
                obj["pt"] >= threshold
                if threshold is not None
                else np.ones(len(obj["pt"]), bool)
            )
            if quantity == "unmatched_density":
                sel &= ~obj["matched"]
            counts = np.histogram(self._x(obj, vs)[sel], bins=edges)[0]
            y = counts / self.n_norm
            err = np.sqrt(counts) / self.n_norm
            series[key] = {
                "y": y,
                "err_lo": err,
                "err_hi": err,
                "num": counts,
                "den": np.full(len(counts), self.n_norm),
            }
        self._add_plot(quantity, vs, threshold, series)

    def _per_event_counts(self, obj: dict, threshold: float) -> np.ndarray:
        sel = obj["pt"] >= threshold
        if self.matching:
            sel &= ~obj["matched"]
        return np.bincount(obj["event"][sel], minlength=self.n_events)

    def _event_quantity(self, counts: np.ndarray, quantity: str, ev_sel=None) -> tuple:
        """mean count per event, or fraction of events with >= `min_objects` (CP;
        x rate factor).
        """
        if ev_sel is not None:
            counts = counts[ev_sel]
        if quantity.endswith("per_event"):
            if ev_sel is None and self.n_norm > len(counts):
                # the events outside the sample have no objects above threshold
                counts = np.concatenate(
                    [counts, np.zeros(self.n_norm - len(counts), counts.dtype)]
                )
            m, e = mean_with_error(counts)
            return m, e, e, int(counts.sum()), len(counts)
        n_pass = int((counts >= self.pur_cfg.min_objects).sum())
        n_den = len(counts) if ev_sel is not None else self.n_norm
        y, lo, hi = ratio_with_cp([n_pass], [n_den])
        scale = constants.RATE_NORM_FACTOR if self.pur_cfg.rate else 1.0
        return y[0] * scale, lo[0] * scale, hi[0] * scale, n_pass, n_den

    def _event_vs_threshold(self, quantity: str) -> None:
        thresholds = self.pur_cfg.vs_edges("threshold")
        series = {}
        for key, obj in self.objects.items():
            cols = [
                self._event_quantity(self._per_event_counts(obj, t), quantity)
                for t in thresholds
            ]
            y, lo, hi, num, den = (np.asarray(c) for c in zip(*cols))
            series[key] = {"y": y, "err_lo": lo, "err_hi": hi, "num": num, "den": den}
            if quantity == "event_fraction" and self.pur_cfg.rate:
                series[key].update(self._garwood(num))
                if self.pur_cfg.rate_lines:
                    series[key]["at_rate"] = self._thresholds_at_rates(obj)
        self._add_plot(quantity, "threshold", None, series)

    def _garwood(self, num: np.ndarray) -> dict[str, np.ndarray]:
        """68% Garwood interval on the rate (kHz) from the event counts."""
        scale = constants.RATE_NORM_FACTOR / self.n_norm
        bounds = np.array([garwood_interval(int(n)) for n in num]).reshape(-1, 2)
        return {"garwood_lo": bounds[:, 0] * scale, "garwood_hi": bounds[:, 1] * scale}

    def _thresholds_at_rates(self, obj: dict) -> dict[str, Optional[dict]]:
        """Per rate line (kHz): the lowest threshold (on a `rate_line_step` grid over
        the threshold range) whose rate is at most the line, with the event count,
        the rate and its 68% Garwood interval. None if no threshold in the range does.
        """
        edges = self.pur_cfg.vs_edges("threshold")
        step = self.pur_cfg.rate_line_step
        grid = np.arange(edges[0], edges[-1] + step / 2, step)
        scale = constants.RATE_NORM_FACTOR / self.n_norm
        n_pass = np.array(
            [
                int((self._per_event_counts(obj, t) >= self.pur_cfg.min_objects).sum())
                for t in grid
            ]
        )
        out: dict[str, Optional[dict]] = {}
        for line in self.pur_cfg.rate_lines:
            ok = np.nonzero(n_pass * scale <= line)[0]
            if len(ok) == 0:
                out[f"{line:g}"] = None
                continue
            i = int(ok[0])
            lo, hi = garwood_interval(int(n_pass[i]))
            out[f"{line:g}"] = {
                "threshold": float(grid[i]),
                "n": int(n_pass[i]),
                "rate": float(n_pass[i] * scale),
                "garwood_lo": lo * scale,
                "garwood_hi": hi * scale,
                # the rate is already below the line at the start of the scan
                "at_scan_start": i == 0,
            }
        return out

    def _event_vs_pu(self, quantity: str, threshold: float) -> None:
        edges = self.pur_cfg.vs_edges("pu")
        assert self.pileup is not None
        ibin = np.digitize(self.pileup, edges) - 1
        centres = 0.5 * (edges[1:] + edges[:-1])
        series = {}
        for key, obj in self.objects.items():
            counts = self._per_event_counts(obj, threshold)
            cols = [
                self._event_quantity(counts, quantity, ibin == i)
                for i in range(len(edges) - 1)
            ]
            y, lo, hi, num, den = (np.asarray(c, dtype=float) for c in zip(*cols))
            fit = linear_fit(centres, y, 0.5 * (lo + hi))
            series[key] = {
                "y": y,
                "err_lo": lo,
                "err_hi": hi,
                "num": num,
                "den": den,
                "fit": fit,
            }
        self._add_plot(quantity, "pu", threshold, series)

    def _compute(self) -> None:
        mode = self.pur_cfg.mode
        vs_all = list(self.pur_cfg.vs)
        if "pu" in vs_all and self.pileup is None:
            vs_all.remove("pu")
        thresholds = self.pur_cfg.thresholds
        ratio_q: Optional[str]
        density_q: Optional[str]
        event_q: tuple[str, ...]
        if mode == "multiplicity":
            ratio_q, density_q, event_q = (
                None,
                "density",
                ("objects_per_event", "event_fraction"),
            )
        else:
            ratio_q = "purity" if mode == "purity" else "fake_fraction"
            density_q = "unmatched_density" if mode == "fake_rate" else None
            event_q = (
                ("unmatched_per_event", "event_fraction") if mode == "fake_rate" else ()
            )

        for vs in vs_all:
            if vs in OBJECT_VS:
                # vs pt: no threshold (the binning shows the pT dependence);
                # others: per threshold
                thrs: Sequence[Optional[float]] = thresholds
                if vs == "pt":
                    thrs = [None]
                for t in thrs:
                    if ratio_q:
                        self._object_ratio(ratio_q, vs, t)
                    if density_q and vs != "pt":
                        self._object_density(density_q, vs, t)
            elif vs == "pu":
                for t in thresholds:
                    if ratio_q:
                        self._object_ratio(ratio_q, vs, t)
                    for q in event_q:
                        self._event_vs_pu(q, t)
            elif vs == "threshold":
                for q in event_q:
                    self._event_vs_threshold(q)

    def create(self) -> None:
        self._load_arrays()
        self._apply_test_obj_cuts()
        if self.matching:
            self._load_gen()
        self._build_objects()
        self._load_seed()
        self._load_reference_bits()
        self._add_seeded_objects()
        self.gen = None
        self.ak_arrays = {}
        if "pu" in self.pur_cfg.vs:
            self._load_pileup(self.n_events)
        self._compute()
