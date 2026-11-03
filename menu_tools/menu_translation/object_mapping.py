"""Mapping between MenuTools object keys and CMSSW P2GT object definitions.

The cut languages of the two frameworks differ (e.g. `{hwQual} >= 3` vs the
bit mask `qualityFlags`), so object cuts are never translated: every object
used in a translated seed must be listed in `configs/<version>/cmssw_object_map.yaml`:

    "L1gmtTkMuon:VLoose":
      cmssw_object: l1tGTtkMuonVLoose
    "L1puppiJetSC4sums:MHT":
      cmssw_object: l1tGTHtSum
      ptScaling: l1tGTScaling_CL2HtSum_MHT  # only if not the object default
      threshold: pt                          # pt (default) or scalarSumPt
    primary_vertex: "L1PV:default"
"""

import glob
import os
from dataclasses import dataclass
from typing import Optional

import yaml

from menu_tools.menu_translation.seed_model import TranslationError


@dataclass(frozen=True)
class CmsswObject:
    object: str
    ptScaling: Optional[str] = None
    threshold: str = "pt"  # "pt" or "scalarSumPt"


class ObjectMapping:
    def __init__(self, mapping: dict[str, CmsswObject], primary_vertex: str):
        self.to_cmssw_map = mapping
        self.primary_vertex = primary_vertex
        self._from_cmssw: dict[tuple, str] = {}
        for key, obj in mapping.items():
            if obj in self._from_cmssw:
                raise TranslationError(
                    f"Ambiguous object mapping: {key} and {self._from_cmssw[obj]}"
                )
            self._from_cmssw[obj] = key

    @classmethod
    def from_file(cls, path: str) -> "ObjectMapping":
        with open(path) as f:
            raw = yaml.safe_load(f)
        primary_vertex = raw.pop("primary_vertex")
        mapping = {}
        for key, value in raw.items():
            if isinstance(value, str):
                value = {"cmssw_object": value}
            mapping[key] = CmsswObject(
                object=value["cmssw_object"],
                ptScaling=value.get("ptScaling"),
                threshold=value.get("threshold", "pt"),
            )
            if mapping[key].threshold not in ("pt", "scalarSumPt"):
                raise TranslationError(f"{key}: unknown threshold type")
        return cls(mapping, primary_vertex)

    def to_cmssw(self, key: str) -> CmsswObject:
        try:
            return self.to_cmssw_map[key]
        except KeyError:
            raise TranslationError(
                f"no CMSSW object mapped for MenuTools object '{key}'"
            )

    def from_cmssw(self, obj: CmsswObject) -> str:
        try:
            return self._from_cmssw[obj]
        except KeyError:
            raise TranslationError(f"no MenuTools object mapped for CMSSW {obj}")

    def undefined_menu_tools_keys(self, objects_dir: str) -> list[str]:
        """Mapped keys whose object/ID is not defined in `objects_dir`/*.y*ml"""
        defined: dict[str, set] = {}
        for path in glob.glob(os.path.join(objects_dir, "*.y*ml")):
            with open(path) as f:
                for name, cfg in (yaml.safe_load(f) or {}).items():
                    defined.setdefault(name, set()).update(
                        (cfg.get("ids") or {}).keys()
                    )
        keys = list(self.to_cmssw_map) + [self.primary_vertex]
        undefined = []
        for key in keys:
            name, obj_id = key.split(":")[:2]
            if obj_id not in defined.get(name, set()):
                undefined.append(key)
        return undefined
