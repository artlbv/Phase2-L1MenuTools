import os
import textwrap

import pytest
import yaml

from menu_tools.menu_translation.cmssw_menu_writer import seed_to_cmssw
from menu_tools.menu_translation.cross_mask_parser import (
    LegCut,
    PairCut,
    parse_cross_mask,
)
from menu_tools.menu_translation.object_mapping import CmsswObject, ObjectMapping
from menu_tools.menu_translation.rate_table_yaml_io import (
    parse_threshold,
    read_rate_table_menu,
    seed_from_yaml,
    seed_to_yaml,
)
from menu_tools.menu_translation.seed_model import (
    Leg,
    Seed,
    Threshold,
    TranslationError,
)

PV = "L1PV:default"
MAPPING = ObjectMapping(
    {
        "L1gmtTkMuon:VLoose": CmsswObject("l1tGTtkMuonVLoose"),
        "L1gmtTkMuon:Loose": CmsswObject("l1tGTtkMuonLoose"),
        "L1puppiJetSC4sums:HT": CmsswObject("l1tGTHtSum", threshold="scalarSumPt"),
        "L1puppiJetSC4sums:MHT": CmsswObject(
            "l1tGTHtSum", ptScaling="l1tGTScaling_CL2HtSum_MHT"
        ),
    },
    PV,
)

DIMUON_YAML = {
    "cross_masks": [
        "((abs(leg2.z0-leg1.z0) < 1))",
        "((leg2.charge*leg3.charge < 0.0) & (leg2.deltaR(leg3) > 0))",
        "(leg2+leg3).mass < 9.0",
        "abs(leg3.eta) < 1.5",
    ],
    "leg1": {"threshold_cut": None, "obj": PV},
    "leg2": {"threshold_cut": "offline_pt >= 15.0", "obj": "L1gmtTkMuon:VLoose"},
    "leg3": {"threshold_cut": "pt > 7", "obj": "L1gmtTkMuon:Loose"},
}


def test_parse_cross_mask():
    cuts, warnings = parse_cross_mask("(leg1.deltaR(leg2) < 1.2) & (abs(leg1.eta) < 2)")
    assert cuts == [PairCut(1, 2, "dR", "<", 1.2), LegCut(1, "absEta", "<", 2.0)]
    assert warnings == []


def test_parse_cross_mask_precedence_warning():
    cuts, warnings = parse_cross_mask(
        "(abs(leg2.z0-leg1.z0) < 1 & (leg2.deltaR(leg3) > 0))"
    )
    assert cuts == [PairCut(2, 1, "dz", "<", 1.0), PairCut(2, 3, "dR", ">", 0.0)]
    assert len(warnings) == 1


def test_parse_cross_mask_unsupported():
    with pytest.raises(TranslationError):
        parse_cross_mask("(leg2.btagScore + leg3.btagScore) > 2.2")
    with pytest.raises(TranslationError):
        parse_cross_mask("leg1.pt / leg2.pt > 2")


def test_parse_threshold():
    assert parse_threshold("offline_pt >= 22.0") == Threshold(22.0, True, True)
    assert parse_threshold("leg1.offline_pt > 28") == Threshold(28.0, True, False)
    assert parse_threshold("pt > 4.4") == Threshold(4.4, False, False)
    assert parse_threshold(None) is None
    with pytest.raises(TranslationError):
        parse_threshold("offline_pt < 10")


def test_seed_from_yaml():
    seed, _ = seed_from_yaml("L1_DoubleTkMu", DIMUON_YAML, PV)
    assert seed.legs == [
        Leg("L1gmtTkMuon:VLoose", Threshold(15.0), pv_dz_max=1.0),
        Leg("L1gmtTkMuon:Loose", Threshold(7.0, False, False), abs_eta_max=1.5),
    ]
    assert seed.correlations == {(1, 2): {"os": True, "minDR": 0.0, "maxInvMass": 9.0}}


def test_yaml_round_trip():
    seed, _ = seed_from_yaml("L1_DoubleTkMu", DIMUON_YAML, PV)
    again, _ = seed_from_yaml("L1_DoubleTkMu", seed_to_yaml(seed, PV), PV)
    assert again == seed


def test_conflicting_cuts():
    config = dict(
        DIMUON_YAML, cross_masks=["(leg2+leg3).mass < 9", "(leg2+leg3).mass < 10"]
    )
    with pytest.raises(TranslationError):
        seed_from_yaml("L1_X", config, PV)


def test_read_rate_table_menu_reports_skipped(tmp_path):
    menu = {
        "L1_Good": DIMUON_YAML,
        "L1_Bad": {"cross_masks": ["leg1.foo > 1"], "leg1": DIMUON_YAML["leg2"]},
    }
    path = tmp_path / "menu.yml"
    path.write_text(yaml.safe_dump(menu))
    seeds, skipped, _ = read_rate_table_menu(str(path), PV)
    assert [s.name for s in seeds] == ["L1_Good"]
    assert list(skipped) == ["L1_Bad"]


def test_seed_to_cmssw():
    seed, _ = seed_from_yaml("L1_DoubleTkMu", DIMUON_YAML, PV)
    text = seed_to_cmssw(seed, MAPPING)
    assert text.startswith("DoubleTkMu = l1tGTDoubleObjectCond.clone(")
    for fragment in [
        'object = gt_ref("l1tGTtkMuonVLoose"),',
        "offlineMinPt = cms.double(15.0),",
        "maxPrimVertDz = cms.double(1.0),",
        "primVertex = cms.uint32(0),",
        "minPt = cms.double(7.0),",
        "minEta = cms.double(-1.5),",
        "os = cms.bool(True),",
        "maxInvMass = cms.double(9.0),",
        "pDoubleTkMu = cms.Path(DoubleTkMu)",
    ]:
        assert fragment in text


def test_seed_to_cmssw_sums():
    ht = Seed("L1_HT", [Leg("L1puppiJetSC4sums:HT", Threshold(450.0))])
    mht = Seed("L1_MHT", [Leg("L1puppiJetSC4sums:MHT", Threshold(135.5))])
    assert "offlineMinScalarSumPt = cms.double(450.0)" in seed_to_cmssw(ht, MAPPING)
    text = seed_to_cmssw(mht, MAPPING)
    assert 'ptScaling = gt_ref("l1tGTScaling_CL2HtSum_MHT")' in text
    assert "offlineMinPt = cms.double(135.5)" in text


def test_unmapped_object():
    seed = Seed("L1_X", [Leg("L1gmtTkMuon:Medium", Threshold(22.0))])
    with pytest.raises(TranslationError):
        seed_to_cmssw(seed, MAPPING)


def test_config_object_map():
    path = os.path.join(
        os.path.dirname(__file__), "../../../configs/V43nano/cmssw_object_map.yaml"
    )
    mapping = ObjectMapping.from_file(path)
    assert mapping.to_cmssw("L1gmtTkMuon:VLoose").object == "l1tGTtkMuonVLoose"
    assert mapping.primary_vertex == PV


# Reading CMSSW configurations needs FWCore.ParameterSet from a CMSSW checkout
CMSSW_SRC = os.environ.get("CMSSW_SRC")

DUMP = textwrap.dedent(
    """
    import FWCore.ParameterSet.Config as cms
    process = cms.Process("DUMP")
    process.l1tGTScaling_GMTTkMuons_VLoose = cms.PSet(
        regionsOffset = cms.vdouble(1, 1, 1), regionsSlope = cms.vdouble(1, 1, 1))
    process.l1tGTtkMuonVLoose = cms.PSet(
        tag = cms.InputTag("l1tGTProducer", "GMTTkMuons"),
        ptScaling = cms.PSet(refToPSet_ = cms.string("l1tGTScaling_GMTTkMuons_VLoose")))
    process.l1tGTtkMuonLoose = process.l1tGTtkMuonVLoose.clone()
    process.DoubleTkMu = cms.EDFilter("L1GTDoubleObjectCond",
        collection1 = cms.PSet(
            object = cms.PSet(refToPSet_ = cms.string("l1tGTtkMuonVLoose")),
            offlineMinPt = cms.double(15), maxPrimVertDz = cms.double(1),
            primVertex = cms.uint32(0)),
        collection2 = cms.PSet(
            object = cms.PSet(refToPSet_ = cms.string("l1tGTtkMuonLoose")),
            regionsMinPt = cms.vdouble(7, 7, 7), minEta = cms.double(-1.5),
            maxEta = cms.double(1.5)),
        maxInvMass = cms.double(9), minDR = cms.double(0), os = cms.bool(True),
        primVertTag = cms.InputTag("l1tGTProducer", "GTTPrimaryVert"))
    process.pDoubleTkMu = cms.Path(process.DoubleTkMu)
    process.l1tGTAlgoBlockProducer = cms.EDProducer("L1GTAlgoBlockProducer",
        algorithms = cms.VPSet(
            cms.PSet(expression = cms.string("pDoubleTkMu")),
            cms.PSet(expression = cms.string("pDoubleTkMu and pDoubleTkMu"))))
    """
)


@pytest.mark.skipif(not CMSSW_SRC, reason="set CMSSW_SRC to a CMSSW source checkout")
def test_read_cmssw_dump(tmp_path):
    from menu_tools.menu_translation.cmssw_menu_reader import read_cmssw_menu

    path = tmp_path / "dump.py"
    path.write_text(DUMP)
    seeds, skipped = read_cmssw_menu(str(path), MAPPING, CMSSW_SRC)
    expected, _ = seed_from_yaml("L1_DoubleTkMu", DIMUON_YAML, PV)
    assert seeds == [expected]
    assert list(skipped) == ["pDoubleTkMu and pDoubleTkMu"]
