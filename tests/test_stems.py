from lightvegemanager.stems import *
from openalea.mtg.mtg import MTG
import openalea.plantgl.all as pgl_all

def test_extract_stems_from_MTG():
    # empty MTG with a stem property
    my_mtg = MTG()
    geom = pgl_all.FaceSet([(0.1,0.1,0.96),      \
                                (0.1,0.6,0.96),  \
                                (0.6,0.6,0.96),  \
                                (0.6,0.1,0.96)], \
                                [range(4)])
    scene = pgl_all.Scene([
                            pgl_all.Shape(geom, pgl_all.Material((250,0,0),1), 19)
                                                    ])
    
    geo_property = {"geometry" : scene, "label" : "StemElement"}                                          
    
    my_mtg._add_vertex_properties(19, geo_property)

    entity_id = 0

    stems = extract_stems_from_MTG(my_mtg, entity_id)

    assert stems == [(19,0)]

def test_manage_stems_for_ratp():
    
    stems_id = [ (19,0), (20,0) ]
    matching_ids = {
        0 : [18, 0],
        1 : [19, 0],
        2 : [20, 0]
    }
    ratp_parameters = {
        "reflectance coefficients" : [ [0.5] ],
        "mu" : [0.5]
    }

    manage_stems_for_ratp(stems_id, matching_ids, ratp_parameters)
    
    assert matching_ids == {
                                0 : [18, 0],
                                1 : [19, 1],
                                2 : [20, 1]
                            }

    assert ratp_parameters == {
                                "reflectance coefficients" : [ [0.5], [0.5] ],
                                "mu" : [0.5, 0.5]
                            }
