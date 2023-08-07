from lightvegemanager.trianglesmesh import (
    triangles_entity,
    globalid_to_elementid,
    globalid_to_triangle,
    compute_area_max,
    compute_minmax_coord,
    compute_trilenght_max,
    chain_triangulations,
    vgx_to_caribu,
    apply_transformations,
    create_heterogeneous_canopy,
)
import pytest
import numpy
import os
import itertools

import openalea.plantgl.all as pgl
from alinea.adel.adel_dynamic import AdelDyn
from alinea.adel.echap_leaf import echap_leaves


def test_triangles_entity():
    scene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        1: [[(1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)], [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]],
        2: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        3: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
    }
    matching_ids = {0: [11, 0], 1: [22, 1], 2: [33, 0], 3: [41, 0]}
    entity_id = 1

    triangles = triangles_entity(scene, entity_id, matching_ids)

    assert triangles[0][0][0] == 1.0


def test_globalid_to_elementid():
    scene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        1: [
            [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)],
        ],
        2: [
            [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)],
        ],
        3: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
    }
    triangles_id = 5
    id = globalid_to_elementid(scene, triangles_id)
    assert id == 3

    triangles_id = 100
    with pytest.raises(Exception) as e_info:
        id = globalid_to_elementid(scene, triangles_id)
    assert str(e_info.typename) == "IndexError"


def test_globalid_triangle():
    scene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        1: [
            [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)],
        ],
        2: [
            [(5.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
        ],
        3: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
    }

    triangles_id = 100
    with pytest.raises(Exception) as e_info:
        t = globalid_to_triangle(scene, triangles_id)
    assert str(e_info.typename) == "IndexError"

    triangles_id = 5
    t = globalid_to_triangle(scene, triangles_id)
    assert t[0][0] == 5.0


def test_compute_area_max():
    scene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        1: [[(0.0, 0.0, 0.0), (0.0, 2.0, 0.0), (2.0, 2.0, 0.0)]],
        2: [[(0.0, 0.0, 0.0), (0.0, 4.0, 0.0), (4.0, 4.0, 0.0)]],
        3: [[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)]],
    }

    a = compute_area_max(scene)
    assert a == 8.0


def test_compute_minmax_coord():
    scene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        1: [[(0.0, 0.0, 0.0), (0.0, 2.0, 0.0), (2.0, 2.0, 0.0)]],
        2: [[(0.0, 0.0, 0.0), (0.0, 4.0, 0.0), (4.0, 4.0, 0.0)]],
        3: [[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)]],
    }

    pmin, pmax = compute_minmax_coord(scene)
    numpy.testing.assert_array_equal(pmin, [0.0, 0.0, 0.0])
    numpy.testing.assert_array_equal(pmax, [4.0, 4.0, 0.0])


def test_compute_trilenght_max():
    scene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        1: [[(0.0, 0.0, 0.0), (0.0, 2.0, 0.0), (2.0, 2.0, 0.0)]],
        2: [[(0.0, 0.0, 0.0), (0.0, 4.0, 0.0), (4.0, 4.0, 0.0)]],
        3: [[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)]],
    }

    lmax = compute_trilenght_max(scene)
    assert lmax == 4.0


def test_vgx_to_caribu():
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "NICatObs1P2.vgx")
    scene = vgx_to_caribu(path, 1)

    assert scene[1][1][0][2] == 95.0273


def test_chain_triangulations():
    # plantGL scene
    pgl_scene = pgl.Scene([pgl.Shape(pgl.Box(), pgl.Material(), 888), pgl.Shape(pgl.Cylinder(), pgl.Material(), 999)])

    # fichier VGX
    vgx_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "NICatObs1P2.vgx")

    # scene CARIBU
    cscene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        11: [[(1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]],
        24: [
            [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)],
        ],
        3: [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
    }

    # table MTG with plantGL scene
    INPUTS_DIRPATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    adel_wheat = AdelDyn(seed=1, scene_unit="m", leaves=echap_leaves(xy_model="Soissons_byleafclass"))
    g = adel_wheat.load(dir=INPUTS_DIRPATH)

    # scene l-egume
    l_scene = {"LA": numpy.ones([2, 10, 4, 4]), "distrib": [[0.5, 0.5], [0.3, 0.7]]}

    scenes = [pgl_scene, vgx_path, cscene, g, l_scene]

    complete_trimesh, matching_ids, legume_grid, id_legume_scene = chain_triangulations(scenes)

    assert legume_grid
    assert id_legume_scene == 4

    expected_matching_ids = {
        # plantGL
        0: [888, 0],
        1: [999, 0],
        # vgx file
        2: [2, 1],
        # dict scene
        3: [0, 2],
        4: [11, 2],
        5: [24, 2],
        6: [3, 2],
        # MTG from an adelwheat example
        7: [19, 3],
        8: [34, 3],
        9: [813, 3],
        10: [814, 3],
        11: [51, 3],
    }

    assert expected_matching_ids == matching_ids
    # verify the mesh object
    assert len(complete_trimesh) == len(expected_matching_ids)
    assert numpy.cumsum([len(v) for v in complete_trimesh.values()])[-1] == 5564
    for t in itertools.chain(*complete_trimesh.values()):
        assert len(t) == 3
        for i in range(3):
            assert len(t[i]) == 3


transformation_1 = {"scenes unit": {0: "m", 1: "dm"}}
expected_1 = {
    0: [[(0.0, 0.0, 0.0), (0.0, 100.0, 0.0), (100.0, 100.0, 100.0)]],
    1: [[(0.0, 0.0, 10.0), (0.0, 10.0, 10.0), (10.0, 10.0, 20.0)]],
}
transformation_2 = {"rescale": {0: 1000}}
expected_2 = {
    0: [[(0.0, 0.0, 0.0), (0.0, 1000.0, 0.0), (1000.0, 1000.0, 1000.0)]],
    1: [[(0.0, 0.0, 1.0), (0.0, 1.0, 1.0), (1.0, 1.0, 2.0)]],
}
transformation_3 = {"translate": {1: (5.0, 5.0, 5.0)}}
expected_3 = {
    0: [[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 1.0)]],
    1: [[(5.0, 5.0, 6.0), (5.0, 6.0, 6.0), (6.0, 6.0, 7.0)]],
}
transformation_4 = {"xyz orientation": {0: "x+ = S", 1: "x+ = W"}}
expected_4 = {
    0: [[(0.0, 0.0, 0.0), (0.0, -1.0, 0.0), (-1.0, -1.0, 1.0)]],
    1: [[(0.0, 0.0, 1.0), (-1.0, 0.0, 1.0), (-1.0, 1.0, 2.0)]],
}
transformation_5 = {
    "scenes unit": {0: "m", 1: "dm"},
    "xyz orientation": {0: "x+ = S", 1: "x+ = W"},
    "rescale": {0: 5000},
    "translate": {1: (5.0, 5.0, 5.0)},
}
expected_5 = {
    0: [[(0.0, 0.0, 0.0), (0.0, -500000.0, 0.0), (-500000.0, -500000.0, 500000.0)]],
    1: [[(-5.0, 5.0, 15.0), (-15.0, 5.0, 15.0), (-15.0, 15.0, 25.0)]],
}


@pytest.mark.parametrize(
    "test_input, expected",
    [
        (transformation_1, expected_1),
        (transformation_2, expected_2),
        (transformation_3, expected_3),
        (transformation_4, expected_4),
        (transformation_5, expected_5),
    ],
)
def test_apply_transformations(test_input, expected):
    global_scene_unit = "cm"
    scene = {
        0: [[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 1.0)]],
        1: [[(0.0, 0.0, 1.0), (0.0, 1.0, 1.0), (1.0, 1.0, 2.0)]],
    }
    matching_ids = {0: [0, 0], 1: [1, 1]}

    apply_transformations(scene, matching_ids, test_input, global_scene_unit)

    numpy.testing.assert_allclose(scene[0], expected[0], atol=1e-10)
    numpy.testing.assert_allclose(scene[1], expected[1], atol=1e-10)
