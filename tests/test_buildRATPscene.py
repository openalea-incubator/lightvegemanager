from lightvegemanager.buildRATPscene import (
    extract_grid_origin,
    build_RATPscene_empty,
    build_RATPscene_from_trimesh,
    legumescene_to_RATPscene,
)
import pytest
import numpy
import os


parameters_empty = {}
parameters_orig2 = {"origin": (0, 0)}
parameters_orig3 = {"origin": (0, 0, 1)}
minmax = [(-1, -1, -1), (5, 5, 5)]


@pytest.mark.parametrize(
    "test_input, expected",
    [(parameters_empty, (-1, -1, 1)), (parameters_orig2, (0, 0, 1)), (parameters_orig3, (0, 0, 1))],
)
def test_extract_grid_origin(test_input, expected):
    result = extract_grid_origin(test_input, minmax)
    assert result == expected


parameters_empty = {"voxel size": (1.0, 1.0, 1.0)}
parameters_nvox = {"voxel size": (1.0, 1.0, 1.0), "number voxels": (5, 5, 5)}
coordinates = [0.0, 0.0, 0.0]


@pytest.mark.parametrize("test_input,expected", [(parameters_empty, (1, 1, 1)), (parameters_nvox, (5, 5, 5))])
def test_build_RATPscene_empty(test_input, expected):
    grid, distrib = build_RATPscene_empty(test_input, minmax, coordinates, True)

    assert distrib == {"global": [[1.0]]}
    assert (grid.dx, grid.dy, grid.dz[0], grid.njx, grid.njy, grid.njz, grid.isolated_box) == (
        test_input["voxel size"][0],
        test_input["voxel size"][1],
        test_input["voxel size"][2],
        expected[0],
        expected[1],
        expected[2],
        0,
    )


def test_legumescene_to_RATPscene():
    la = numpy.zeros([1, 10, 4, 4])
    for i in range(4):
        for j in range(4):
            la[0][9][i][j] = 2.0
    legumescene = {"LA": la, "distrib": numpy.array([[0.2, 0.2, 0.4, 0.1, 0.1]])}
    parameters = {"voxel size": (1.0, 1.0, 1.0)}

    grid, distrib, nb0 = legumescene_to_RATPscene(legumescene, parameters, coordinates, False, True)

    assert numpy.array_equal(distrib["global"], numpy.array([[0.2, 0.2, 0.4, 0.1, 0.1]]))
    assert nb0 == 9
    assert (grid.dx, grid.dy, grid.dz[0], grid.njx, grid.njy, grid.njz, grid.isolated_box, grid.s_canopy) == (
        parameters["voxel size"][0],
        parameters["voxel size"][1],
        parameters["voxel size"][2],
        4,
        4,
        1,
        0,
        32.0,
    )


parameters_1 = {
    "voxel size": (0.5, 0.5, 0.5),
    "angle distrib algo": "compute global",
    "nb angle classes": 12,
    "grid slicing": "ground = 0.",
    "tesselation level": 3,
    "soil reflectance": [0.5, 0.0],
}
more_inputs_1 = (None, 2)
expected_1 = {}
expected_1["nxyz"] = (7, 2, 2)
expected_1["nveg"] = 3
expected_1["nent"] = 2
expected_1["area"] = 0.47874045
expected_1["orig"] = (0.2, 0.1, -0.1)
expected_1["nb classes"] = 12
expected_1["nb triangles"] = 48
expected_1["rs"] = [0.0, 0.0]

parameters_2 = {
    "voxel size": (1.0, 1.0, 1.0),
    "angle distrib algo": "compute voxel",
    "nb angle classes": 18,
    "tesselation level": 0,
}
more_inputs_2 = ([(222, 0)], 2)
expected_2 = {}
expected_2["nxyz"] = (4, 1, 1)
expected_2["nveg"] = 3
expected_2["nent"] = 3
expected_2["area"] = 0.40063795
expected_2["orig"] = (0.2, 0.1, -0.1)
expected_2["nb classes"] = 18
expected_2["nb triangles"] = 3
expected_2["rs"] = [0.0, 0.0]

parameters_3 = {
    "voxel size": (1.0, 1.0, 1.0),
    "origin": (0.0, 0.0, 0.0),
    "number voxels": (4, 4, 4),
    "angle distrib algo": "file",
    "angle distrib file": os.path.join(
        os.path.dirname(os.path.dirname(__file__)), "data", "luzerne_angle_distrib.data"
    ),
    "tesselation level": 0,
}
more_inputs_3 = (None, 2)
expected_3 = {}
expected_3["nxyz"] = (4, 4, 4)
expected_3["nveg"] = 3
expected_3["nent"] = 1
expected_3["area"] = 0.47874045
expected_3["orig"] = (0.0, 0.0, 0.0)
expected_3["nb classes"] = 9
expected_3["nb triangles"] = 3
expected_3["rs"] = [0.0, 0.0]

parameters_4 = {
    "voxel size": "dynamic",
    "angle distrib algo": "compute global",
    "nb angle classes": 9,
    "soil reflectance": [0.5, 0.0],
    "tesselation level": 0,
}
more_inputs_4 = (None, 2)
expected_4 = {}
expected_4["nxyz"] = (2, 1, 1)
expected_4["nveg"] = 2
expected_4["nent"] = 2
expected_4["area"] = 0.47874045
expected_4["orig"] = (0.2, 0.1, -0.1)
expected_4["nb classes"] = 9
expected_4["nb triangles"] = 3
expected_4["rs"] = [0.5, 0.0]


@pytest.mark.parametrize(
    "test_input_1, test_input_2, test_input_3, expected",
    [
        (parameters_1, False, more_inputs_1, expected_1),
        (parameters_2, False, more_inputs_2, expected_2),
        (parameters_3, False, more_inputs_3, expected_3),
        (parameters_4, True, more_inputs_4, expected_4),
    ],
)
def test_build_RATPscene_from_trimesh(test_input_1, test_input_2, test_input_3, expected):
    """On vérifie que les options ont bien fonctionné"""
    triangles = {
        0: [[(0.2, 0.1, 0.1), (0.7, 0.1, 0.1), (0.7, 0.7, 0.6)]],
        1: [[(1.3, 0.1, 0.1), (1.7, 0.1, 0.1), (1.7, 0.7, 0.6)]],
        2: [[(3.4, 0.1, 0.1), (3.7, 0.1, 0.1), (3.7, 0.7, 0.7)]],
    }
    minmax = [(0.2, 0.1, 0.1), (3.7, 0.7, 0.5)]
    matching_ids = {
        0: (123, 0),
        1: (222, 0),
        2: (111, 1),
    }
    if "angle distrib file" in test_input_1:
        matching_ids[2] = (111, 0)

    if isinstance(test_input_3[0], list):
        matching_ids = {
            0: (123, 0),
            1: (222, 2),
            2: (111, 1),
        }
    infinite = True
    triLmax = 0.5

    ratpgrid, matching_tri_vox, distrib = build_RATPscene_from_trimesh(
        triangles,
        minmax,
        triLmax,
        matching_ids,
        test_input_1,
        coordinates,
        test_input_2,
        infinite,
        stems_id=test_input_3[0],
        nb_input_scenes=test_input_3[1],
    )

    # vérifie nb vox créé
    assert (ratpgrid.njx, ratpgrid.njy, ratpgrid.njz) == expected["nxyz"]

    # vérifie nb vox plein
    assert ratpgrid.nveg == expected["nveg"]

    # vérifie nb ent (prise en charge de stiges)
    assert ratpgrid.nent == expected["nent"]
    assert pytest.approx(ratpgrid.s_canopy) == expected["area"]

    # vérifie grid origin
    numpy.testing.assert_allclose(
        (float(ratpgrid.xorig), float(ratpgrid.yorig), float(ratpgrid.zorig)), expected["orig"]
    )

    # vérifie nb de class
    assert len(distrib["global"][0]) == expected["nb classes"]
    assert len(distrib["global"]) == ratpgrid.nent

    # vérifie si distrib voxel a rempli tous les vox
    if "voxel" in distrib:
        assert len(distrib["voxel"][0][0]) == expected["nb classes"]
        assert len(distrib["voxel"]) == ratpgrid.nveg

    # vérifie si la tesselation a eu lieu
    assert len(matching_tri_vox) == expected["nb triangles"]

    # prise en charge du réfléchi
    numpy.testing.assert_array_equal(ratpgrid.rs, expected["rs"])
