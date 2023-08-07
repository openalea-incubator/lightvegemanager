from lightvegemanager.voxelsmesh import (
    compute_grid_size_from_trimesh,
    tesselate_trimesh_on_grid,
    fill_ratpgrid_from_trimesh,
    fill_ratpgrid_from_legumescene,
    reduce_layers_from_trimesh,
)
import pytest
import numpy
from alinea.pyratp import grid

pmax_1 = [15.0, 15.0, 15.0]
expected_1 = (15, 5, 16)
pmax_2 = [15.1, 15.1, 15.2]
expected_2 = (16, 6, 17)
pmax_3 = [0.0, 0.0, -1.0]
expected_3 = (2, 2, 2)
grid_slice = "ground = 0."
expected_4 = (15, 5, 15)


@pytest.mark.parametrize(
    "test_input, expected",
    [
        (pmax_1, expected_1),
        (pmax_2, expected_2),
        (pmax_3, expected_3),
        (grid_slice, expected_4),
    ],
)
def test_compute_grid_size_from_trimesh(test_input, expected):
    pmin = [0.0, 0.0, -1.0]
    pmax = [15.0, 15.0, 15.0]
    dv = [1.0, 3.0, 1.0]
    if type(test_input) != str:
        nx, ny, nz = compute_grid_size_from_trimesh(pmin, test_input, dv)
    else:
        nx, ny, nz = compute_grid_size_from_trimesh(pmin, pmax, dv, test_input)

    assert (nx, ny, nz) == expected


def test_tesselate_trimesh_on_grid():
    pars = {"latitude": 0, "longitude": 0, "timezone": 0, "nent": 1, "rs": (0, 0), "idecaly": 0, "orientation": 0}
    pars.update(
        {"njx": 2, "njy": 2, "njz": 2, "dx": 1.0, "dy": 1.0, "dz": [1.0] * 2, "xorig": 0, "yorig": 0, "zorig": 0}
    )
    mygrid = grid.Grid.initialise(**pars)
    scene = {
        0: [[(0.0, 0.0, 0.5), (0.0, 0.75, 0.5), (0.75, 0.75, 0.5)]],
        1: [[(0.0, 0.0, 0.5), (0.0, 1.5, 0.5), (1.5, 1.5, 0.5)]],
    }

    new_scene = tesselate_trimesh_on_grid(scene, mygrid, 3)
    assert len(new_scene[0]) == 1
    assert len(new_scene[1]) == 13


def test_fill_ratpgrid_from_trimesh():
    pars = {"latitude": 0, "longitude": 0, "timezone": 0, "nent": 2, "rs": (0, 0), "idecaly": 0, "orientation": 0}
    pars.update(
        {"njx": 2, "njy": 2, "njz": 2, "dx": 1.0, "dy": 1.0, "dz": [1.0] * 2, "xorig": 0, "yorig": 0, "zorig": 0}
    )
    mygrid = grid.Grid.initialise(**pars)
    scene = {
        0: [[(0.0, 0.0, 0.5), (0.0, 0.75, 0.5), (0.75, 0.75, 0.5)]],
        1: [[(1.1, 1.1, 1.5), (1.1, 1.85, 1.5), (1.85, 1.85, 1.5)]],
    }
    matching_ids = {
        0: [0, 0],
        1: [1, 1],
    }

    mygrid, match_tri_vox = fill_ratpgrid_from_trimesh(scene, matching_ids, mygrid)

    expected_s_vt = [0.28125, 0.28125]
    expected_leafareadensity = [[0.28125, 0.28125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0] * 8]
    expected_nume = [[1, 2, 0, 0, 0, 0, 0, 0], [0.0] * 8]

    numpy.testing.assert_array_equal(mygrid.s_vt, expected_s_vt)
    numpy.testing.assert_array_equal(mygrid.leafareadensity, expected_leafareadensity)
    numpy.testing.assert_array_equal(mygrid.nume, expected_nume)
    numpy.testing.assert_array_equal([mygrid.numx[0], mygrid.numy[0], mygrid.numz[0]], [1, 2, 2])
    numpy.testing.assert_array_equal([mygrid.numx[1], mygrid.numy[1], mygrid.numz[1]], [2, 1, 1])


def test_fill_ratpgrid_from_legumescene():
    pars = {"latitude": 0, "longitude": 0, "timezone": 0, "nent": 2, "rs": (0, 0), "idecaly": 0, "orientation": 0}
    pars.update(
        {"njx": 2, "njy": 2, "njz": 2, "dx": 2.0, "dy": 2.0, "dz": [1.0] * 2, "xorig": 0, "yorig": 0, "zorig": 0}
    )
    mygrid = grid.Grid.initialise(**pars)
    la = numpy.zeros([2, 3, 2, 2])
    la[0][2][0][0] = 0.7
    la[1][2][0][0] = 0.7
    legume_scene = {"LA": la}
    nb0 = 1

    fill_ratpgrid_from_legumescene(legume_scene, mygrid, nb0)
    expected_s_vt = [0.7 + 7e-14, 0.7 + 7e-14]
    expected_leafareadensity = [
        [1e-14 / 4, 0.7 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0],
        [1e-14 / 4.0, 0.7 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0, 1e-14 / 4.0],
    ]
    expected_nume = [[1] * 8, [2] * 8]

    # conversion from float32 to float64 involves mantis errors (around 1e-8 difference)
    numpy.testing.assert_allclose(mygrid.s_vt, expected_s_vt)
    numpy.testing.assert_allclose(mygrid.leafareadensity, expected_leafareadensity)
    numpy.testing.assert_array_equal(mygrid.nume, expected_nume)
    numpy.testing.assert_array_equal([mygrid.numx[1], mygrid.numy[1], mygrid.numz[1]], [2, 1, 2])


@pytest.mark.parametrize("test_input, expected", [(None, 3), ([0, 1], 2)])
def test_reduce_layers_from_trimesh(test_input, expected):
    scene = {
        0: [[(0.0, 0.0, 0.5), (0.0, 0.75, 0.5), (0.75, 0.75, 0.5)]],
        1: [[(1.1, 1.1, 1.5), (1.1, 1.85, 1.5), (1.85, 1.85, 1.5)]],
    }
    pmax = (1.0, 1.0, 1.0)
    matching_ids = {0: [0, 0], 1: [1, 1]}
    nxyz = [1, 1, 6]
    dxyz = [0.5, 0.5, 0.5]

    skylayer = reduce_layers_from_trimesh(scene, pmax, dxyz, nxyz, matching_ids, test_input)

    assert skylayer == expected
