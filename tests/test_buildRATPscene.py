import pytest
from lightvegemanager.buildRATPscene import *


parameters_empty = {}
parameters_orig2 = { "origin" : (0,0)}
parameters_orig3 = { "origin" : (0,0,1)}
minmax = [(-1, -1, -1), (5, 5, 5)]
@pytest.mark.parametrize("test_input,expected", [
                                            (parameters_empty, (-1, -1, 1)), 
                                            (parameters_orig2, (0, 0, 1)), 
                                            (parameters_orig3, (0, 0, 1))])
def test_extract_grid_origin(test_input, expected):
    result = extract_grid_origin(test_input, minmax)
    assert result == expected

parameters_empty = {"voxel size" : (1., 1., 1.)}
parameters_nvox = {"voxel size" : (1., 1., 1.), "number voxels" : (5, 5, 5)}
coordinates = [0., 0., 0.]
@pytest.mark.parametrize("test_input,expected", [
                                            (parameters_empty, (1,1,1)), 
                                            (parameters_nvox, (5, 5, 5))])
def test_build_RATPscene_empty(test_input, expected):
    grid, distrib = build_RATPscene_empty(test_input, minmax, coordinates, True)

    assert distrib == {"global" : [[1.]]}
    assert (grid.dx, grid.dy, grid.dz[0],
                grid.njx, grid.njy, grid.njz, 
                grid.isolated_box) == (test_input["voxel size"][0], test_input["voxel size"][1], test_input["voxel size"][2],
                                    expected[0], expected[1], expected[2],
                                    0)


def test_legumescene_to_RATPscene():
    la = numpy.zeros([1, 10, 4, 4])
    for i in range(4):
        for j in range(4):
            la[0][9][i][j] = 2.
    legumescene = { "LA" : la, "distrib" : numpy.array([[0.2, 0.2, 0.4, 0.1, 0.1]])}
    parameters = {"voxel size" : (1., 1., 1.)}

    grid, distrib, nb0 = legumescene_to_RATPscene(legumescene, parameters, coordinates, False, True)

    assert numpy.array_equal(distrib["global"], numpy.array([[0.2, 0.2, 0.4, 0.1, 0.1]]))
    assert nb0 == 9
    assert (grid.dx, grid.dy, grid.dz[0], 
                grid.njx, grid.njy, grid.njz, 
                grid.isolated_box,
                grid.s_canopy) == (parameters["voxel size"][0], parameters["voxel size"][1], parameters["voxel size"][2],
                                    4, 4, 1,
                                    0,
                                    32.)

def build_RATPscene_from_trimesh()