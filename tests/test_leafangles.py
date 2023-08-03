from lightvegemanager.basicgeometry import triangle_area
from lightvegemanager.leafangles import read_distrib_file, compute_distrib_globale, compute_distrib_voxel
import os
import numpy


def test_read_distrib_file():
    expected = [0.1382, 0.1664, 0.1972, 0.1925, 0.1507, 0.0903, 0.0425, 0.0172, 0.005]

    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "luzerne_angle_distrib.data")
    result = read_distrib_file(path, 1)

    numpy.testing.assert_array_equal(expected, result[0])


# global variables for the rest of the tests
trimesh = {
    0: [
        [[0, 0, 0], [0, 1, 1], [1, 1, 1]],  # e = 45°, a = 0.7071
        [[1, 0, 0], [1, 2, 2], [2, 2, 2]],  # e = 45°, a = 1.4142
        [[0, -1, 0], [0, 1, 1], [2, 1, 1]],  # e = 26.56°, a = 2.2361
        [[0, -1, 10], [0, 1, 11], [2, 1, 11]],  # e = 26.56°, a = 2.2361
        [[-2, -1, 10], [-2, 1, 11], [0, 1, 11]],  # e = 26.56°, a = 2.2361
    ]
}
total_area = sum([triangle_area(t) for t in trimesh[0]])

# identification dict
matching_ids = {0: [0, 0]}
numberofclasses = 3


def test_compute_distrib_globale():
    result = compute_distrib_globale(trimesh, matching_ids, numberofclasses)

    # expected
    expected = [
        [
            3 * triangle_area(trimesh[0][2]) / total_area,
            (triangle_area(trimesh[0][0]) + triangle_area(trimesh[0][1])) / total_area,
            0.0,
        ]
    ]

    numpy.testing.assert_array_equal(expected[0], result[0])


def test_compute_distrib_voxel():
    numberofvoxels = 2
    matching_triangle_voxel = {0: 0, 1: 1, 2: 0, 3: 1, 4: 1}

    result = compute_distrib_voxel(trimesh, matching_ids, numberofclasses, numberofvoxels, matching_triangle_voxel)

    area_v1 = triangle_area(trimesh[0][0]) + triangle_area(trimesh[0][2])
    area_v2 = triangle_area(trimesh[0][1]) + 2 * triangle_area(trimesh[0][2])
    expected = [
        [[triangle_area(trimesh[0][2]) / area_v1, triangle_area(trimesh[0][0]) / area_v1, 0.0]],
        [[2 * triangle_area(trimesh[0][2]) / area_v2, triangle_area(trimesh[0][1]) / area_v2, 0.0]],
    ]

    for i in range(2):
        numpy.testing.assert_array_equal(expected[i][0], result[i][0])
