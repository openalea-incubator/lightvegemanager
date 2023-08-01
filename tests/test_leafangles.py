from lightvegemanager.leafangles import read_distrib_file, compute_distrib_globale, compute_distrib_voxel
import pytest
import os


def test_read_distrib_file():
    expected = [0.1382, 0.1664, 0.1972, 0.1925, 0.1507, 0.0903, 0.0425, 0.0172, 0.005]

    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "luzerne_angle_distrib.data")
    result = read_distrib_file(path, 1)

    for x, y in zip(expected, result[0]):
        assert x == y


test_read_distrib_file()

# def test_compute_distrib_globale():


# def test_compute_distrib_voxel():
