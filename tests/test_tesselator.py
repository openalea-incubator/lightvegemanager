from math import e
from lightvegemanager.tesselator import (
    tesselate,
    tesselate2,
    iterate_triangles,
    iterate_trianglesingrid,
)
from alinea.pyratp import grid

pars = {"latitude": 0, "longitude": 0, "timezone": 0, "nent": 1, "rs": (0, 0), "idecaly": 0, "orientation": 0}
pars.update({"njx": 2, "njy": 2, "njz": 2, "dx": 1.0, "dy": 1.0, "dz": [1.0] * 2, "xorig": 0, "yorig": 0, "zorig": 0})
mygrid = grid.Grid.initialise(**pars)
triangle_in = ((0.25, 0.25, 0.5), (0.25, 0.75, 0.5), (0.75, 0.75, 0.5))
triangle_between = ((0.25, 0.25, 0.5), (0.25, 1.75, 0.5), (1.75, 1.75, 0.5))


def test_tesselate():
    triangles_processed = tesselate(mygrid, triangle_in)
    assert len(triangles_processed) == 1

    triangles_processed = tesselate(mygrid, triangle_between)
    assert len(triangles_processed) == 4


def test_tesselate2():
    triangles_processed = tesselate2(triangle_in)
    assert len(triangles_processed) == 4

    triangles_processed = tesselate2(triangle_between)
    assert len(triangles_processed) == 4


def test_iterate_trianglesingrid():
    levelmax = 4

    triangles_shape = []
    iterate_trianglesingrid(triangle_in, mygrid, 0, levelmax, triangles_shape)
    assert len(triangles_shape) == 1

    triangles_shape = []
    iterate_trianglesingrid(triangle_between, mygrid, 0, levelmax, triangles_shape)
    assert len(triangles_shape) == 34


def test_iterate_triangles():
    levelmax = 4

    triangles_shape = []
    iterate_triangles(triangle_in, 0, levelmax, triangles_shape)
    assert len(triangles_shape) == 64

    triangles_shape = []
    iterate_triangles(triangle_between, 0, levelmax, triangles_shape)
    assert len(triangles_shape) == 64
