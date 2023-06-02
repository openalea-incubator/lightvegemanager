from lightvegemanager.basicgeometry import *

def test_crossproduct():
    v1 = (2, 3, 4)
    v2 = (3, 0, 1)
    assert crossproduct(v1, v2) == (3, 10, -9)

def test_middle():
    v1 = (2, 3, 4)
    v2 = (3, 0, 1)
    assert middle(v1, v2) == (2.5, 1.5, 2.5)

def test_triangle_normal():
    triangle = [(0, 0, 0),(0, 1, 0),(1, 1, 1)]
    assert triangle_normal(triangle) == (math.sqrt(2)/2, 0.0, -math.sqrt(2)/2)

def test_triangle_elevation() :
    triangle = [(0, 0, 0),(0, 1, 0),(1, 1, 1)]
    assert triangle_elevation(triangle) == 45.

def test_triangle_area() :
    triangle = [(0, 0, 0),(0, 1, 0),(1, 1, 1)]
    assert triangle_area(triangle) == math.sqrt(2)/2

def test_triangle_barycenter():
    triangle = [(0, 0, 0),(0, 1, 0),(1, 1, 1)]
    assert triangle_barycenter(triangle) == (1/3, 2/3, 1/3)

def test_rescale():
    triangles = [[(0,0,0), (0,1,0), (0,1,1)], \
                    [(0,0,1), (0,0,0), (1,0,1)]]
    assert rescale(triangles, 3) == [[(0,0,0), (0,3,0), (0,3,3)], \
                                        [(0,0,3), (0,0,0), (3,0,3)]]

def test_translate():
    triangles = [[(0,0,0), (0,1,0), (0,1,1)], 
                    [(0,0,1), (0,0,0), (1,0,1)]]
    tvec = (3, 0, 1)
    assert translate(triangles, tvec) == [[(3,0,1), (3,1,1), (3,1,2)], \
                                            [(3,0,2), (3,0,1), (4,0,2)]]
    
def test_zrotate():
    triangles = [[(0,0,0), (0,1,0), (0,1,1)], \
                    [(0,0,1), (0,0,0), (1,0,1)]]
    r3_2 = math.sqrt(3)/2
    assert zrotate(triangles, 60) == [[(0.0, 0.0, 0), \
                                    (-r3_2, 0.5000000000000001, 0),\
                                    (-r3_2, 0.5000000000000001, 1)],\
                                    [(0.0, 0.0, 1), (0.0, 0.0, 0), (0.5000000000000001, r3_2, 1)]]