'''
fonction math pour comparer des vecteurs etc... selon les besoins

'''

import numpy
import math

def crossproduct(v1, v2) :
    x = v1[1] * v2[2] - v2[1] * v1[2]
    y = v1[2] * v2[0] - v2[2] * v1[0]
    z = v1[0] * v2[1] - v2[0] * v1[1]
    return (x,y,z)

def middle(v1, v2) :
    p = tuple([a+b for a,b in zip(v1, v2)])
    return (p[0]/2, p[1]/2, p[2]/2)

def triangle_normal(triangle) :
    side1 = [x-y for x, y in zip(triangle[1], triangle[0])]
    side2 = [x-y for x, y in zip(triangle[2], triangle[1])]
    side3 = [x-y for x, y in zip(triangle[0], triangle[2])]

    v12 = crossproduct(side1, side2)
    v23 = crossproduct(side2, side3)
    v31 = crossproduct(side3, side1)
    n = [x+y+z for x,y,z in zip (v12, v23, v31)]
    norm = math.sqrt(sum([c**2 for c in n]))
    return (n[0]/norm, n[1]/norm, n[2]/norm)

def triangle_elevation(triangle) :
    n = triangle_normal(triangle)
    # recupère l'élévation de la normale
    e = math.acos(abs(n[2]))
    # passe en degré
    e *= 180/math.pi
    # on reste entre 0 et 90
    if e > 90 and e < 180 : e = 90-(e%90)
    else : e = e%90
    
    return e

def triangle_area(triangle) :
    '''aire d'un triangle
        triangle : [(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]
        copie de _surf dans caributriangleset.py
    '''
    a, b, c = tuple(map(numpy.array, triangle))
    x, y, z = numpy.cross(b - a, c - a).tolist()
    return numpy.sqrt(x ** 2 + y ** 2 + z ** 2) / 2.0

def triangle_barycenter(triangle) :
    return tuple([s/3 for s in [x+y+z for x,y,z in zip(*triangle)]])

def rescale(triangles, h) :
    '''
    triangles est une liste de [(x1,y1,z1), (x2,y2,z2), (x3,y3,z3)]
    '''
    return [[tuple([x*h for x in p]) for p in t] for t in triangles] 

def translate(triangles, tvec) :
    '''
    tvec = (t1, t2, t3)
    '''
    return [[tuple([x+y for x,y in zip(p,tvec)]) for p in t] \
                                                    for t in triangles]

def zrotate(triangles, omegadeg) :
    '''
    angle en deg
    '''
    omegarad = omegadeg * math.pi/180
    newtriangles = []
    for t in triangles :
        newt = []
        for p in t :
            x = math.cos(omegarad)*p[0] - math.sin(omegarad)*p[1]
            y = math.sin(omegarad)*p[0] + math.cos(omegarad)*p[1]
            newt.append(tuple([x,y,p[2]]))
        newtriangles.append(newt)
    return newtriangles

