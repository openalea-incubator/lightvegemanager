'''
contient le main déroulé de la construction de la scène
'''

import numpy
import bisect

from src.LightVegeManager_trianglesmesh import *


def read_distrib_file(path, numberofentities) :
    f_angle = open(path, 'r')
    distrib = []
    for i in range(numberofentities):
        line = f_angle.readline()
        distrib.append([float(x) for x in line.split(',')[1:]])
    f_angle.close()

    return distrib

def compute_distrib_globale(trimesh, matching_ids, numberofclasses) :
    # dimensions [entité][classe d'angle]
    distrib = []

    angles = list(numpy.linspace(90/numberofclasses, 90, numberofclasses))
    
    # nombre d'entités
    numberofentities = max([v[1] for v in matching_ids.values()]) + 1
    for e in range(numberofentities) :
        tri_entity = triangles_entity(trimesh, e, matching_ids)
        area_entity = sum([triangle_area(t) for t in tri_entity])

        classes = [0] * numberofclasses
        if area_entity > 0 :
            for t in tri_entity :
                # id tq angles[id-1] < triangle_elevation(t) <= angles[id]
                id_classe = bisect.bisect_left(angles, triangle_elevation(t))
                classes[id_classe] += triangle_area(t)
            
            # ramène en relatif
            classes = [a/area_entity for a in classes]
        distrib.append(classes)

    return distrib

def compute_distrib_voxel(trimesh, 
                            matching_ids, 
                            numberofclasses, 
                            numberofvoxels, 
                            matching_tri_vox) :
    # on va jusqu'à 91° pour prendre en compte les plans
    angles = list(numpy.linspace(90/numberofclasses, 90, numberofclasses))
    
    # nombre d'entités
    numberofentities = max([v[1] for v in matching_ids.values()]) + 1

    # dimensions [#voxel][#entité][#classe angle]
    distrib = numpy.zeros([numberofvoxels,numberofentities,numberofclasses])

    # boucles suivent l'ordre des dimensions de distrib
    for k in range(numberofvoxels) :
        # liste de triangles dans le voxels
        tri_in_vox = [int(id) for id, v in matching_tri_vox.items() \
                                                            if int(v) == k]
        # on les range par entité        
        for n in range(numberofentities) :
            shape_ent = [ke for ke,v in matching_ids.items() if v[1] == n]
            # sélectionne les triangles appartenant à l'entité n
            tri_in_vox_in_ent = [i for i in tri_in_vox \
                        if globalid_to_elementid(trimesh, i) in shape_ent]
            
            # si l'entité est présente dans le voxel
            if tri_in_vox_in_ent :
                area_vox_ent = 0
                for id in tri_in_vox_in_ent :
                    t = globalid_to_triangle(trimesh, id)
                    id_classe = bisect.bisect_left(angles, \
                                                    triangle_elevation(t))
                    
                    at = triangle_area(t)
                    distrib[k][n][id_classe] += at
                    area_vox_ent += at
                
                distrib[k][n] = distrib[k][n] / area_vox_ent
    return distrib
            