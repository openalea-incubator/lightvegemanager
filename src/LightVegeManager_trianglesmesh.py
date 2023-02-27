'''
contient le création et gestion d'une triangulation

traduit les scenes en scene caribu, liste à extend à la scene finale

scene (dict): a {primitive_id: [triangles,]} dict.A triangle is a
                    list of 3-tuples points coordinates
        {
            0 : [
                    [(1,1,1), (0,0,0), (1,2,2)], 
                    [(10,10,10), (0,0,0), (10,20,20)]
                ]
        }

'''
import itertools
import numpy
import bisect

import openalea.plantgl.all as pgl
from alinea.caribu import plantgl_adaptor
from openalea.mtg.mtg import MTG

from src.LightVegeManager_basicgeometry import *

def triangles_entity(cscene, entity_id, matching_ids) :
    '''
    return [triangle,] triangles as [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]
    '''
    select_triangles = [t for k,t in cscene.items() \
                                        if matching_ids[k][1] == entity_id]
    return list(itertools.chain(*select_triangles))

def globalid_to_elementid(cscene, triangleid) :
    '''
    triangleid: from 0 to nb totale de triangles
    retourne l'indice de la shape auquelle appartient le triangle
    '''
    cumul_triangles_per_ele = numpy.cumsum([len(v) for v in cscene.values()])
    if triangleid >= cumul_triangles_per_ele[-1] :
        raise IndexError("id %i > number of triangles: %i" % (triangleid, \
                                                cumul_triangles_per_ele[-1]))

    # renvoie i tq cumul_triangles[i-1] <= triangleid < cumul_triangles[i]
    return bisect.bisect(cumul_triangles_per_ele.tolist(), triangleid)

def globalid_to_triangle(cscene, triangleid) :
    cumul_triangles_per_ele = numpy.cumsum([len(v) for v in cscene.values()])
    if triangleid >= cumul_triangles_per_ele[-1] :
        raise IndexError("id %i > number of triangles: %i" % (triangleid, \
                                                cumul_triangles_per_ele[-1]))
    ele = globalid_to_elementid(cscene, triangleid)
    return cscene[ele][triangleid - sum(cumul_triangles_per_ele[0:ele])]

def compute_area_max(cscene) :
    amax = -999999
    for t in itertools.chain(*[v for v in cscene.values()]) :
        a = triangle_area(t)
        if a > amax : amax = a
    return amax

def compute_minmax_coord(cscene) :
    (xmax, ymax, zmax) = (-999999 for i in range(3))
    (xmin, ymin, zmin) = (999999 for i in range(3))
    for t in itertools.chain(*[v for v in cscene.values()]) :
        for p in t :
            if p[0] > xmax :
                xmax = p[0]
            if p[0] < xmin:
                xmin = p[0]
            if p[1] > ymax :
                ymax = p[1]
            if p[1] < ymin:
                ymin = p[1]
            if p[2] > zmax :
                zmax = p[2]
            if p[2] < zmin:
                zmin = p[2]
    
    return [xmin, ymin, zmin], [xmax, ymax, zmax]

def compute_trilenght_max(cscene) :
    lenmax = -999999
    for t in itertools.chain(*[v for v in cscene.values()]) :
        (xmax, ymax, zmax) = (-999999 for i in range(3))
        (xmin, ymin, zmin) = (999999 for i in range(3))
        for p in t :
            if p[0] > xmax :
                xmax = p[0]
            if p[0] < xmin:
                xmin = p[0]
            if p[1] > ymax :
                ymax = p[1]
            if p[1] < ymin:
                ymin = p[1]
            if p[2] > zmax :
                zmax = p[2]
            if p[2] < zmin:
                zmin = p[2]
        # recherche de la différence de la longueur max d'un triangle 
        # max(xmax-xmin, ymax-ymin, zmax-zmin)
        m = max(xmax - xmin, ymax - ymin, zmax - zmin)
        if m > lenmax : lenmax = m
    
    return lenmax

def chain_triangulations(scenes) :
    complete_trimesh = {}
    matching_ids = {}
    legume_grid = False
    id_legume_scene = None
    
    element_count = 0
    for entity, scene in enumerate(scenes) :
        cscene = {}

        # scene planteGL
        if isinstance(scene, pgl.Scene):
            # les clés sont shape.id
            cscene = plantgl_adaptor.scene_to_cscene(scene)
        
        # fichier VGX
        elif isinstance(scene, str) and scene.split(".")[-1] == "vgx":
            # une seule clé = id_element courant
            cscene = vgx_to_caribu(scene, element_count)

        # scene déjà au format CARIBU
        elif isinstance(scene, list) and isinstance(scene[0], tuple) : 
            cscene = scene

        # table MTG avec un identifiant "geometry"
        elif isinstance(scene, MTG) :
            # lit la scene plantgl stockée dans le MTG key = shape.id
            cscene = plantgl_adaptor.mtg_to_cscene(scene)

        # grille de voxels au format l-egume (RiRi)
        elif isinstance(scene, dict) and len(scene) == 2 : 
            legume_grid = True
            id_legume_scene = entity

        # lisse les element_id et les ajoute à la triangulation globale
        for key, val in cscene.items() :
            matching_ids[element_count] = [key, entity]
            val = [list(a) for a in val] # conversion triangle en liste
            complete_trimesh[element_count] = val
            element_count += 1

    return complete_trimesh, matching_ids, legume_grid, id_legume_scene

def vgx_to_caribu(file, id_element) :
    f = open(file, 'r')
    lines = f.readlines()
    scene = {}
    triangles = []
    for l in lines[1:]:                    
        l_list = l.split("\t")

        # on considère comme feuille les éléments R (RGB) != 42
        if l_list[10] != '42' :
            tr = [(float(l_list[13]), float(l_list[14]), float(l_list[15])), \
                    (float(l_list[16]), float(l_list[17]), float(l_list[18])),\
                    (float(l_list[19]), float(l_list[20]), float(l_list[21]))]
            
            if triangle_area(tr) > 0:
                triangles.append(tr)
        
    f.close()
    
    scene[id_element] = triangles
    return scene

def apply_transformations(cscene, matching_ids, transformations, cscene_unit) :
    '''
    dans le dico transformation
    ["scenes unit] = {id_scene: unit}
    ["rescale"] = {id_scene: h}
    ["translate"] = {id_scene: vector}
    ["xyz orientation"] = {id_scene: orientation}
    
    return modification dans cscene
    
    '''

    # mise à l'échelle entre les unités de mesure
    if "scenes unit" in transformations :
        units = {'mm': 0.001, 
                'cm': 0.01, 
                'dm': 0.1, 
                'm': 1, 
                'dam': 10, 
                'hm': 100,
                'km': 1000}

        if cscene_unit not in units :
            raise ValueError("Unknown final scene unit: select one in this \
                                    list ['mm','cm','dm','m','dam','hm','km']")
        
        for key, value in transformations["scenes unit"].items() :
            scene_unit = transformations["scenes unit"]
            if scene_unit not in units :
                raise ValueError("Unknown scene %i unit: select one in this \
                            list ['mm','cm','dm','m','dam','hm','km']" % (key))

            if (scene_unit != cscene_unit) :
                h = units[scene_unit]/units[cscene_unit]
                for id in [k for k,v in matching_ids.items() if v[1] == key ] :
                    cscene[id] = rescale(cscene[id], h)

    if "rescale" in transformations :
        for key, value in transformations["rescale"].items() :
            for id in [k for k,v in matching_ids.items() if v[1] == key ] :
                cscene[id] = rescale(cscene[id], value)

    if "translate" in transformations :
        for key, value in transformations["translate"].items() :
            for id in [k for k,v in matching_ids.items() if v[1] == key ] :
                cscene[id] = translate(cscene[id], value)

    # on ramène la scène à la convention x+ = N
    if "xyz orientation" in transformations :
        for key, value in transformations["xyz orientation"].items() :
            if value == "x+ = S":
                rot = 180
            elif value == "x+ = W":
                rot = 90
            elif value == "x+ = E":
                rot = -90
            elif value == "x+ = N":
                rot = 0
            else :
                raise ValueError("Unknown xyz orientation in scene %i" % (key))
            for id in [k for k,v in matching_ids.items() if v[1] == key ] :
                cscene[id] = zrotate(cscene[id], rot)