'''
    trianglesmesh
    **************

    It builds and handles triangulation mesh.

    The main triangulation format in LightVegeManager is the CARIBU format:
    scene (dict): a {primitive_id: [triangle,]} dict. A triangle is a list of 3-tuples points coordinates.
    example:
    
    .. code-block:: python
    
        scene = {   # element 0
                    0 : [ 
                            [(1,1,1), (0,0,0), (1,2,2)],  # triangle 0 of element 0
                            [(10,10,10), (0,0,0), (10,20,20)] # triangle 1 of element 0
                        ] ,
                    
                    # element 1
                    1 : [
                            [(12,12,12), (10,20,30), (21,22,22)], 
                            [(10,10,10), (10,20,30), (20,20,20)]
                        ]
                }

    .. seealso:: For more details :ref:`Inputs description <inputs>`

'''
import itertools
import numpy
import bisect

import openalea.plantgl.all as pgl
from alinea.caribu import plantgl_adaptor
from openalea.mtg.mtg import MTG

from lightvegemanager.basicgeometry import *

def triangles_entity(cscene, entity_id, matching_ids) :
    """Return a list of triangles belonging to the specy ``entity_id``

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param entity_id: specy indice in matching_ids
    :type entity_id: int
    :param matching_ids: 
        dict that matches new element indices in cscene with specy indice and
        input element indice, 
        .. code:: matching_ids = { new_element_id : (input_element_id, specy_id)}

        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :return: list of triangles belonging to ``entity_id`` in ``cscene``, a triangle is ``[(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]``
    :rtype: list
    """
    select_triangles = [t for k,t in cscene.items() \
                                        if matching_ids[k][1] == entity_id]
    return list(itertools.chain(*select_triangles))

def globalid_to_elementid(cscene, triangleid) :
    """Return the element index in cscene from a global triangle indice, 
    as if triangles were in a list without a sorting by element.

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param triangleid: indice from 0 to total number of triangles in cscene
    :type triangleid: int
    :raises IndexError: if triangleid > total number of triangles in cscene
    :return: element indice where triangleid belongs in cscene
    :rtype: int
    """    
    cumul_triangles_per_ele = numpy.cumsum([len(v) for v in cscene.values()])
    if triangleid >= cumul_triangles_per_ele[-1] :
        raise IndexError("id %i > number of triangles: %i" % (triangleid, \
                                                cumul_triangles_per_ele[-1]))

    # returns i as if cumul_triangles_per_ele[i-1] <= triangleid < cumul_triangles_per_ele[i]
    return bisect.bisect(cumul_triangles_per_ele.tolist(), triangleid)

def globalid_to_triangle(cscene, triangleid) :
    """Return the corresponding triangle in cscene from a global triangle indice,
    as if triangles were in a list without a sorting by element.

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param triangleid: indice from 0 to total number of triangles in cscene
    :type triangleid: int
    :raises IndexError: if triangleid > total number of triangles in cscene
    :return: Corresponding triangle in cscene
    :rtype: list of tuple
    """    
    cumul_triangles_per_ele = numpy.cumsum([len(v) for v in cscene.values()])
    if triangleid >= cumul_triangles_per_ele[-1] :
        raise IndexError("id %i > number of triangles: %i" % (triangleid, \
                                                cumul_triangles_per_ele[-1]))
    ele = globalid_to_elementid(cscene, triangleid)
    return cscene[ele][triangleid - sum(cumul_triangles_per_ele[0:ele])]

def compute_area_max(cscene) :
    """Maximum triangle area in a triangulation

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :return: maximum triangle area in ``cscene``
    :rtype: float
    """    
    amax = -999999
    for t in itertools.chain(*[v for v in cscene.values()]) :
        a = triangle_area(t)
        if a > amax : amax = a
    return amax

def compute_minmax_coord(cscene) :
    """Maximum and minimum point in all triangulation mesh

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :return: minimum point, maximum point
    :rtype: list, list
    """    
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
    """Search for maximum side of an axis oriented voxel where all triangles could fit in

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :return: maximum side of an axis oriented voxel where all triangles in cscene could fit in
    :rtype: float
    """    
    lenmax = -999999
    # iterate on all triangles
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
        m = max(xmax - xmin, ymax - ymin, zmax - zmin)
        if m > lenmax : lenmax = m
    
    return lenmax

def chain_triangulations(scenes) :
    """Aggregates all input geometric scenes in one global triangulation
        
        Current known formats:

            * plantGL scene
            * VGX file
            * CARIBU triangulation format (which is used in LightVegeManager for the global scene)
            * MTG table with ``"geometry"`` identifier
            * l-egume grid: dict of two entries, leaf area through a voxels grid and leaf angle distribution, for each specy

    Only l-egume grid can stores multiple species, otherwise each entry scene must represent one specy
    :param scenes: ``geometric["scenes"]`` from LightVegeManager inputs, list of geometric scenes. Scenes can be in different format in the list.
    :type scenes: list
    :return: 
        it returns 4 objects

            * ``complete_trimesh``: global triangulation which aggregates all input scenes. It is in CARIBU format
            * ``matching_ids``: dict which stores correspondances between new element indice, element indice in the input and specy indice
            * ``legume_grid``: boolean specifing if at least one l-egume grid is among the input scenes
            * ``id_legume_scene``: indice in the input scenes of a l-egume grid

    :rtype: dict of list, dict of list, bool, int
    """    
    complete_trimesh = {}
    matching_ids = {}
    legume_grid = False
    id_legume_scene = None
    
    element_count = 0
    for entity, scene in enumerate(scenes) :
        cscene = {}

        # scene planteGL
        if isinstance(scene, pgl.Scene):
            # keys/element indice are shape.id
            cscene = plantgl_adaptor.scene_to_cscene(scene)
        
        # fichier VGX
        elif isinstance(scene, str) and scene.split(".")[-1] == "vgx":
            # only one element id equal to element_count
            cscene = vgx_to_caribu(scene, element_count)

        # scene already in format CARIBU
        elif isinstance(scene, list) and isinstance(scene[0], tuple) : 
            cscene = scene

        # MTG table with "geometry" identifier
        elif isinstance(scene, MTG) :
            # mtg["geometry"] is plantGL scene
            cscene = plantgl_adaptor.mtg_to_cscene(scene)

        # voxels grid in l-egume format
        elif isinstance(scene, dict) and len(scene) == 2 : 
            legume_grid = True
            id_legume_scene = entity

        # creates a new numerotation for elements id adds them to the final triangulation
        for key, val in cscene.items() :
            matching_ids[element_count] = [key, entity]
            val = [list(a) for a in val] # converts triangles in list
            complete_trimesh[element_count] = val
            element_count += 1

    return complete_trimesh, matching_ids, legume_grid, id_legume_scene

def vgx_to_caribu(file, id_element) :
    """Reads VGX file and converts them in CARIBU scene format

    :param file: path name of the file
    :type file: string
    :param id_element: element indice where the triangulation will be stored
    :type id_element: int
    :return: triangulation in CARIBU format ``{id_element : [triangle,]}``
    :rtype: dict of list
    """    
    f = open(file, 'r')
    lines = f.readlines()
    scene = {}
    triangles = []
    for l in lines[1:]:                    
        l_list = l.split("\t")

        # we consider as leaf elements line with R (RGB) != 42
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
    """Applies geometric transformations on some part of the triangulation set in input datas
    Each transformation applies on all triangles from the same input scene.

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param matching_ids: 
        dict that matches new element indices in cscene with specy indice and
        input element indice, 
        
        .. code-block:: 
        
            matching_ids = { new_element_id : (input_element_id, specy_id)}

        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :param transformations: dict containing which geometric to apply on which element
    Possible transformations:

        * ``["scenes unit] = {id_scene: unit}``
        * ``["rescale"] = {id_scene: h}``
        * ``["translate"] = {id_scene: vector}``
        * ``["xyz orientation"] = {id_scene: orientation}``

    ``id_scene`` is the index in the input geometric scenes list

    :type transformations: dict of dict
    :param cscene_unit: measure unit of cscene
    :type cscene_unit: string
    :raises ValueError: unit is not in the list ``['mm','cm','dm','m','dam','hm','km']``
    :raises ValueError: unit is not in the list ``['mm','cm','dm','m','dam','hm','km']``
    :raises ValueError: xyz orientation is not one of ``"x+ = S", "x+ = E", "x+ = W", "x+ = N"``
    :return: ``cscene`` with transformations applied to inputs scene parts in the global mesh
    :rtype: dict of list
    """    
    # Rescale based on measure unit
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
        
        for key, scene_unit in transformations["scenes unit"].items() :
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