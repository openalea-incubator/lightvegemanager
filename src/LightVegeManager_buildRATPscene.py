'''
Création d'une scène RATP avec une grille de voxels et une distribution
d'angle
'''

from PyRATP.pyratp import grid

from src.LightVegeManager_trianglesmesh import *
from src.LightVegeManager_voxelsmesh import *
from src.LightVegeManager_leafangles import *

def extract_grid_origin(parameters, minmax) :
    # origine de la grille
    if "origin" not in parameters : 
        [xorig, yorig, zorig] = minmax[0]
        zorig = -zorig
    else : 
        if len(parameters["origin"]) == 2 :
            xorig, yorig, zorig = parameters["origin"][0],  \
                                    parameters["origin"][1], \
                                    -minmax[0][2]
        elif len(parameters["origin"]) == 3 :
            [xorig, yorig, zorig] = parameters["origin"]
    
    return xorig, yorig, zorig

def build_RATPscene_from_trimesh(trimesh, 
                                    minmax, 
                                    triLmax,
                                    matching_ids, 
                                    parameters,
                                    coordinates,
                                    reflected,
                                    infinite,
                                    stems_id=None,
                                    nb_input_scenes=0) :
    # nombre d'entités
    numberofentities = max([v[1] for v in matching_ids.values()]) + 1
    
    ## construction de la distribution d'angle ##
    distrib = {}
    distrib_algo = parameters["angle distrib algo"]

    # calcul par entité
    if distrib_algo == "compute global" or \
        distrib_algo == "compute voxel" :
        distrib["global"] = compute_distrib_globale(    \
                                            trimesh,     \
                                            matching_ids, \
                                            parameters["nb angle classes"])
    # lecture d'un fichier
    elif distrib_algo == "file" :
        read_distrib_file(parameters["angle distrib file"], numberofentities)
    
    ## initialisation de la grille ##
    # taille des voxels
    if parameters["voxel size"] == "dynamic" :
        dv = [5*triLmax for i in range(3)]
    else:
        dv = parameters["voxel size"]

    # origine de la grille
    xorig, yorig, zorig = extract_grid_origin(parameters, minmax) 
    
    # nombre de voxels
    if "number voxels" in parameters :
        [nx, ny, nz] = parameters["number voxels"]

    else :
        grid_slicing = None
        if "grid slicing" in parameters : 
            grid_slicing = parameters["grid slicing"]
        
        nx, ny, nz = compute_grid_size_from_trimesh(minmax[0],
                                                    minmax[1],
                                                    dv, 
                                                    grid_slicing)

    # init grille
    # si pas de rayonnement réfléchi on annules la réflexion du sol
    if reflected : soil_reflectance = parameters["soil reflectance"]
    else : soil_reflectance = [0., 0.]
    ratpgrid = grid.Grid.initialise(nx, ny, nz, 
                                    dv[0], dv[1], dv[2],
                                    xorig, yorig, zorig, 
                                    coordinates[0], 
                                    coordinates[1], 
                                    coordinates[2], 
                                    numberofentities, 
                                    soil_reflectance, 
                                    toric=infinite)
        
    ## remplissage de la grille ##
    # tesselation des triangles dans les dimensions de la grille
    if parameters["tesselation level"] > 0 :
        tesselate_trimesh_on_grid(trimesh, ratpgrid, \
                                            parameters["tesselation level"])

    # remplissage
    ratpgrid, matching_tri_vox = fill_ratpgrid_from_trimesh(
                                                    trimesh, 
                                                    matching_ids, 
                                                    ratpgrid, 
                                                    stems_id,
                                                    nb_input_scenes)

    ## distribution des angles par voxel
    if distrib_algo == "compute voxel" :
        distrib["voxel"] = compute_distrib_voxel(trimesh, 
                                                matching_ids, 
                                                parameters["nb angle classes"], 
                                                int(ratpgrid.nveg), 
                                                matching_tri_vox)
    
    return ratpgrid, matching_tri_vox, distrib

def build_RATPscene_empty(parameters, minmax, coordinates, infinite) :
    # nombre de voxels
    if "number voxels" in parameters :
        [nx, ny, nz] = parameters["number voxels"]
    else :
        nx,ny,nz = 1,1,1
    
    # taille des voxels
    dv = parameters["voxel size"]
    
    # origine de la grille
    xorig, yorig, zorig = extract_grid_origin(parameters, minmax) 
    
    # init de la grille
    ratpgrid = grid.Grid.initialise(nx, ny, nz,
                                    dv[0], dv[1], dv[2],
                                    xorig, yorig, zorig,
                                    coordinates[0],
                                    coordinates[1],
                                    coordinates[2],
                                    1, 
                                    [0., 0.],
                                    toric=infinite)
    distrib = {"global" : [[1.]]}

    return ratpgrid, distrib
    
def legumescene_to_RATPscene(legumescene,
                                parameters,
                                coordinates,
                                reflected,
                                infinite) :
    '''
    
    legumescene est unique et regroupe toutes les entités
    legumescene = {
        "LA" : np.ndarray avec [nent,nz,ny,nx]
        "distrib" : np.ndarray [nent, nclasses]
    }

    '''
    
    # distribution d'angle
    distrib = {"global" : legumescene["distrib"]}

    # initialisation de la grille
    numberofentities = legumescene["LA"].shape[0]
    dv = parameters["voxel size"]
    dvolume = numpy.prod(dv)
    xorig, yorig, zorig = 0,0,0
    nx = legumescene["LA"].shape[3]
    ny = legumescene["LA"].shape[2]
    nz = legumescene["LA"].shape[1]

    # enlève les couches vides de l-egume au dessus du couvert
    nb0 = 0
    laicum = numpy.sum(legumescene["LA"], axis=0)
    laicumvert = numpy.sum(laicum, axis=(1, 2))
    for i in range(len(laicumvert)):
        if laicumvert[i] == 0.:
            nb0 += 1
        else:
            break
    nz = nz - nb0
    
    # initialisation de la grille
    # si pas de rayonnement réfléchi on annules la réflexion du sol
    if reflected : soil_reflectance = parameters["soil reflectance"]
    else : soil_reflectance = [0., 0.]
    ratpgrid = grid.Grid.initialise(nx, ny, nz, 
                                    dv[0], dv[1], dv[2],
                                    xorig, yorig, zorig, 
                                    coordinates[0], 
                                    coordinates[1], 
                                    coordinates[2], 
                                    numberofentities, 
                                    soil_reflectance, 
                                    toric=infinite)


    # remplissage de la grille
    fill_ratpgrid_from_legumescene(legumescene, ratpgrid, nb0, dvolume)

    return ratpgrid, distrib, nb0

