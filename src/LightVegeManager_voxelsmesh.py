'''
contient le création et gestion d'une grille de voxel
'''

import itertools

from PyRATP.pyratp import grid

from src.LightVegeManager_basicgeometry import *
from src.LightVegeManager_tesselator import *
from src.LightVegeManager_trianglesmesh import *

def compute_grid_size_from_trimesh(pmin, pmax, dv, grid_slicing=None) :
    [dx, dy, dz] = dv

    # on ajuste les min-max si la scène est plane pour avoir un espace 3D
    for i in range(3) :
        if pmin[i] == pmax[i]:
            pmax[i] += dv[i]
            pmin[i] -= dv[i]

    # calcul
    nx = int((pmax[0] - pmin[0]) // dx)
    ny = int((pmax[1] - pmin[1]) // dy)
    if grid_slicing is not None :
        if grid_slicing == "ground = 0." :
            nz = int((pmax[2] - 0.) // dz)
    else :
        nz = int((pmax[2] - pmin[2]) // dz)
    if (pmax[0] - pmin[0]) % dx > 0 : nx += 1
    if (pmax[1] - pmin[1]) % dy > 0 : ny += 1
    if (pmax[2] - pmin[2]) % dz > 0 : nz += 1

    return nx, ny, nz

def tesselate_trimesh_on_grid(trimesh, ratpgrid, levelmax) :
    # traite les triangles de la shape
    new_trimesh = {}
    for id_ele, triangles in trimesh.items() :
        new_tr_scene = []
        for t in triangles :
            level = 0
            iterate_trianglesingrid(t, ratpgrid, level, levelmax, new_tr_scene)
        new_trimesh[id_ele] = new_tr_scene

    return new_trimesh

def fill_ratpgrid_from_trimesh(trimesh, 
                                matching_ids, 
                                ratpgrid, 
                                stems_id,
                                nb_input_scenes=0) :
    # pour chaque triangle, indice entité, x, y, z, aire, azote
    entity, barx, bary, barz, a, n = [],[], [], [], [], []
    for id, ts in trimesh.items() :
        for t in ts :
            bar = triangle_barycenter(t)
            barx.append(bar[0])
            bary.append(bar[1])
            barz.append(bar[2])
            
            # id : id de la shape, val : [id shape en input, id de l'entité]
            # c'est une tige on divise par 2 le LAD
            if isinstance(stems_id, list) :
                if (matching_ids[id][0], 
                    matching_ids[id][1] - nb_input_scenes) in stems_id or \
                    (matching_ids[id][0], matching_ids[id][1]) in stems_id :
                    a.append(triangle_area(t) / 2)
            else:
                a.append(triangle_area(t))
            
            n.append(0.)
            entity.append(matching_ids[id][1])

    return grid.Grid.fill_1(entity, barx, bary, barz, a, n, ratpgrid)

def fill_ratpgrid_from_legumescene(legumescene, ratpgrid, nb0, dvolume) :
    # remplissage des voxels avec la grille 
    n_vox_per_nent = []
    for ne in range(ratpgrid.nent):
        k=0
        for ix in range(ratpgrid.njx):
            for iy in range(ratpgrid.njy):
                for iz in range(ratpgrid.njz):
                    legume_iz = iz + nb0

                    # on force les voxels vides dans les couches non vides à 
                    # être interprétés par RATP
                    S_voxel = max(1e-14, \
                                    legumescene["LA"][ne][legume_iz][iy][ix])

                    #ajouter 1 pour utilisation f90
                    ratpgrid.kxyz[ratpgrid.njy-(iy+1), ix, iz] = k + 1 
                    ratpgrid.numx[k] = (ratpgrid.njy-(iy+1)) + 1
                    ratpgrid.numy[k] = ix + 1
                    ratpgrid.numz[k] = iz + 1
                    ratpgrid.nume[ne,k] = ne + 1
                    ratpgrid.nje[k] = max(ne + 1, ratpgrid.nje[k])
                    ratpgrid.nemax = max(ratpgrid.nemax, ratpgrid.nje[k])

                    ratpgrid.leafareadensity[ne,k] += S_voxel / dvolume
                    ratpgrid.s_vt_vx[ne,k] += S_voxel
                    ratpgrid.s_vx[k] += S_voxel
                    ratpgrid.s_vt[ne] += S_voxel
                    ratpgrid.s_canopy += S_voxel
                        
                    k=k+1
        n_vox_per_nent.append(k)

    ratpgrid.nveg=max(n_vox_per_nent)
    ratpgrid.nsol=ratpgrid.njx*ratpgrid.njy   # Numbering soil surface areas
    for jx in range(ratpgrid.njx):
        for jy in range(ratpgrid.njy):
            ratpgrid.kxyz[jx,jy,ratpgrid.njz] = ratpgrid.njy * jx + jy + 1

    for k in range(ratpgrid.nveg):
        for je in range(ratpgrid.nje[k]):
            if je == 0:
                # !!! on considère l-egume avec une hauteur fixe de voxel
                # Incrementing total canopy volume
                ratpgrid.volume_canopy[ratpgrid.nent] = \
                                ratpgrid.volume_canopy[ratpgrid.nent] + dvolume  
            if  ratpgrid.s_vt_vx[je,k] > 0. :
                ratpgrid.volume_canopy[ratpgrid.nume[je,k] - 1] = \
                    ratpgrid.volume_canopy[ratpgrid.nume[je,k] - 1] + dvolume
                ratpgrid.voxel_canopy[ratpgrid.nume[je,k] - 1] =  \
                    ratpgrid.voxel_canopy[ratpgrid.nume[je,k] - 1] + 1


def reduce_layers_from_trimesh(trimesh, 
                                pmax,
                                dxyz,
                                nxyz,
                                matching_ids, 
                                ids=None) :
    # réduction si le nombre de couches remplies < nombre de couches prévues
    if ids is None :
        skylayer = (pmax[2]) // dxyz[2]
        if skylayer < nxyz[2] and  pmax[2] > 0 : 
            skylayer = int(nxyz[2] - 1 - skylayer)
    
        # autrement on garde le nombre de voxels prévus
        else : skylayer = 0 
    
    # on calcule le z max de la scène en omettant certaines espèces
    else:
        # on regarde d'abord si l'indice des capteurs est dans la scène 
        # d'entrée (cas si scene plantgl vide)
        i = 0
        specie_empty = matching_ids[i][1] not in ids
        while ((matching_ids[i][1] not in ids) and \
                    (i+1 < len(matching_ids))) :
            i+=1
            if matching_ids[i][1] in ids : specie_empty = False
        
        if specie_empty : skylayer = nxyz[2] - 1
        else :
            # zmax de la scène
            zmax = -999999  
            if trimesh :
                for i in ids :
                    triangles = triangles_entity(trimesh, i, matching_ids)
                    for z in list(itertools.chain(*[[p[2] for p in t] \
                                                    for t in triangles])) :
                            if z > zmax :
                                zmax = z

            skylayer = zmax // dxyz[2]
            if skylayer < nxyz[2] and  zmax > 0 : 
                skylayer = int(nxyz[2] - 1 - skylayer)
    
    return skylayer