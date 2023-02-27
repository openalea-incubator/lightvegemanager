'''
gestion des entrées des modèles de lumières

'''

import openalea.plantgl.all as pgl

from src.LightVegeManager_trianglesmesh import *
from src.LightVegeManager_voxelsmesh import *

def Prepare_CARIBU(trimesh, 
                    geometry, 
                    matching_ids, 
                    minmax, 
                    parameters,
                    infinite, 
                    idsensors) :
    # propriétés optiques
    stems = []
    if 'stems id' in geometry : stems = geometry["stems id"]
    opt = CARIBU_opticals(matching_ids, 
                            parameters, 
                            stems)

    # debug, sensors, domain
    sensors_caribu = None
    if "sensors" in parameters and parameters["sensors"][0] == "grid" :
        dxyz = parameters["sensors"][1]
        nxyz = parameters["sensors"][2]
        orig = parameters["sensors"][3]
        arg = (dxyz,
                nxyz,
                orig,
                minmax[1],
                trimesh,
                matching_ids,
                idsensors,
                infinite)
        sensors_caribu, \
        sensors_plantgl, \
        Pmax_capt = create_caribu_legume_sensors(*arg)

    # pattern scene infini par défaut
    if "domain" not in geometry :
        # cas où il y a uniquement l-egume dans les scènes d'entrée
        if "sensors" in parameters and parameters["sensors"][0] == "grid" :
            reduction_domain = 0
            d = ((orig[0] + reduction_domain, 
                    orig[1] + reduction_domain), 
                        (nxyz[0] * dxyz[0] - reduction_domain, 
                            nxyz[1] * dxyz[1] - reduction_domain))
        else:
            d = ((minmax[0][0], minmax[0][1]), 
                    (minmax[1][0], minmax[1][1]))

        geometry["domain"] = d

    # On active l'option verbose de CARIBU
    debug = False
    if "debug" in parameters and parameters["debug"] :  debug = True
    
    return opt, sensors_caribu, debug
    

def CARIBU_opticals(matching_ids, parameters, stems_id) :
    opt = {}
    for band, coef in parameters["caribu opt"].items() :
        opt[band] = {}
    for id, val in matching_ids.items():
        for band, coef in parameters["caribu opt"].items() :
            # id : id de la shape, val : [id shape en input, id de l'entité]
            # c'est une tige il n'y a pas de transmission
            if (val[0], val[1]) in stems_id:
                #: (reflectance,) of the stems ou organe opaque
                opt[band][id] = (coef[0],)  
            else:
                # A 2-tuple encode a symmetric translucent material defined
                # by a reflectance and a transmittance
                # A 4-tuple encode an asymmetric translucent material defined
                # the reflectance and transmittance
                # of the upper and lower side respectively 
                #: (reflectance, transmittance) of the adaxial side of the 
                # leaves, élément translucide symétrique
                opt[band][id] = coef

    return opt

def create_caribu_legume_sensors(dxyz, 
                                nxyz, 
                                orig, 
                                pmax, 
                                trimesh, 
                                matching_ids, 
                                id_sensors, 
                                infinite) :
    # réduction si le nombre de couches remplies < nombre de couches prévues
    skylayer = reduce_layers_from_trimesh(trimesh, 
                                            pmax,
                                            dxyz,
                                            nxyz,
                                            matching_ids, 
                                            id_sensors)
            
    # scene plantGL qui acceuillera les capteurs
    s_capt = pgl.Scene()
    
    if infinite : 
        # capeur bas oriente vers le haut
        points = [(0, 0, 0), 
                    (dxyz[0], 0, 0), 
                    (dxyz[0], dxyz[1], 0), 
                    (0, dxyz[1], 0)]
    else :
        # capeur bas oriente vers le haut
        points = [(0, 0, 0),  
                    (0, dxyz[1], 0), 
                    (dxyz[0], dxyz[1], 0), 
                    (dxyz[0], 0, 0)]  
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]

    # carré taille d'un voxel
    carre = pgl.QuadSet(points, indices, normals, indices)
    
    ID_capt = 0
    dico_translat = {}
    for ix in range(nxyz[0]):
        for iy in range(nxyz[1]):
            for iz in range(nxyz[2] - skylayer):
                # vecteur de translation
                tx = orig[0] + ix * dxyz[0]
                ty = orig[1] + iy * dxyz[1]
                tz = orig[2] + iz * dxyz[2]

                # retient la translation
                dico_translat[ID_capt] = [tx, ty, tz]

                # ajoute un voxel à la scene des capteurs
                Vox = pgl.Translated(geometry=carre, translation=(tx, ty, tz))
                s_capt.add(pgl.Shape(geometry=Vox, id=ID_capt))
                ID_capt += 1

    
    Dico_Sensors = {}
    liste_hmax_capt = []
    liste_x_capt = []
    liste_y_capt = []
    # mise en forme de triangulation CARIBU
    for x in s_capt:

        # pgl_to_caribu(x)

        # Preparation capteurs
        pt_lst = x.geometry.geometry.pointList
        idx_lst = x.geometry.geometry.indexList

        for i in range(0, len(idx_lst)):
            x11 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
            y11 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
            z11 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

            liste_hmax_capt.append(z11)
            liste_x_capt.append(x11)
            liste_y_capt.append(y11)

            x12 = pt_lst[idx_lst[i][1]][0] + dico_translat[x.id][0]
            y12 = pt_lst[idx_lst[i][1]][1] + dico_translat[x.id][1]
            z12 = pt_lst[idx_lst[i][1]][2] + dico_translat[x.id][2]

            liste_hmax_capt.append(z12)
            liste_x_capt.append(x12)
            liste_y_capt.append(y12)

            x13 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
            y13 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
            z13 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

            liste_hmax_capt.append(z13)
            liste_x_capt.append(x13)
            liste_y_capt.append(y13)

            tple1 = [(x11, y11, z11), (x12, y12, z12), (x13, y13, z13)]
            triangle = []
            triangle.append(tple1)

            x21 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
            y21 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
            z21 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

            liste_hmax_capt.append(z21)
            liste_x_capt.append(x21)
            liste_y_capt.append(y21)

            x22 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
            y22 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
            z22 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

            liste_hmax_capt.append(z22)
            liste_x_capt.append(x22)
            liste_y_capt.append(y22)

            x23 = pt_lst[idx_lst[i][3]][0] + dico_translat[x.id][0]
            y23 = pt_lst[idx_lst[i][3]][1] + dico_translat[x.id][1]
            z23 = pt_lst[idx_lst[i][3]][2] + dico_translat[x.id][2]

            liste_hmax_capt.append(z23)
            liste_x_capt.append(x23)
            liste_y_capt.append(y23)

            tple2 = [(x21, y21, z21), (x22, y22, z22), (x23, y23, z23)]
            triangle.append(tple2)

            Dico_Sensors[x.id] = triangle
    
    return Dico_Sensors, s_capt, [(min(liste_x_capt) + max(liste_x_capt)) / 2,
                                (min(liste_y_capt) + max(liste_y_capt)) / 2,
                                max(liste_hmax_capt)]

def run_caribu(c_scene, direct_active, infinite, sensors) :
    if "sensors" is None :
        raw, aggregated = c_scene.run(direct=direct_active, infinite=infinite)
    else :
        raw, aggregated = c_scene.run(direct=direct_active, infinite=infinite, 
                                                            sensors=sensors)

    return raw, aggregated

