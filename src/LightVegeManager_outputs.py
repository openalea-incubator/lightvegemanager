'''
gère et reformate les sorties des modèles de lumières

'''

import pandas
import itertools

from src.LightVegeManager_trianglesmesh import *

def out_ratp_empty_grid(day, hour) :
    return  pandas.DataFrame({'VegetationType':[0],
                                'Iteration':[1],
                                'Day':day,
                                'Hour':hour,
                                'Voxel':[0],
                                'Nx':[0],
                                'Ny':[0],
                                'Nz':[0],
                                'ShadedPAR':[0],
                                'SunlitPAR':[0],
                                'ShadedArea':[0],
                                'SunlitArea': [0],
                                'Area': [0],
                                'PARa': [0],
                                'Intercepted': [0],
                                'Transmitted': [0]
                            })

def out_ratp_voxels(ratpgrid, res, parunit) :
    # récupère les sorties de RATP
    # np.array en une dimension, de taille nbvoxels x nbiteration
    VegetationType, \
    Iteration,      \
    day,            \
    hour,           \
    VoxelId,        \
    ShadedPAR,      \
    SunlitPAR,      \
    ShadedArea,     \
    SunlitArea,     \
    xintav,         \
    Ptransmitted = res.T

    if parunit == "W.m-2" :
        # ('PAR' is expected in  Watt.m-2 in RATP input, whereas output is in 
        # micromol => convert back to W.m2 (cf shortwavebalance, line 306))
        # on reste en micromol !
        # On enregistre tout dans une dataframe pandas
        ShadedPAR = ShadedPAR / 4.6
        SunlitPAR = SunlitPAR / 4.6
    
    para_list=[]
    for i in range(len(ShadedPAR)):
        if (ShadedArea[i] + SunlitArea[i]) > 0 :
            shaded = ShadedPAR[i] * ShadedArea[i]
            sunlit = SunlitPAR[i] * SunlitArea[i]
            area = ShadedArea[i] + SunlitArea[i]
            para_list.append((shaded + sunlit )/area)
        else:
            para_list.append(0.)
    
    # vérifie qu'on a pas de erel négatif
    erel_list=[]
    for i in range(len(xintav)):
        if xintav[i] >= 1e-6 :
            erel_list.append(xintav[i])
        else:
            erel_list.append(0.)

    # liste des indices
    numx=[]
    numy=[]
    numz=[]
    for v in VoxelId :
        numx.append(ratpgrid.numx[int(v)-1])
        numy.append(ratpgrid.numy[int(v)-1])
        numz.append(ratpgrid.numz[int(v)-1])

    dfvox =  pandas.DataFrame({'VegetationType':VegetationType,
                                    'Day':day,
                                    'Hour':hour,
                                    'Voxel':VoxelId,
                                    'Nx':numx,
                                    'Ny':numy,
                                    'Nz':numz,
                                    'ShadedPAR':ShadedPAR,
                                    'SunlitPAR':SunlitPAR,
                                    'ShadedArea':ShadedArea,
                                    'SunlitArea': SunlitArea,
                                    'Area': ShadedArea + SunlitArea,
                                    'PARa': para_list,
                                    'Intercepted': erel_list, 
                                    'Transmitted': Ptransmitted
                                })
                
    # ne prend pas le sol
    return dfvox[dfvox['VegetationType'] > 0]

def out_ratp_triangles(trimesh,
                        matching_ele_ent,
                        matching_tr_vox,
                        voxels_outputs) :
    # création du tableau résultats des triangles
    entity = {}
    for id, match in matching_ele_ent.items():
        entity[id] = match[1] + 1
    
    index = range(len(matching_tr_vox))
    vox_id = [matching_tr_vox[str(i)] + 1 for i in index]
    triangles = list(itertools.chain(*trimesh.values()))
    sh_id = [globalid_to_elementid(trimesh, i) for i in range(len(triangles))]
    s = [triangle_area(t) for t in triangles]

    # nouvelle data frame avec les triangles en index
    dftriangles = pandas.DataFrame({'Triangle': index,
                                'Organ': sh_id, 
                                'Voxel':vox_id, 
                                'VegetationType':[entity[id] for id in sh_id], 
                                'primitive_area':s})

    # supposé copie dans les index avec des colonnes en commun
    # colonnes en commun : VegetationType, VoxelId
    trianglesoutputs = pandas.merge(dftriangles, voxels_outputs)        
    
    # tri les lignes par ordre de triangles
    return  trianglesoutputs.sort_values('Triangle')

def out_ratp_elements(matching_ele_ent, 
                        reflected, 
                        reflectance_coef, 
                        trianglesoutputs) :
    # enregistre les valeurs par shape et plantes
    nshapes = len(matching_ele_ent)
    s_shapes = []
    s_area=[]
    s_para=[]
    s_pari=[]
    s_intercepted=[]
    s_transmis=[]
    s_day=[]
    s_hour=[]
    s_ent=[]
    s_parsun=[]
    s_parsha=[]
    s_areasun=[]
    s_areasha=[]
    for id in range(nshapes):
        # itérations commencent à 1
        nent = matching_ele_ent[id][1]
        dffil = trianglesoutputs[(trianglesoutputs.Organ == id)]
        
        t_areas = dffil["primitive_area"]
        sum_area = sum(t_areas)
        
        s_hour.append(dffil["Hour"].values[0])
        s_day.append(dffil["Day"].values[0])
        s_area.append(sum_area)

        # le rayonnemenbt réfléchi est calculé dans RATP
        if reflected :
            para = sum(t_areas * dffil['PARa']) / sum_area
            s_para.append(para)
        
        # sinon on le rajoute manuellement
        else:
            # PAR incident
            pari = sum(t_areas * dffil['PARa']) / sum_area
            s_pari.append(pari)
            s_para.append(pari - (pari * reflectance_coef[nent][0]))
        
        s_parsun.append(sum(t_areas * dffil['SunlitPAR']) / sum_area)
        s_parsha.append(sum(t_areas * dffil['ShadedPAR']) / sum_area)
        s_areasun.append(sum(t_areas * dffil['SunlitArea']) / sum_area)
        s_areasha.append(sum(t_areas * dffil['ShadedArea']) / sum_area)
        s_intercepted.append(sum(t_areas * dffil['Intercepted']) /sum_area)
        s_transmis.append(sum(t_areas * dffil['Transmitted']) /sum_area)
        s_ent.append(dffil["VegetationType"].values[0])
        s_shapes.append(matching_ele_ent[id][0])
    
    return pandas.DataFrame({
                                "Day" : s_day,
                                "Hour" : s_hour,
                                "Organ" : s_shapes,
                                "VegetationType" : s_ent,
                                "Area" : s_area,
                                "PARa" : s_para,
                                "Intercepted" : s_intercepted,
                                "Transmitted" : s_intercepted,
                                "SunlitPAR" : s_parsun,
                                "SunlitArea" : s_areasun,
                                "ShadedPAR" : s_parsha,
                                "ShadedArea" : s_areasha
                            })

def out_caribu_mix( rdrs,
                    c_scene_sun, c_scene_sky,
                    raw_sun, aggregated_sun, 
                    raw_sky, aggregated_sky,
                    issensors,
                    issoilmesh) :

    raw, aggregated, sensors, soil_energy = raw_sun, aggregated_sun, {}, {}
    for band in raw.keys() : 
        for ray in ['Eabs', 'Ei'] :
            for idt, v_sun in raw[band][ray].items() :
                v_sky = raw_sky[band][ray][idt]
                raw[band][ray][idt] = rdrs * v_sky + (1 - rdrs) * v_sun

            for ide, v_sun in aggregated[band][ray].items() :
                v_sky = aggregated_sky[band][ray][ide]
                aggregated[band][ray][ide] = rdrs * v_sky + (1 - rdrs) * v_sun

        if issensors :
            sensors[band] = {}
            for id, v_sun in aggregated_sun[band]['sensors']['Ei'].items() :
                v_sky = aggregated_sky[band]['sensors']['Ei'][id]
                sensors[band][id] = rdrs * v_sky + (1 - rdrs) * v_sun

    if issoilmesh : 
        q1, e1 = c_scene_sky.getSoilEnergy()
        q2, e2 = c_scene_sun.getSoilEnergy()
        q = rdrs * q1 + (1 - rdrs) * q2
        e = rdrs * e1 + (1 - rdrs) * e2

        soil_energy['Qi'] = q
        soil_energy['Einc'] = e

    return raw, aggregated, sensors, soil_energy

def out_caribu_nomix(c_scene, aggregated, issensors, issoilmesh) :
    sensors, soil_energy = {}, {}
    
    if issensors :
        for band in aggregated.keys() :
            sensors[band] =  aggregated[band]['sensors']['Ei']

    if issoilmesh :
        q, e = c_scene.getSoilEnergy()
        
        soil_energy['Qi'] = q
        soil_energy['Einc'] = e

    return sensors, soil_energy

def out_caribu_elements(day, hour, trimesh, matching_ids, aggregated, sun_up) :

    (s_shapes, 
    s_area, 
    s_day, 
    s_hour, 
    s_ent) = ([0] * len(matching_ids) for i in range(5))
    
    for key,val in matching_ids.items():
        s_shapes[key] = val[0]
        s_area[key] = sum([triangle_area(t) for t in trimesh[key]])
        s_day[key] = day
        s_hour[key] = hour
        s_ent[key] = matching_ids[key][1]

    dico_shape = {
        "Day" : s_day,
        "Hour" : s_hour,
        "Organ" : s_shapes,
        "VegetationType" : s_ent,
        "Area" : s_area
    }
    
    s_Eabs = {}
    s_Ei = {}
    for band, dict_val in aggregated.items() : 
        s_Eabs[band] = [0]*len(matching_ids)
        s_Ei[band] = [0]*len(matching_ids)
        for key,val in matching_ids.items():
            if sun_up: s_Eabs[band][key] = dict_val['Eabs'][key]
            if sun_up: s_Ei[band][key] = dict_val['Ei'][key]
    
        dico_shape[band + " Eabs"] = s_Eabs[band]
        dico_shape[band + " Ei"] = s_Ei[band]

    return pandas.DataFrame(dico_shape)

def out_caribu_triangles(day, hour, trimesh, matching_ids, raw, sun_up) :
    nb_triangles = sum([len(triangles) for triangles in trimesh.values()])
    
    (s_shapes, 
    s_tr,
    s_area, 
    s_day, 
    s_hour, 
    s_ent) = ([0] * nb_triangles for i in range(6))

    id_tr = 0
    for ele, triangles in trimesh.items():
        for t in triangles :
            s_shapes[id_tr] = ele
            s_tr[id_tr] = id_tr
            s_area[id_tr] = triangle_area(t)
            s_day[id_tr] = day
            s_hour[id_tr] = hour
            s_ent[id_tr] = matching_ids[ele][1]
            id_tr += 1

    dico_tr = {
        "Day" : s_day,
        "Hour" : s_hour,
        "Triangle" : s_tr,
        "Organ" : s_shapes,
        "VegetationType" : s_ent,
        "Area" : s_area
    }
    
    s_Eabs = {}
    s_Ei = {}
    for band, dict_val in raw.items() : 
        if sun_up: s_Eabs[band] = \
            list(itertools.chain(*[dict_val['Eabs'][ele] for ele in trimesh.keys()]))
        if sun_up: s_Ei[band] =  \
            list(itertools.chain(*[dict_val['Ei'][ele] for ele in trimesh.keys()]))
    
        dico_tr[band + " Eabs"] = s_Eabs[band]
        dico_tr[band + " Ei"] = s_Ei[band]

    return pandas.DataFrame(dico_tr)