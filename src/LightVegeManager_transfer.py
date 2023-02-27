'''
gère le transfert vers les modèles de plantes

'''

import numpy
import scipy

def transfer_ratp_legume(m_lais, 
                            energy, 
                            ratp_grid,
                            voxels_outputs,
                            nb0,
                            epsilon) :
    # transfert des sorties
    res_abs_i = numpy.zeros((m_lais.shape[0], 
                                m_lais.shape[1], 
                                m_lais.shape[2], 
                                m_lais.shape[3]))
    # si le voxel est vide, on considère le transmis comme ce qui sort 
    # d'une de ses surfaces
    dS = ratp_grid.dx * ratp_grid.dy
    res_trans = numpy.ones((m_lais.shape[1], 
                            m_lais.shape[2], 
                            m_lais.shape[3]))
    res_trans = res_trans * (energy * dS)
            
    for ix in range(m_lais.shape[3]):
        for iy in range(m_lais.shape[2]):
            for iz in range(ratp_grid.njz):
                legume_iz = iz + nb0
                
                condition_x = (voxels_outputs.Nx==m_lais.shape[2] - iy)
                vox_data = voxels_outputs[ condition_x & 
                                            (voxels_outputs.Ny==ix+1) & 
                                            (voxels_outputs.Nz==iz+1)]
                
                a = min(sum(vox_data["Transmitted"]), dS)
                res_trans[legume_iz, iy, ix] = energy * a
                                        

                s_entity = 0
                for k in range(m_lais.shape[0]) : 
                    s_entity += m_lais[k][legume_iz][iy][ix]
                
                if s_entity > 0. :
                    for ie in range(m_lais.shape[0]) :
                        if len(vox_data) > 0 :
                            v_dat = vox_data[vox_data.VegetationType == ie+1]
                            v = v_dat["Intercepted"].values[0]
                            if v > epsilon :
                                res_abs_i[ie, legume_iz, iy, ix] = energy * v
                                                
                            # si le voxel est non vide, on fixe quand même une 
                            # valeur min
                            else:
                                res_abs_i[ie, legume_iz, iy, ix]

def transfer_caribu_legume(energy,
                            skylayer,
                            id, 
                            elements_outputs, 
                            sensors_outputs, 
                            sensors_dxyz, 
                            sensors_nxyz, 
                            m_lais, 
                            list_invar, 
                            list_lstring,
                            list_dicFeuilBilanR, 
                            infinite,
                            epsilon) :
    
    # réécriture de calc_paraF dans ShootMorpho.py pour correspondre aux 
    # triangles
    # parip : somme des organes sur tous les statuts
    # parap : idem sauf si plante senescente alors == 0

    ## Calcul du rayonnement absorbé par plante et par espèce
    for k in range(len(list_invar)) :
        # on initialise la somme sur chaque plante à 0
        nplantes = len(list_invar[k]['Hplante'])
        list_invar[k]['parap'] = scipy.array([0.] * nplantes)
        list_invar[k]['parip'] = scipy.array([0.] * nplantes)

        if id == None:
            filter = elements_outputs.VegetationType == k
            ent_organs_outputs = elements_outputs[filter]
        elif type(id) == list or type(id) == tuple :
            filter = elements_outputs.VegetationType == id[k]
            ent_organs_outputs = elements_outputs[filter]

        # # scene non vide
        for i in range(len(ent_organs_outputs)) :
            organe_id = int(ent_organs_outputs.iloc[i]["Organ"])

            # PAR en W/m²
            par_intercept = ent_organs_outputs.iloc[i]['par Ei'] * energy
            S_leaf = ent_organs_outputs.iloc[i]['Area']
            
            id_plante = list_lstring[k][organe_id][0]
            p_s = par_intercept * S_leaf
            a = float(list_invar[k]['parip'][id_plante])
            list_invar[k]['parip'][id_plante] = a + p_s
                
            # on enlève les feuilles senescentes 
            if list_lstring[k][organe_id][9] != 'sen' :
                a =  float(list_invar[k]['parap'][id_plante])
                list_invar[k]['parap'][id_plante] = a + p_s

        
        # on pose une valeur > 0 pour les plantes avec des feuilles
        for p in range(len(list_invar[k]['parip'])) :
            if list_invar[k]['parip'][p] == 0. and \
                                    list_dicFeuilBilanR[k]["surf"][p] > 0. :
                list_invar[k]['parip'][p] = epsilon
            
        # conversion
        c = (3600 * 24)/1000000
        list_invar[k]['parap'] *= c
        list_invar[k]['parip'] *= c
    
    ## calcul de res_trans : rayonnement transmis à travers la grille de voxel
    res_trans = numpy.ones((m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))

    # si scène non vide
    if elements_outputs:
        # traitements des sensors différents si scène infinie ou non
        if infinite :         
            ID_capt = 0
            for ix in range(sensors_nxyz[0]):
                for iy in range(sensors_nxyz[1]):
                    for iz in range(sensors_nxyz[2] - skylayer):
                        a =  min(sensors_outputs['par'][ID_capt], 1)
                        res_trans[(sensors_nxyz[2]-1) - iz][iy][ix] = a
                        ID_capt += 1
        
        else :
            raise ValueError("CARIBU Sensor + no infinite -> Doesn't work yet")

    # surface d'une face d'un voxel
    dS = sensors_dxyz[0] * sensors_dxyz[1]
    res_trans = res_trans * energy * dS

    return res_trans

    