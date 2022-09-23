from ftplib import error_reply
import sys
import os
import getopt
import time

import numpy as np
import pandas as pd

try :
    from src.LightVegeManager import *
    from src.l_egume_template import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
    from src.LightVegeManager import *
    from src.l_egume_template import *

'''
Comparaison de la lumière entre l-egume (RiRi) et RATP

'''

def simulation(foldin, foldout, active, passive, ratpgeo, writegeo=False,triangles="n"):
    
    start = time.time()
    fusms = "liste_usms_exemple.xls"
    ongletBatch = "exemple"

    # prépare une liste de l-system à lancer
    lsystem_simulations, names_simulations = initialisation_legume(foldin, foldout, fusms, ongletBatch)

    index_simulation = 0
    sim_id = names_simulations[index_simulation]
    lstring = lsystem_simulations[sim_id].axiom
    nb_iter = lsystem_simulations[sim_id].derivationLength #lire dans derivation_length #335 #30
    lsystem_simulations[sim_id].opt_external_coupling = 1 # met a un l'option external coupling

    # copie du lsystem pour récupérer la taille des voxels
    lsys_temp = lsystem_simulations[sim_id]
    lstring_temp = lsys_temp.derive(lstring, 0, 1)
    # récupère toutes les variables du lsystem
    tag_loop_inputs = lsystem_simulations[sim_id].tag_loop_inputs
    invar, outvar, invar_sc, ParamP, \
        station, carto, meteo_j, mng_j,  \
        DOY, cutNB, start_time, nbplantes,  \
        surfsolref, m_lais, dicFeuilBilanR,  \
        surf_refVOX, triplets, ls_dif, S, par_SN,  \
        lims_sol, ls_roots, stateEV, Uval, b_,  \
        ls_mat_res, vCC, ls_ftswStress, ls_NNIStress,  \
        ls_TStress, lsApex, lsApexAll, dicOrgans,  \
        deltaI_I0, nbI_I0, I_I0profilLfPlant, I_I0profilPetPlant,  \
        I_I0profilInPlant, NlClasses, NaClasses, NlinClasses,  \
        opt_stressW, opt_stressN, opt_stressGel, opt_residu, dxyz = tag_loop_inputs

    if active=="ratp" or passive=="ratp" :
        ## INITIALISATION LIGHTVEGEMANAGER
        environment = {}
        ratp_parameters_legume = {}
        ratp_parameters_plantgl = {}

        # Paramètres pré-simulation
        environment["names"] = ["l-egume"]
        environment["coordinates"] = [46.43,0,0] # latitude, longitude, timezone
        environment["sky"] = "turtle46" # ["file", "runscripts/legume/sky_5.data"] # "turtle46" # turtle à 46 directions par défaut
        environment["diffus"] = True
        environment["direct"] = False
        environment["reflected"] = False
        environment["reflectance coefficients"] = [[0., 0.]]
        environment["infinite"] = True

        if ratpgeo == "grid":
            ## PARAMETRES RATP scene grille l-egume ##
            ratp_parameters_legume["voxel size"] = [d*0.01 for d in dxyz]
            ratp_parameters_legume["soil reflectance"] = [0., 0.]
            ratp_parameters_legume["mu"] = [1.]
            ratp_parameters_legume["origin"] = [0., 0., 0.]

            lghtratp = LightVegeManager(environment=environment,
                                            lightmodel="ratp",
                                            lightmodel_parameters=ratp_parameters_legume,
                                            main_unit="m")
        elif ratpgeo == "plantgl":
            ## PARAMETRES RATP scene PlantGL ##
            ratp_parameters_plantgl["voxel size"] = [d*0.01 for d in dxyz]
            ratp_parameters_plantgl["xy max"] = lsys_temp.cote
            ratp_parameters_plantgl["soil reflectance"] = [0., 0.]
            ratp_parameters_plantgl["mu"] = [1.]
            ratp_parameters_plantgl["tesselation level"] = 2
            ratp_parameters_plantgl["angle distrib algo"] = "file"
            ratp_parameters_plantgl["angle distrib file"] = "runscripts/legume/luzerne_angle_distrib.data"

            lghtratp = LightVegeManager(environment=environment,
                                            lightmodel="ratp",
                                            lightmodel_parameters=ratp_parameters_plantgl,
                                            main_unit="m")

        l_epsi_passive = []

    # tableau de comparaison entre riri et ratp
    trans_min_diff_vox = []
    trans_max_diff_vox = []
    trans_mean_diff_vox = []
    trans_median_diff_vox = []
    trans_sd_diff_vox = []
    abs_min_diff_vox = []
    abs_max_diff_vox = []
    abs_mean_diff_vox = []
    abs_median_diff_vox = []
    abs_sd_diff_vox = []

    # début de la simulation
    for i in range(nb_iter+1):
        print('time step: ',i)

        if i == 90:
            print("spot")

        lstring = lsystem_simulations[sim_id].derive(lstring, i, 1)

        # récupère toutes les variables du lsystem
        tag_loop_inputs = lsystem_simulations[sim_id].tag_loop_inputs
        invar, outvar, invar_sc, ParamP, \
            station, carto, meteo_j, mng_j,  \
            DOY, cutNB, start_time, nbplantes,  \
            surfsolref, m_lais, dicFeuilBilanR,  \
            surf_refVOX, triplets, ls_dif, S, par_SN,  \
            lims_sol, ls_roots, stateEV, Uval, b_,  \
            ls_mat_res, vCC, ls_ftswStress, ls_NNIStress,  \
            ls_TStress, lsApex, lsApexAll, dicOrgans,  \
            deltaI_I0, nbI_I0, I_I0profilLfPlant, I_I0profilPetPlant,  \
            I_I0profilInPlant, NlClasses, NaClasses, NlinClasses,  \
            opt_stressW, opt_stressN, opt_stressGel, opt_residu, dxyz = tag_loop_inputs    

        ## Paramètres météo ## 
        doy = lsystem_simulations[sim_id].meteo["DOY"][i]
        hour = 12
        energy = 0.48*meteo_j['RG']*10000/(3600*24)#flux PAR journalier moyen en W.m-2 / RG en j.cm-2

        if active=="ratp" or passive=="ratp" :
            ## Paramètres scene ##
            # Transfert grille l-egume vers RATP
            if ratpgeo == "grid":
                # Surface foliaire et distribution des angles de l-egume
                scene_legume = {}
                scene_legume["LA"] = m_lais
                scene_legume["distrib"] = ls_dif

                geometry = {}
                geometry["scenes"] = [scene_legume]

             # # Scene PlantGL feuilles
            elif ratpgeo == "plantgl":
                # # récupère la scène PlantGL
                scene_plantgl = lsystem_simulations[sim_id].sceneInterpretation(lstring)

                geometry = {}
                geometry["scenes"] = [scene_plantgl]
                geometry["transformations"] = {}
                geometry["transformations"]["scenes unit"] = ["cm"]
                geometry["transformations"]["xyz orientation"] = ["y+ = y-"]

            lghtratp.init_scenes(geometry)
            lghtratp.run(PARi=energy, day=doy, hour=hour, truesolartime=True, parunit="RG")   

            if writegeo:
                if ratpgeo == "grid":
                    lghtratp.VTKinit(foldout+"grid_"+str(i)+"_",printtriangles=False)
                    if triangles == "y":
                        scene_plantgl = lsystem_simulations[sim_id].sceneInterpretation(lstring)
                        lghtratp.PlantGL_to_VTK(scene_plantgl, i, foldout+"ratpgrid_")

                elif ratpgeo == "plantgl":
                    lghtratp.VTKinit(foldout+"plantgl_"+str(i)+"_")
        
        ############
        # step light transfer coupling
        ############
        if active=="legume":
            # PAR / Blue voxel
            tag_light_inputs = [m_lais / surf_refVOX, triplets, ls_dif, energy * surf_refVOX]  # input tag

            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, optsky=station['optsky'], opt=station['sky'])

        if active=="ratp" or passive=="ratp":
            # transfert des sorties
            res_trans_2, res_abs_i_2 = lghtratp.to_l_egume(m_lais, energy)

            if active=="ratp":
                res_trans, res_abs_i = res_trans_2, res_abs_i_2
            elif passive=="ratp":
                # calcul du epsi
                nblignes = len(dicFeuilBilanR['nump']) #- 1
                temp_paraf = []
                if nblignes > 0:
                    for i in range(nblignes):
                        surf = dicFeuilBilanR['surf'][i]
                        vox = [dicFeuilBilanR['Vox0'][i], dicFeuilBilanR['Vox1'][i], dicFeuilBilanR['Vox2'][i]]
                        sVOX = m_lais[0][vox[2]][vox[1]][vox[0]]
                        paraF = res_abs_i_2[0][vox[2]][vox[1]][vox[0]] * surf / sVOX * 3600. * 24 / 1000000.
                        temp_paraf.append(paraF)
                # ls_epsi
                transmi_sol = np.sum(res_trans_2[-1][:][:]) / (energy * surfsolref)
                epsi = 1. - transmi_sol  # bon
                ls_epsi = epsi * np.array(temp_paraf)/ (np.sum(np.array(temp_paraf)) + np.sum(np.array(temp_paraf)) + 10e-15)
                l_epsi_passive.append(sum(ls_epsi))
                print('epsi ratp passive', sum(ls_epsi))

                # comparaison erreur relative
                trans_err_diff = []
                abs_err_diff = []
                for ix in range(res_trans.shape[2]):
                    for iy in range(res_trans.shape[1]):
                        for iz in range(res_trans.shape[0]):
                            if m_lais[0][iz][iy][ix] > 0:
                                trans_err_diff.append(100*(res_trans[iz][iy][ix] - res_trans_2[iz][iy][ix])/res_trans[iz][iy][ix])
                                abs_err_diff.append(100*(res_abs_i[0][iz][iy][ix] - res_abs_i_2[0][iz][iy][ix])/res_abs_i[0][iz][iy][ix])

                trans_min_diff_vox.append(np.min(trans_err_diff))
                trans_max_diff_vox.append(np.max(trans_err_diff))
                trans_mean_diff_vox.append(np.mean(trans_err_diff))
                trans_median_diff_vox.append(np.median(trans_err_diff))
                # trans_sd_diff_vox.append(np.sd(trans_err_diff))
                abs_min_diff_vox.append(np.min(abs_err_diff))
                abs_max_diff_vox.append(np.max(abs_err_diff))
                abs_mean_diff_vox.append(np.mean(abs_err_diff))
                abs_median_diff_vox.append(np.median(abs_err_diff))
                # abs_sd_diff_vox = []


        iteration_legume_withoutlighting(lsystem_simulations[sim_id], res_trans, res_abs_i, tag_loop_inputs, energy)

    print((''.join((sim_id, " - done"))))
    print("simulation time : ", time.time() - start, " s")

    # impression
    if passive=="ratp":
        df_ratp_passive = pd.DataFrame({
                                            'var' : ['epsi']*(len(l_epsi_passive)-2),
                                            'steps' : list(range(lsystem_simulations[sim_id].DOYdeb, lsystem_simulations[sim_id].DOYend-1)),
                                            'epsi' : l_epsi_passive[:-2]
                                        })
        df_ratp_passive.to_csv(foldout+"outputs_ratp_passive.csv", index=False)

        df_vox_compare = pd.DataFrame({
                                            'trans min' : trans_min_diff_vox,
                                            'trans max' : trans_max_diff_vox,
                                            'trans mean' : trans_mean_diff_vox,
                                            'trans median' : trans_median_diff_vox,
                                            'abs min' : abs_min_diff_vox,
                                            'abs max' : abs_max_diff_vox,
                                            'abs mean' : abs_mean_diff_vox,
                                            'abs median' : abs_median_diff_vox
                                        })
        df_vox_compare.to_csv(foldout+"diff_voxels_ratp_passive.csv", index=False)

    lsystem_simulations[sim_id].clear()


if __name__ == "__main__":
    #definition d'arguments avec getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i: o: a: p: r: w:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    # valeur par défaut
    foldin = "l-egume/legume/input/"
    foldout = "outputs/legume/"
    active = "legume" # legume ou ratp
    passive = "ratp" # legume ou ratp
    ratpgeo = "grid" # grid ou plantgl
    writegeo = "n" # "y" ou "no"
    triangles = "y"

    # récupère les arguments en entrée
    for opt, arg in opts:
        if opt in ("-i"):
            foldin = str(arg)
        elif opt in ("-o"):
            foldout = str(arg)
        elif opt in ("-a"):
            active = str(arg)
        elif opt in ("-p"):
            passive = str(arg)
        elif opt in ("-r"):
            ratpgeo = str(arg)
        elif opt in ("-w"):
            writegeo = str(arg)
    
    print("=== BEGIN ===")
    print("--- Simulation l-egume")
    if active=="legume":
        if passive=="ratp":
            if ratpgeo=="grid":
                print("=== === L-EGUME + LVM : RATP PASSIVE (GRID) === ===")
            elif ratpgeo=="plantgl":
                print("=== === L-EGUME + LVM : RATP PASSVE (PLANTGL) === ===")
        else:
            print("=== === L-EGUME  === ===")
        
    elif active=="ratp":
        if ratpgeo=="grid":
                print("=== === LVM : RATP  (GRID) === ===")
        elif ratpgeo=="plantgl":
            print("=== === LVM : RATP  (PLANTGL) === ===")
        
    if writegeo=="y":
        simulation(foldin, foldout, active, passive, ratpgeo, True, triangles)
    else : simulation(foldin, foldout, active, passive, ratpgeo)
    
    print("=== END ===")