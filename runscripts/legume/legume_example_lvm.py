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

def simulation(foldin, foldout, active, passive, ratpgeo, skytype=2, writegeo=False,triangles="n"):
    
    
    fusms = "liste_usms_exemple.xls"
    ongletBatch = "exemple"

    # prépare une liste de l-system à lancer
    lsystem_simulations, names_simulations = initialisation_legume(foldin, foldout, fusms, ongletBatch)

    # initialisation pour chaque lsystem
    nb_iter = lsystem_simulations[names_simulations[0]].derivationLength
    lstring = []
    for n in names_simulations :
        lstring.append(lsystem_simulations[n].axiom)
        #lire dans derivation_length #335 #30
        lsystem_simulations[n].opt_external_coupling = 1 # met a un l'option external coupling

    if active=="ratp" or passive=="ratp" :
        # copie du lsystem pour récupérer la taille des voxels
        lsys_temp = lsystem_simulations[names_simulations[0]]
        lstring_temp = lsys_temp.derive(lstring[0], 0, 1)
        # récupère toutes les variables du lsystem
        tag_loop_inputs = lsystem_simulations[names_simulations[0]].tag_loop_inputs
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
        ## INITIALISATION LIGHTVEGEMANAGER
        environment = {}
        ratp_parameters_legume = {}
        ratp_parameters_plantgl = {}

        # nombre d'entité à prendre en compte
        nent = len(names_simulations)

        # Paramètres pré-simulation
        environment["names"] = ["l-egume"]
        environment["coordinates"] = [46.43,0,0] # latitude, longitude, timezone
        if skytype == 1 : environment["sky"] = ["file", "runscripts/legume/sky_5.data"]
        elif skytype == 2 : environment["sky"] = "turtle46"
        elif skytype == 3 : environment["sky"] = [10, 10, "soc"]
        environment["diffus"] = True
        environment["direct"] = False
        environment["reflected"] = False
        environment["reflectance coefficients"] = [[0., 0.]] * nent
        environment["infinite"] = True

        if ratpgeo == "grid":
            ## PARAMETRES RATP scene grille l-egume ##
            ratp_parameters_legume["voxel size"] = [d*0.01 for d in dxyz] # en m
            ratp_parameters_legume["soil reflectance"] = [0., 0.]
            ratp_parameters_legume["mu"] = [1.] * nent
            ratp_parameters_legume["origin"] = [0., 0., 0.]
            # ratp_parameters_legume["transmitted"] = True

            lghtratp = LightVegeManager(environment=environment,
                                            lightmodel="ratp",
                                            lightmodel_parameters=ratp_parameters_legume,
                                            main_unit="m")
        # elif ratpgeo == "plantgl":
        #     ## PARAMETRES RATP scene PlantGL ##
        #     ratp_parameters_plantgl["voxel size"] = [d*0.01 for d in dxyz]
        #     ratp_parameters_plantgl["number voxels"] = [m_lais.shape[3], m_lais.shape[2], m_lais.shape[1]]
        #     ratp_parameters_legume["origin"] = [0., 0., 0.]
        #     ratp_parameters_plantgl["xy max"] = lsys_temp.cote
        #     ratp_parameters_plantgl["soil reflectance"] = [0., 0.]
        #     ratp_parameters_plantgl["mu"] = [1.]
        #     ratp_parameters_plantgl["tesselation level"] = 2
        #     ratp_parameters_plantgl["angle distrib algo"] = "file"
        #     ratp_parameters_plantgl["angle distrib file"] = "runscripts/legume/luzerne_angle_distrib.data"

        #     lghtratp = LightVegeManager(environment=environment,
        #                                     lightmodel="ratp",
        #                                     lightmodel_parameters=ratp_parameters_plantgl,
        #                                     main_unit="m")
        
        # tableaux de sorties spécifiques 
        epsi_passive = []
        para_passive = []
        for k,n in enumerate(names_simulations) :
            lstring_temp = lsystem_simulations[n].derive(lstring[k], 0, 1)
            epsi_passive.append([[] for i in range(lsystem_simulations[n].nbplantes)])
            para_passive.append([[] for i in range(lsystem_simulations[n].nbplantes)])
        diff_voxel_para = []
        diff_voxel_part = []
        surf_vox = []

        # temps de calcul
        time_legume = 0
        time_ratp_tot = 0
        time_ratp_run = 0
        step_time_ratp = []
        step_time_ratp_run = []
        step_time_leg = []

        # somme respectifs
        S_para_legume = 0
        S_part_legume = 0
        S_para_ratp = 0
        S_part_ratp = 0

    if active == "riri":
        # copie du lsystem pour récupérer la taille des voxels
        lsys_temp = lsystem_simulations[names_simulations[0]]
        lstring_temp = lsys_temp.derive(lstring[0], 0, 1)
        # récupère toutes les variables du lsystem
        tag_loop_inputs = lsystem_simulations[names_simulations[0]].tag_loop_inputs
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
        ## INITIALISATION LIGHTVEGEMANAGER
        environment = {}
        riri_parameters_legume = {}

        # nombre d'entité à prendre en compte
        nent = len(names_simulations)

        # Paramètres pré-simulation
        environment["names"] = ["l-egume"]
        environment["sky"] = [station['optsky'], station['sky']]

        ## PARAMETRES RATP scene grille l-egume ##
        riri_parameters_legume["voxel size"] = [d*0.01 for d in dxyz] # en m

        lghtriri = LightVegeManager(environment=environment,
                                        lightmodel="riri",
                                        lightmodel_parameters=ratp_parameters_legume,
                                        main_unit="m")

    # début de la simulation
    for i in range(nb_iter+1):
        print('time step: ',i)

        for k, n in enumerate(names_simulations):
            lstring[k] = lsystem_simulations[n].derive(lstring[k], i, 1)

        # récupère toutes les variables du lsystem (on retient seulement les variables communes)
        tag_loop_inputs = lsystem_simulations[names_simulations[0]].tag_loop_inputs
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
        
        if len(names_simulations) > 1 :
                    #gere difference de dsitib par especes
            list_ls_dif = []
            list_m_lais = []
            for n in names_simulations:
                list_ls_dif.append(lsystem_simulations[n].tag_loop_inputs[17])
                list_m_lais.append(lsystem_simulations[n].tag_loop_inputs[13])
            
            same_entity = False
            
            for dis in list_ls_dif[1:]:
                same_entity = (list_ls_dif[0][0] == dis[0]).all()
            if same_entity : 
                ls_dif = list_ls_dif[0][0]
                m_lais = list_m_lais[0]
                for m in list_m_lais[1:]: m_lais = m_lais + m
            else:
                ls_dif = list_ls_dif[0]
                for d in list_ls_dif[1:]: ls_dif = ls_dif + d
                m_lais = np.array([m[0] for m in list_m_lais])

        ## Paramètres météo ## 
        doy = DOY
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
            # elif ratpgeo == "plantgl":
            #     # # récupère la scène PlantGL
            #     scene_plantgl = lsystem_simulations[names_simulations[0]].sceneInterpretation(lstring)

            #     geometry = {}
            #     geometry["scenes"] = [scene_plantgl]
            #     geometry["transformations"] = {}
            #     geometry["transformations"]["scenes unit"] = ["cm"]
            #     geometry["transformations"]["xyz orientation"] = ["y+ = y-"]

            start=time.time()
            lghtratp.init_scenes(geometry)
            lghtratp.run(PARi=energy, day=doy, hour=hour, truesolartime=True, parunit="RG")
            t_ratp_tot = (time.time() - start)
            time_ratp_tot += t_ratp_tot
            time_ratp_run += lghtratp.modelruntime

            step_time_ratp.append(t_ratp_tot)
            step_time_ratp_run.append(lghtratp.modelruntime)

            if writegeo:
                if ratpgeo == "grid":
                    lghtratp.VTKout(foldout+"grid",iteration=i, triangles=False, voxels=True, outvariables=["intercepted", "transmitted"])
                    if triangles == "y":
                        scenes_plantgl = [lsystem_simulations[names_simulations[0]].sceneInterpretation(lstring[0])]
                        for k,n in enumerate(names_simulations[1:]) : scenes_plantgl.append(lsystem_simulations[n].sceneInterpretation(lstring[k]))
                        lghtratp.PlantGL_to_VTK(scenes_plantgl, i, foldout+"ratpgrid_", "cm")

                elif ratpgeo == "plantgl":
                    lghtratp.VTKinit(foldout+"plantgl_"+str(i)+"_")
        
        ############
        # step light transfer coupling
        ############
        if active=="legume":
            # PAR / Blue voxel
            tag_light_inputs = [m_lais / surf_refVOX, triplets, ls_dif, energy * surf_refVOX]  # input tag

            start=time.time()
            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, optsky=station['optsky'], opt=station['sky'])
            t_legume = (time.time() - start)
            time_legume += t_legume
            step_time_leg.append(t_legume)
            
        elif active=="riri":
            scene_legume = {}
            scene_legume["LA"] = m_lais
            scene_legume["distrib"] = ls_dif
            scene_legume["triplets"] = triplets

            geometry = {}
            geometry["scenes"] = [scene_legume]

            lghtriri.run(PARi=energy, day=doy, hour=hour, truesolartime=True, parunit="RG")

            res_trans, res_abs_i = lghtriri.legume_transmitted_light, lghtriri.legume_intercepted_light

        if active=="ratp" or passive=="ratp":
            # transfert des sorties
            res_trans_2, res_abs_i_2 = lghtratp.to_l_egume(m_lais, energy)

            if active=="ratp":
                res_trans, res_abs_i = res_trans_2, res_abs_i_2
            
            # calcul du epsi sur les résultats de RATP (non utilisé dans la simulation)
            elif passive=="ratp":
                # ls_epsi
                transmi_sol = np.sum(res_trans_2[-1][:][:]) / (energy * surfsolref)
                epsi = 1. - transmi_sol  # bon

                # calcul paramètres par entité
                list_diff_para=[]
                surf_ent=[]
                for k in range(len(names_simulations)) :  
                    dicFeuilBilanR = lsystem_simulations[names_simulations[k]].tag_loop_inputs[14]    
                    invar_temp =   lsystem_simulations[names_simulations[k]].tag_loop_inputs[0]   
                    dicFeuilBilanR = sh.calc_paraF(dicFeuilBilanR, m_lais, res_abs_i_2, force_id_grid = k)
                    sh.calc_para_Plt(invar_temp, dicFeuilBilanR)
                    ls_epsi = epsi * invar_temp['parip'] / (np.sum(invar_temp['parip']) + np.sum(invar_temp['parip']) + 10e-15)
                    print('RATP passive',names_simulations[k],'epsi', sum(ls_epsi))
                    
                    for p in range(lsystem_simulations[names_simulations[k]].nbplantes):
                        epsi_passive[k][p].append(ls_epsi[p])
                        para_passive[k][p].append(invar_temp['parip'][p])       

                    diffpara_ent = []
                    list_surf_vox = []
                    for iz in range(lghtratp.legume_empty_layers, m_lais.shape[1]):
                        layer_diffpara_ent = []
                        layer_surf_vox = []
                        for ix in range(m_lais.shape[3]):
                            for iy in range(m_lais.shape[2]):
                                layer_diffpara_ent.append((res_abs_i[k][iz][iy][ix]-res_abs_i_2[k][iz][iy][ix]))
                                layer_surf_vox.append(m_lais[k][iz][iy][ix])

                                S_para_legume += res_abs_i[k][iz][iy][ix]
                                S_para_ratp += res_abs_i_2[k][iz][iy][ix]
                        diffpara_ent.append(layer_diffpara_ent)
                        list_surf_vox.append(layer_surf_vox)
                    list_diff_para.append(diffpara_ent)
                    surf_ent.append(list_surf_vox)
                diffpart=[]
                for iz in range(lghtratp.legume_empty_layers, m_lais.shape[1]):
                    layer_diffpart = []
                    for ix in range(m_lais.shape[3]):
                        for iy in range(m_lais.shape[2]):
                            layer_diffpart.append((res_trans[iz][iy][ix]-res_trans_2[iz][iy][ix]))
                            S_part_legume += res_trans[iz][iy][ix]
                            S_part_ratp += res_trans_2[iz][iy][ix]
                    diffpart.append(layer_diffpart)

                diff_voxel_para.append(list_diff_para)
                diff_voxel_part.append(diffpart)
                surf_vox.append(surf_ent)

        iteration_legume_withoutlighting(i, lsystem_simulations, names_simulations, 
                                            m_lais, res_trans, res_abs_i, 
                                            meteo_j, surf_refVOX, surfsolref, 
                                            energy)

    for n in names_simulations : print((''.join((n, " - done"))))
    print("simulation time : ", time.time() - start, " s")

    # impression
    if passive=="ratp":
        deb = lsystem_simulations[names_simulations[0]].DOYdeb
        fin = lsystem_simulations[names_simulations[0]].DOYend-1
        for k,n in enumerate(names_simulations):
            dict_ratp_passive = {
                                    'var' : ['epsi']*(len(epsi_passive[k][0])-2) + ['PARi']*(len(epsi_passive[k][0])-2),
                                    'steps' : list(range(deb, fin)) + list(range(deb, fin))
                                }
            for p in range(lsystem_simulations[names_simulations[k]].nbplantes) :
                dict_ratp_passive['plante'+str(p)] = epsi_passive[k][p][0:fin-deb] + para_passive[k][p][0:fin-deb]
            pd.DataFrame(dict_ratp_passive).to_csv(foldout+"outputs_ratp_passive_"+str(n)+".csv", index=False)

        for i in range(len(diff_voxel_para)):
            dic_part = {}
            for k,n in enumerate(names_simulations):
                dic_para = {}
                dic_surf = {}
                for j in range(len(diff_voxel_para[i][k])) :
                    dic_para["Layer"+str(j)] = diff_voxel_para[i][k][j]
                    dic_surf["Layer"+str(j)] = surf_vox[i][k][j]
                pd.DataFrame(dic_para).to_csv(foldout+"diff_para_"+str(n)+"_"+str(i)+".csv", index=False)
                pd.DataFrame(dic_surf).to_csv(foldout+"surf_vox_"+str(n)+"_"+str(i)+".csv", index=False)
            for j in range(len(diff_voxel_para[i][k])) :
                dic_part["Layer"+str(j)] = diff_voxel_part[i][j]
            pd.DataFrame(dic_part).to_csv(foldout+"diff_part_"+str(i)+".csv", index=False)
        pd.DataFrame({
            "legume" : step_time_leg, 
            "ratp run" : step_time_ratp_run, 
            "ratp tot" : step_time_ratp                   
            }).to_csv(foldout+"cputime_steps.csv", index=False)

    pd.DataFrame({
                    "legume" : [time_legume], 
                    "ratp run" : [time_ratp_run], 
                    "ratp tot" : [time_ratp_tot],
                    "PARa l-egume" : [S_para_legume],
                    "PARt l-egume" : [S_part_legume],
                    "PARa RATP" : [S_para_ratp],
                    "PARt RATP" : [S_part_ratp]
                    }).to_csv(foldout+"global_outputs.csv", index=False)


    # désallocation des lsystem
    for n in names_simulations : lsystem_simulations[n].clear()


if __name__ == "__main__":
    #definition d'arguments avec getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i: o: a: p: r: w: s:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    # valeur par défaut
    foldin = "l-egume/legume/input/"
    foldout = "outputs/legume_ratp/photomorpho_ramif_leg_sky100_1ent/"
    active = "ratp" # legume ou ratp
    passive = "legume" # legume ou ratp
    ratpgeo = "grid" # grid ou plantgl
    skytype = 3 # type de ciel : 1=sky5, 2=sky46, 3=sky100
    writegeo = "n" # "y" ou "no"
    triangles = "y"

    # # récupère les arguments en entrée
    # for opt, arg in opts:
    #     if opt in ("-i"):
    #         foldin = str(arg)
    #     elif opt in ("-o"):
    #         foldout = str(arg)
    #     elif opt in ("-a"):
    #         active = str(arg)
    #     elif opt in ("-p"):
    #         passive = str(arg)
    #     elif opt in ("-r"):
    #         ratpgeo = str(arg)
    #     elif opt in ("-w"):
    #         writegeo = str(arg)
    #     elif opt in ("-s"):
    #         skytype = str(arg)
    
    # print("=== BEGIN ===")
    # print("--- Simulation l-egume")
    # if active=="legume":
    #     if passive=="ratp":
    #         if ratpgeo=="grid":
    #             print("=== === L-EGUME + LVM : RATP PASSIVE (GRID) === ===")
    #         elif ratpgeo=="plantgl":
    #             print("=== === L-EGUME + LVM : RATP PASSVE (PLANTGL) === ===")
    #     else:
    #         print("=== === L-EGUME  === ===")
        
    # elif active=="ratp":
    #     if ratpgeo=="grid":
    #             print("=== === LVM : RATP  (GRID) === ===")
    #     elif ratpgeo=="plantgl":
    #         print("=== === LVM : RATP  (PLANTGL) === ===")
        
    # if writegeo=="y":
    #     simulation(foldin, foldout, active, passive, ratpgeo, skytype, True, triangles)
    # else : simulation(foldin, foldout, active, passive, ratpgeo, skytype)


    foldout = "outputs/legume_ratp/photomorpho_noramif_ratp_skt/"
    skytype = 1 # type de ciel : 1=sky5, 2=sky46, 3=sky100
    simulation(foldin, foldout, active, passive, ratpgeo, skytype)

    # foldout = "outputs/legume_ratp/nophotomorpho_ramif_ratp_sky46_vtk/"
    # skytype = 2 # type de ciel : 1=sky5, 2=sky46, 3=sky100
    # simulation(foldin, foldout, active, passive, ratpgeo, skytype)

    # foldout = "outputs/legume_ratp/nophotomorpho_ramif_ratp_sky100/"
    # skytype = 3 # type de ciel : 1=sky5, 2=sky46, 3=sky100
    # simulation(foldin, foldout, active, passive, ratpgeo, skytype)

    
    print("=== END ===")