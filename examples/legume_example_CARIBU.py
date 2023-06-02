import sys
import os
import getopt
import time

import numpy as np
import pandas as pd

try :
    from  LightVegeManager import *
    from examples.l_egume_template import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
    from  LightVegeManager import *
    from examples.l_egume_template import *

'''
Comparaison de la lumière entre l-egume (RiRi) et CARIBU

'''

def simulation(foldin, foldout, active, passive, writegeo=False):
    globalstart = time.time()

    # nom du test à lancer    
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

    if active=="caribu" or passive=="caribu" :
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
        caribu_parameters = {}

        # Paramètres pré-simulation
        environment["names"] = ["l-egume"]
        environment["coordinates"] = [46.43,0,0] # latitude, longitude, timezone
        environment["sky"] = "turtle46" # [4, 5, "soc"] 
        environment["diffus"] = True
        environment["direct"] = False
        environment["reflected"] = False
        environment["infinite"] = True

        environment["caribu opt"] = {} 
        environment["caribu opt"]["par"] = (0.10, 0.07)
       
        caribu_parameters["sun algo"] = "caribu"
        nxyz = [m_lais.shape[3], m_lais.shape[2], m_lais.shape[1]]
        dxyz = [x * 0.01 for x in dxyz] # conversion de cm à m
        orig = [0.,0.,0.]
        caribu_parameters["sensors"] = ["grid", dxyz, nxyz, orig, foldout, "vtk"]
        caribu_parameters["debug"] = False

        lghtcaribu = LightVegeManager(environment=environment,
                                    lightmodel="caribu",
                                    lightmodel_parameters=caribu_parameters, 
                                    main_unit="m")

        # tableaux de sorties spécifiques 
        epsi_passive = []
        para_passive = []
        for k,n in enumerate(names_simulations) :
            lstring_temp = lsystem_simulations[n].derive(lstring[k], 0, 1)
            epsi_passive.append([[] for i in range(lsystem_simulations[n].nbplantes)])
            para_passive.append([[] for i in range(lsystem_simulations[n].nbplantes)])

        # temps de calcul
        time_legume = 0
        time_ratp_tot = 0
        time_ratp_run = 0
        step_time_ratp = []
        step_time_ratp_run = []
        step_time_leg = []


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

        if active=="caribu" or passive=="caribu" :
            ## Paramètres scene ##     
            scenes_plantgl = []
            if i > 0 :      
                # récupère la scène PlantGL
                for k, n in enumerate(names_simulations):
                    scenes_plantgl.append(lsystem_simulations[n].sceneInterpretation(lstring[k]))

            geometry = {}
            geometry["scenes"] = scenes_plantgl
            # domaine dans l'unité de la scène finale (ici en m)
            geometry["domain"] = ((0,0), (m_lais.shape[3] * dxyz[0] * 0.01, m_lais.shape[2] * dxyz[1] * 0.01))
            geometry["transformations"] = {}
            geometry["transformations"]["scenes unit"] = ["cm"] * len(names_simulations) # ne concerne que geometry["scenes"]

            start=time.time()
            lghtcaribu.init_scenes(geometry)
            lghtcaribu.run(energy=energy, day=doy, hour=hour, truesolartime=True, parunit="RG")
            
            t_ratp_tot = (time.time() - start)
            time_ratp_tot += t_ratp_tot
            time_ratp_run += lghtcaribu.modelruntime

            step_time_ratp.append(t_ratp_tot)
            step_time_ratp_run.append(lghtcaribu.modelruntime)

            if writegeo:
                lghtcaribu.VTKout(foldout+"triangle",iteration=i, triangles=True, voxels=False, outvariables=["par Eabs", "par Ei"])

        
        ############
        # step light transfer coupling
        ############
        if active=="legume":
            # PAR / Blue voxel
            tag_light_inputs = [m_lais / surf_refVOX, triplets, ls_dif, energy * surf_refVOX]  # inp+ut tag

            start=time.time()
            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, optsky=station['optsky'], opt=station['sky'])

            t_legume = (time.time() - start)
            time_legume += t_legume
            step_time_leg.append(t_legume)
            
        if active=="caribu" or passive=="caribu":
            # rassemble les paramètres propres à chaque lsystem
            list_invar2, list_dicFeuilBilanR = [], []
            for n in names_simulations : 
                list_invar2.append(lsystem_simulations[n].tag_loop_inputs[0])
                list_dicFeuilBilanR.append(lsystem_simulations[n].tag_loop_inputs[14])

            # transfert des sorties
            res_trans_2 = lghtcaribu.to_l_egume(m_lais=m_lais,
                                                energy=energy,
                                                list_lstring = lstring, 
                                                list_dicFeuilBilanR=list_dicFeuilBilanR, 
                                                list_invar=list_invar2)
                  
            if active=="caribu":
                list_invar = list_invar2
                res_trans = res_trans_2
                res_abs_i = None

            # calcul du epsi sur les résultats de RATP (non utilisé dans la simulation)
            if passive=="caribu":
                list_invar = None

                transmi_sol = np.sum(res_trans_2[-1][:][:]) / (energy * surfsolref)
                epsi = 1. - transmi_sol
                # calul des interception feuille et ls_epsi plante
                pari_canopy = 0
                for k in range(len(names_simulations)) :
                    pari_canopy += np.sum(list_invar2[k]['parip'])
                
                ls_epsi = []
                for k in range(len(names_simulations)) :
                    ls_epsi.append(epsi * list_invar2[k]['parip'] / (pari_canopy + 10e-15))
                    print('CARIBU passive',names_simulations[k],'epsi', sum(ls_epsi[-1]))

                # Enregistre les outputs
                for k in range(len(names_simulations)) : 
                    for p in range(lsystem_simulations[names_simulations[k]].nbplantes):
                        epsi_passive[k][p].append(ls_epsi[k][p])
                        para_passive[k][p].append(list_invar2[k]['parip'][p])

        iteration_legume_withoutlighting(i, lsystem_simulations, names_simulations, 
                                            meteo_j, energy, surf_refVOX, surfsolref, 
                                            m_lais, res_trans, res_abs_i, list_invar)

    for n in names_simulations : print((''.join((n, " - done"))))
    print("simulation time : ", time.time() - globalstart, " s")

    # impression
    if passive=="caribu":
        deb = lsystem_simulations[names_simulations[0]].DOYdeb
        fin = lsystem_simulations[names_simulations[0]].DOYend-1
        for k,n in enumerate(names_simulations):
            dict_ratp_passive = {
                                    'var' : ['epsi']*(len(epsi_passive[k][0])-2) + ['PARiPlante']*(len(epsi_passive[k][0])-2),
                                    'steps' : list(range(deb, fin)) + list(range(deb, fin))
                                }
            for p in range(lsystem_simulations[names_simulations[k]].nbplantes) :
                dict_ratp_passive['plante'+str(p)] = epsi_passive[k][p][0:fin-deb] + para_passive[k][p][0:fin-deb]
            pd.DataFrame(dict_ratp_passive).to_csv(foldout+"outputs_ratp_passive_"+str(n)+".csv", index=False)

        pd.DataFrame({
            "legume" : step_time_leg, 
            "caribu run" : step_time_ratp_run, 
            "caribu tot" : step_time_ratp                   
            }).to_csv(foldout+"cputime_steps.csv", index=False)

    pd.DataFrame({
                    "legume" : [time_legume], 
                    "caribu run" : [time_ratp_run], 
                    "caribu tot" : [time_ratp_tot]
                    }).to_csv(foldout+"global_outputs.csv", index=False)


    # désallocation des lsystem
    for n in names_simulations : lsystem_simulations[n].clear()


if __name__ == "__main__":
    #definition d'arguments avec getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i: o: a: p: w:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    # valeur par défaut
    foldin = "l-egume/legume/input/"
    foldout = "outputs/legume_caribu/photomorph_ramif_3shoots_cari_pass/"
    active = "legume" # legume ou caribu
    passive = "caribu" # legume ou caribu
    writegeo = True

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
        elif opt in ("-w"):
            writegeo_str = str(arg)
            if writegeo_str == "y" : writegeo = True
            elif writegeo_str == "n" : writegeo = False

    print("foldin : ",foldin)
    print("foldout : ",foldout)
    print("active : ",active)
    print("passive : ",passive)
    print("writegeo : ",writegeo)

    simulation(foldin, foldout, active, passive, writegeo)
    
    print("=== END ===")