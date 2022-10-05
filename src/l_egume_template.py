import os, sys

from openalea.lpy import *

import numpy as np

import legume
# import des modules à partir du chemin local de legume
path_ = os.path.dirname(os.path.abspath(legume.__file__))
sys.path.insert(0, path_)
import IOxls
import IOtable
import run_legume_usm as runl
import ShootMorpho as sh
import daily_loop as loop
import RIRI5 as riri

# prépare une liste de l-system à lancer
def initialisation_legume(foldin, foldout, fusms, ongletBatch):
    
    #lecture de la liste des usm
    mn_path = os.path.join(foldin, fusms)
    
    usms = IOxls.xlrd.open_workbook(mn_path)
    ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))

    # cree la liste de L-systems et liste des noms
    lsystem_simulations = {}
    names_simulations = []
    for i in range(len(ls_usms['ID_usm'])):
        if int(ls_usms['torun'][i]) == 1:  # si 1 dans la colonne 'torun' l'ajoute a la liste
            mylsys = runl.lsystemInputOutput_usm(fusms, foldin=foldin, ongletBatch=ongletBatch, i=i, path_OUT=foldout)
            name = list(mylsys)[0]
            names_simulations.append(name)
            lsystem_simulations[name] = mylsys[name]

    return lsystem_simulations, names_simulations



def iteration_legume_withoutlighting(iter, lsystem_simulations, names_simulations, 
                                        m_lais_simu, res_trans, res_abs_i, 
                                        meteo_j, surf_refVOX, surfsolref, 
                                        energy):

    # rassemble les paramètres propres à chaque lsystem
    list_invar, list_outvar, list_invar_sc, list_ParamP, \
    list_cutNB, list_nbplantes, list_start_time, \
    list_dicFeuilBilanR,  \
    list_par_SN,  \
    list_ls_roots, \
    list_ls_mat_res, list_vCC, list_ls_ftswStress, list_ls_NNIStress,  \
    list_ls_TStress, list_lsApex, list_lsApexAll, list_dicOrgans,  \
    list_deltaI_I0, list_nbI_I0, list_I_I0profilLfPlant, list_I_I0profilPetPlant,  \
    list_I_I0profilInPlant, list_NlClasses, list_NaClasses, list_NlinClasses,  \
    list_opt_stressW, list_opt_stressN = ([] for i in range(28))
    for n in names_simulations : 
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
        opt_stressW, opt_stressN, opt_stressGel, opt_residu, dxyz = lsystem_simulations[n].tag_loop_inputs

        list_invar.append(invar)
        list_dicFeuilBilanR.append(dicFeuilBilanR)
        list_ParamP.append(ParamP)
        list_outvar.append(outvar)
        list_invar_sc.append(invar_sc)
        list_cutNB.append(cutNB)
        list_vCC.append(vCC)
        list_start_time.append(start_time)
        list_nbplantes.append(nbplantes)
        list_ls_ftswStress.append(ls_ftswStress)
        list_ls_NNIStress.append(ls_NNIStress)
        list_ls_TStress.append(ls_TStress)
        list_ls_roots.append(ls_roots)
        list_lsApex.append(lsApex)
        list_lsApexAll.append(lsApexAll)
        list_opt_stressW.append(opt_stressW)
        list_opt_stressN.append(opt_stressN)
        list_par_SN.append(par_SN)
        list_I_I0profilLfPlant.append(I_I0profilLfPlant)
        list_I_I0profilPetPlant.append(I_I0profilPetPlant)
        list_I_I0profilInPlant.append(I_I0profilInPlant)
        list_NlClasses.append(NlClasses)
        list_NaClasses.append(NaClasses)
        list_NlinClasses.append(NlinClasses)
        list_dicOrgans.append(dicOrgans)
        list_deltaI_I0.append(deltaI_I0)
        list_nbI_I0.append(nbI_I0)
        list_ls_mat_res.append(ls_mat_res)
    
    # pour les variables communes
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
    opt_stressW, opt_stressN, opt_stressGel, opt_residu, dxyz = lsystem_simulations[names_simulations[0]].tag_loop_inputs
    
    m_lais = m_lais_simu
    # R_FR voxel (calcul de zeta)
    tag_light_inputs2 = [res_trans / (energy * surf_refVOX)]  # input tag
    # tag_light_inputs2 = [res_trans]  # input tag
    res_rfr = riri.rfr_calc_relatif(*tag_light_inputs2)

    # ls_epsi
    transmi_sol = np.sum(res_trans[-1][:][:]) / (energy * surfsolref)
    epsi = 1. - transmi_sol  # bon

    # calul des interception feuille et ls_epsi plante
    list_ls_epsi = []
    for k in range(len(names_simulations)) :        
        list_dicFeuilBilanR[k] = sh.calc_paraF(list_dicFeuilBilanR[k], m_lais, res_abs_i, force_id_grid = k)
        sh.calc_para_Plt(list_invar[k], list_dicFeuilBilanR[k])

        list_ls_epsi.append(epsi * list_invar[k]['parip'] / (np.sum(list_invar[k]['parip']) + np.sum(list_invar[k]['parip']) + 10e-15))
        print('main', names_simulations[k], 'epsi', sum(list_ls_epsi[-1]))

    ##########
    # Step Potential plant growth
    ##########
    list_ls_demandeN_bis = []
    list_temps = []
    for k in range(len(names_simulations)) : 
        list_invar[k], list_outvar[k], ls_demandeN_bis, temps = loop.daily_growth_loop(list_ParamP[k], list_invar[k], list_outvar[k], list_ls_epsi[k], meteo_j, mng_j,
                                                                        list_nbplantes[k], surfsolref, list_ls_ftswStress[k],
                                                                        list_ls_NNIStress[k], list_ls_TStress[k], list_lsApex[k], list_lsApexAll[k],
                                                                        list_opt_stressW[k], list_opt_stressN[k], opt_stressGel)

        list_ls_demandeN_bis.append(ls_demandeN_bis)
        list_temps.append(temps)

    ##########
    # step soil
    ##########

    # gere l'aggregation des entrees par plante
    list_nb = [len(list_ls_epsi[0].tolist())]
    ls_epsi = list_ls_epsi[0].tolist()
    ls_demandeN_bis = list_ls_demandeN_bis[0].tolist()
    ls_roots = list_ls_roots[0]
    ParamP = list_ParamP[0]
    for k in range(1,len(names_simulations)) : 
            list_nb.append(len(list_ls_epsi[k].tolist()))
            ls_epsi = ls_epsi + list_ls_epsi[k].tolist()
            ls_demandeN_bis = ls_demandeN_bis + list_ls_demandeN_bis[k].tolist()
            ls_roots = ls_roots + list_ls_roots[k]
            ParamP = ParamP + list_ParamP[k]

    
    # step soil en commun
    tag_inputs_soil_step = [S, par_SN, surfsolref, stateEV, Uval, b_, meteo_j, mng_j, ParamP, ls_epsi, ls_roots, ls_demandeN_bis, opt_residu]  # input tag
    res_soil_step = loop.step_bilanWN_sol(*tag_inputs_soil_step)
    S, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol = res_soil_step

    # gere desagregattion par esp des sorties
    list_ls_ftsw = []
    list_ls_transp = []
    list_ls_Act_Nuptake_plt = []
    list_temps_sol = []
    for k in range(len(names_simulations)) :
        a = sum(list_nb[:k])
        b = sum(list_nb[:k+1])

        list_ls_ftsw.append(ls_ftsw[a:b])
        list_ls_transp.append(ls_transp[a:b])
        list_ls_Act_Nuptake_plt.append(ls_Act_Nuptake_plt[a:b])
        list_temps_sol.append(temps_sol[a:b])

    ##########
    # setp update plant stress variables
    ##########
    for k in range(len(names_simulations)) :
        tag_inputs_stress = [list_ParamP[k], list_invar[k], list_invar_sc[k], list_temps[k], DOY, list_nbplantes[k], surfsolref, list_ls_epsi[k], list_ls_ftsw[k], list_ls_transp[k],
                                list_ls_Act_Nuptake_plt[k], list_ls_demandeN_bis[k], list_ls_ftswStress[k], list_ls_TStress[k], list_dicOrgans[k], list_dicFeuilBilanR[k], list_lsApex[k],
                                list_start_time[k], list_cutNB[k], list_deltaI_I0[k], list_nbI_I0[k], list_I_I0profilLfPlant[k], list_I_I0profilPetPlant[k],
                                list_I_I0profilInPlant[k], list_NlClasses[k], list_NaClasses[k], list_NlinClasses[k], list_outvar[k]]

        list_invar[k], list_invar_sc[k], list_outvar[k], list_I_I0profilInPlant[k], list_ls_ftswStress[k], list_ls_NNIStress[k], list_ls_TStress[k] = loop.Update_stress_loop(*tag_inputs_stress)
    
    ##########
    # step update soil residues senescence
    ##########

    ls_mat_res = list_ls_mat_res[0]

    # refait initialisation des residues au step 1 avec ensemble des plante (ParamP commun)
    if iter == 1 and opt_residu == 1:
        CC = sh.init_plant_residues_fromParamP(S, opt_residu, ParamP)

    if opt_residu == 1:  # option residu activee: mise a jour des cres
         # gere l'aggregation des entrees par espce
        vCC = list_vCC[0]
        for v in list_vCC[1:] : vCC = vCC + v
        
        invar_merge = list_invar[0]
        for inv in list_invar[1:] :
            invar_merge = IOtable.merge_dict_list([invar_merge, inv], ls_keys=['dMSenRoot', 
                                                                                'dMSenFeuil', 
                                                                                'dMSenTige', 
                                                                                'dMSenNonRec', 
                                                                                'dMSenPiv', 
                                                                                'dMSmortGel_aer', 
                                                                                'dMSmortPlant_aer', 
                                                                                'dMSmortPlant_racfine', 
                                                                                'dMSmortPlant_pivot', 
                                                                                'isGelDam'])

        tag_inputs_residue_updt = [ls_mat_res, S, ls_roots, par_SN['PROFHUMs'], ParamP, invar_merge, opt_stressGel]
        ls_mat_res = loop.distrib_residue_mat_frominvar(*tag_inputs_residue_updt) #update la matrice des residus (propre a l-egume/VGL)
        S = loop.merge_residue_mat(ls_mat_res, vCC, S) #update du sol

    #########
    # reinjecte les sorties midiee dans le lsystem
    #########
    for k,n in enumerate(names_simulations):
        lsystem_simulations[n].invar = list_invar[k]
        lsystem_simulations[n].outvar = list_outvar[k]
        lsystem_simulations[n].invar_sc = list_invar_sc[k]

        lsystem_simulations[n].S = S
        lsystem_simulations[n].stateEV = stateEV
        lsystem_simulations[n].ls_mat_res = ls_mat_res

        lsystem_simulations[n].res_trans = res_trans
        lsystem_simulations[n].res_abs_i = np.array([res_abs_i[k]])
        lsystem_simulations[n].res_rfr = res_rfr

        lsystem_simulations[n].ls_ftswStress = list_ls_ftswStress[k]
        lsystem_simulations[n].ls_NNIStress = list_ls_NNIStress[k]
        lsystem_simulations[n].ls_TStress = list_ls_TStress[k]
        lsystem_simulations[n].I_I0profilInPlant = list_I_I0profilInPlant[k]
