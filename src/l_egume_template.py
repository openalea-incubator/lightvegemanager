import os, sys

from openalea.lpy import *

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



def iteration_legume_withoutlighting(lsystem, res_trans, res_abs_i, tag_loop_inputs):

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
        opt_stressW, opt_stressN, opt_stressGel, opt_residu = tag_loop_inputs
    
    # R_FR voxel (calcul de zeta)
    tag_light_inputs2 = [res_trans / (meteo_j['I0'] * surf_refVOX)]  # input tag
    # tag_light_inputs2 = [res_trans]  # input tag
    local_res_rfr = riri.rfr_calc_relatif(*tag_light_inputs2)

    res_rfr = local_res_rfr  # mise a jour variables globales

    # calul des interception feuille et ls_epsi plante
    dicFeuilBilanR = sh.calc_paraF(dicFeuilBilanR, m_lais, res_abs_i)
    ls_epsi, invar = loop.step_epsi(invar, res_trans, dicFeuilBilanR, meteo_j, surfsolref)

    ##########
    # Step Potential plant growth
    ##########

    invar, outvar, ls_demandeN_bis, temps = loop.daily_growth_loop(ParamP, invar, outvar, ls_epsi, meteo_j, mng_j,
                                                                    nbplantes, surfsolref, ls_ftswStress,
                                                                    ls_NNIStress, ls_TStress, lsApex, lsApexAll,
                                                                    opt_stressW, opt_stressN, opt_stressGel)

    ##########
    # step soil
    ##########

    tag_inputs_soil_step = [S, par_SN, surfsolref, stateEV, Uval, b_, meteo_j, mng_j, ParamP, ls_epsi, ls_roots, ls_demandeN_bis, opt_residu] # input tag

    res_soil_step = loop.step_bilanWN_sol(*tag_inputs_soil_step)
    S, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol = res_soil_step  # unpacks results from a list and updates global variables

    ##########
    # setp update plant stress variables
    ##########

    tag_inputs_stress = [ParamP, invar, invar_sc, temps, DOY, nbplantes, surfsolref, ls_epsi, ls_ftsw, ls_transp,
                            ls_Act_Nuptake_plt, ls_demandeN_bis, ls_ftswStress, ls_TStress, dicOrgans, dicFeuilBilanR, lsApex,
                            start_time, cutNB, deltaI_I0, nbI_I0, I_I0profilLfPlant, I_I0profilPetPlant,
                            I_I0profilInPlant, NlClasses, NaClasses, NlinClasses, outvar]

    invar, invar_sc, outvar, I_I0profilInPlant, ls_ftswStress, ls_NNIStress, ls_TStress = loop.Update_stress_loop(*tag_inputs_stress)

    ##########
    # step update soil residues senescence
    ##########
    tag_inputs_residue_updt = [ls_mat_res, vCC, S, ls_roots, par_SN['PROFHUMs'], ParamP, invar, opt_residu,opt_stressGel] # input tag

    res_residue_step = loop.update_residue_mat(*tag_inputs_residue_updt)
    ls_mat_res, S = res_residue_step  # unpacks results from a list and updates global variables

    #########
    # reinjecte les sorties midiee dans le lsystem
    #########
    lsystem.invar = invar
    lsystem.outvar = outvar
    lsystem.invar_sc = invar_sc

    lsystem.S = S
    lsystem.stateEV = stateEV

    lsystem.res_trans = res_trans
    lsystem.res_abs_i = res_abs_i
    lsystem.res_rfr = res_rfr

    lsystem.ls_ftswStress = ls_ftswStress
    lsystem.ls_NNIStress = ls_NNIStress
    lsystem.ls_TStress = ls_TStress
    lsystem.I_I0profilInPlant = I_I0profilInPlant
