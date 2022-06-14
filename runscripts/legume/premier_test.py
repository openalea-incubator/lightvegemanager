from src.LightVegeManager import *
from src. import *

## --> Installer la dernière version de lpy (3.9.2) pour la cohabitation cnwheat et legume

### ETAPES
# chemin des inputs
# ouvre le excel usm et l'enregistre dans un dico
# fait une liste de lsystem pour chaque situation de l'usm to_run==1 
# 

def runlsystem_bystep(n):
    """run du niem l-system by step dans une liste (names)"""

    lsys = testsim[names[n]]
    lstring = lsys.axiom
    nb_iter = lsys.derivationLength #lire dans derivation_length #335 #30
    lsys.opt_external_coupling = 1 # met a un l'option external coupling

    for i in range(nb_iter+1):
        print('iter ',i,n)
        lstring = lsys.derive(lstring, i, 1)

        ## daily loop
        tag_loop_inputs = lsys.tag_loop_inputs
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

        ############
        # step light transfer coupling
        ############

        # PAR / Blue voxel
        tag_light_inputs = [m_lais / surf_refVOX, triplets, ls_dif, meteo_j['I0'] * surf_refVOX]  # input tag

        # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
        local_res_trans, local_res_abs_i = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, optsky=station['optsky'], opt=station['sky'])

        res_trans, res_abs_i = local_res_trans, local_res_abs_i  # mise a jour variables globales

        # R_FR voxel (calcul de zeta)
        tag_light_inputs2 = [res_trans / (meteo_j['I0'] * surf_refVOX)]  # input tag
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
        lsys.invar = invar
        lsys.outvar = outvar
        lsys.invar_sc = invar_sc

        lsys.S = S
        lsys.stateEV = stateEV

        lsys.res_trans = res_trans
        lsys.res_abs_i = res_abs_i
        lsys.res_rfr = res_rfr

        lsys.ls_ftswStress = ls_ftswStress
        lsys.ls_NNIStress = ls_NNIStress
        lsys.ls_TStress = ls_TStress
        lsys.I_I0profilInPlant = I_I0profilInPlant

    testsim[names[n]].clear()
    print((''.join((names[n], " - done"))))

if __name__ == '__main__':
    runlsystem_bystep(0)