from src.LightVegeManager import *
from src.l_egume_template import *

## --> Installer la dernière version de lpy (3.9.2) pour la cohabitation cnwheat et legume

### ETAPES
# chemin des inputs
# ouvre le excel usm et l'enregistre dans un dico
# fait une liste de lsystem pour chaque situation de l'usm to_run==1 

foldin = "C:/Users/mwoussen/cdd/codes/vegecouplelight/l-egume/legume/input/"
foldout = "C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume/"
fusms = "liste_usms_exemple.xls"
ongletBatch = "exemple"

# prépare une liste de l-system à lancer
lsystem_simulations, names_simulations = initialisation_legume(foldin, foldout, fusms, ongletBatch)

index_simulation = 0
sim_id = names_simulations[index_simulation]
lstring = lsystem_simulations[sim_id].axiom
nb_iter = lsystem_simulations[sim_id].derivationLength #lire dans derivation_length #335 #30
lsystem_simulations[sim_id].opt_external_coupling = 1 # met a un l'option external coupling

## INITIALISATION LIGHTVEGEMANAGER
geometry = {}
environment = {}
ratp_parameters = {}
caribu_parameters = {}

# Paramètres pré-simulation
environment["names"] = ["l-egume"]
environment["coordinates"] = [46.43,0,0] # latitude, longitude, timezone
environment["sky"] = "turtle46" # turtle à 46 directions par défaut
environment["diffus"] = True
environment["direct"] = False
environment["reflected"] = False
environment["reflectance coefficients"] = [[0.1, 0.05]]
environment["infinite"] = True

## PARAMETRES RATP
dx, dy, dz = 4, 4, 2 # m
ratp_parameters["voxel size"] = [dx, dy, dz]
ratp_parameters["soil reflectance"] = [0., 0.]
ratp_parameters["mu"] = [1.]

lghtratp = LightVegeManager(environment=environment,
                                lightmodel="ratp",
                                lightmodel_parameters=ratp_parameters)

# début de la simulation
for i in range(nb_iter+1):
    print('time step: ',i)
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
        opt_stressW, opt_stressN, opt_stressGel, opt_residu = tag_loop_inputs

    
    # preparation lightvegemanager
    scene_legume = [m_lais / surf_refVOX, ls_dif, meteo_j['I0'] * surf_refVOX]
    geometry["scenes"] = [scene_legume]
    lghtratp.init_scenes(geometry)  
    
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
    lsystem_simulations[sim_id].invar = invar
    lsystem_simulations[sim_id].outvar = outvar
    lsystem_simulations[sim_id].invar_sc = invar_sc

    lsystem_simulations[sim_id].S = S
    lsystem_simulations[sim_id].stateEV = stateEV

    lsystem_simulations[sim_id].res_trans = res_trans
    lsystem_simulations[sim_id].res_abs_i = res_abs_i
    lsystem_simulations[sim_id].res_rfr = res_rfr

    lsystem_simulations[sim_id].ls_ftswStress = ls_ftswStress
    lsystem_simulations[sim_id].ls_NNIStress = ls_NNIStress
    lsystem_simulations[sim_id].ls_TStress = ls_TStress
    lsystem_simulations[sim_id].I_I0profilInPlant = I_I0profilInPlant

lsystem_simulations[sim_id].clear()
print((''.join((sim_id, " - done"))))