import sys
import os
import getopt

import numpy as np

try :
    from src.LightVegeManager import *
    from src.l_egume_template import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
    from src.LightVegeManager import *
    from src.l_egume_template import *

'''
Comparaison sur un temps avec un couvert dense entre CARIBU et RATP
Simu 1 : RATP seul
Simu 2 : l-egume seul
Simu 3 : l-egume + RATP sur la géométrie de l-egume

NOTES :
    vérifier le chemin de INPUTS_FOLDER et celui de runstring pour lancer le script par défaut

'''

def simulation(foldin, foldout, active, passive, ratpgeo, writegeo=False):
    
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
        opt_stressW, opt_stressN, opt_stressGel, opt_residu = tag_loop_inputs
    dxyz = lsys_temp.dxyz # récupère nouvelle taille de voxel
    nxyz = lsys_temp.na # récupère le nombre de voxels

    if active=="ratp" or passive=="ratp" :
        ## INITIALISATION LIGHTVEGEMANAGER
        environment = {}
        ratp_parameters_legume = {}
        ratp_parameters_plantgl = {}

        # Paramètres pré-simulation
        environment["names"] = ["l-egume"]
        environment["coordinates"] = [46.43,0,0] # latitude, longitude, timezone
        environment["sky"] ="turtle46" # ["file", "runscripts/legume/sky_5.data"] # "turtle46" # turtle à 46 directions par défaut
        environment["diffus"] = True
        environment["direct"] = False
        environment["reflected"] = False
        environment["reflectance coefficients"] = [[0., 0.]]
        environment["infinite"] = False

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

        if active=="ratp" or passive=="ratp" :
            ## Paramètres météo ## 
            doy = lsystem_simulations[sim_id].meteo["DOY"][i]
            hour = 12
            energy = meteo_j['I0']     

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
                elif ratpgeo == "plantgl":
                    lghtratp.VTKinit(foldout+"plantgl_"+str(i)+"_")
        
        ############
        # step light transfer coupling
        ############
        if active=="legume":
            # PAR / Blue voxel
            tag_light_inputs = [m_lais / surf_refVOX, triplets, ls_dif, meteo_j['I0'] * surf_refVOX]  # input tag

            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, optsky=station['optsky'], opt=station['sky'])

        elif active=="ratp":
            # transfert des sorties
            res_abs_i = np.zeros((m_lais.shape[0], nxyz[2], nxyz[1], nxyz[0]))
            res_trans = np.zeros((nxyz[2], nxyz[1], nxyz[0]))
            for iz in range(m_lais.shape[1]):
                for iy in range(m_lais.shape[2]):
                    for ix in range(m_lais.shape[3]):
                        if m_lais[0][iz][iy][ix] > 0. :
                            vox_data = lghtratp.voxels_outputs[(lghtratp.voxels_outputs.Nx==ix+1) & 
                                                                        (lghtratp.voxels_outputs.Ny==iy+1) & 
                                                                        (lghtratp.voxels_outputs.Nz==iz+1)]
                        
                            res_trans[iz, iy, ix] = energy * sum(vox_data["transmitted"])
                            for ie in range(m_lais.shape[0]) :
                                if len(vox_data) > 0 :
                                    res_abs_i[ie, iz, iy, ix] = energy * vox_data[vox_data.VegetationType == ie+1]["xintav"].values[0]

        iteration_legume_withoutlighting(lsystem_simulations[sim_id], res_trans, res_abs_i, tag_loop_inputs)

    print((''.join((sim_id, " - done"))))

    lsystem_simulations[sim_id].clear()

if __name__ == "__main__":
    #definition d'arguments avec getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i: o: a: p: r: w:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    # valeur par défaut
    foldin = "C:/Users/mwoussen/cdd/codes/vegecouplelight/l-egume/legume/input/"
    foldout = "C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume/"
    active = "legume"
    passive = "legume"
    ratpgeo = "grid"
    writegeo = "no"

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
        simulation(foldin, foldout, active, passive, ratpgeo, True)
    simulation(foldin, foldout, active, passive, ratpgeo)
    
    print("=== END ===")