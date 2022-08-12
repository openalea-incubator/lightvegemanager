import os
import shutil
import sys
import time
import progressbar
import getopt
import random

try :
    from src.LightVegeManager import *
    from src.FSPMWheat_template import *
    from src.l_egume_template import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
    from src.LightVegeManager import *
    from src.FSPMWheat_template import *
    from src.l_egume_template import *

random.seed(1234)
np.random.seed(1234)

def Create_Folders(parentfolderpath, newfolder):
    # dossier des données brutes
    dirName = os.path.join(parentfolderpath, newfolder)
   
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")

def main_init_cnwheat(SENESCWHEAT_TIMESTEP, 
                        ELONGWHEAT_TIMESTEP, 
                        GROWTHWHEAT_TIMESTEP, 
                        CNWHEAT_TIMESTEP, 
                        outfolderpath, 
                        inputs=""):
    
    outfolderpath = os.path.join(outfolderpath, "cnwheat")
    
    PLANT_DENSITY = {1: 250.}
    tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5}
    N_fertilizations={2016: 357143, 2520: 1000000}

    # supprime le contenant du dossier outputs si non vide
    if os.path.exists(outfolderpath):
        for filename in os.listdir(outfolderpath):
            file_path = os.path.join(outfolderpath, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

    # dossier de sortie
    Create_Folders(outfolderpath, "brut")
    Create_Folders(outfolderpath, "postprocessing")

    # Path of the directory which contains the inputs of the model
    INPUTS_FOLDER = inputs
    
    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    METEO_FILENAME = 'meteo_Ljutovac2002.csv'
    PHYTOT_FILENAME = 'phytoT.csv'

    # Path of the directory which contains the outputs of the model
    OUTPUTS_DIRPATH = outfolderpath+"/brut/"

    # Name of the CSV files which will contain the outputs of the model
    AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
    ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
    HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
    ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
    SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

    # Path of the directory which contains the postprocessing
    POSTPROCESSING_DIRPATH = outfolderpath+"/postprocessing/"
    
    # Name of the CSV files which will contain the postprocessing of the model
    AXES_POSTPROCESSING_FILENAME = 'axes_postprocessing.csv'
    ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
    HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
    ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
    SOILS_POSTPROCESSING_FILENAME = 'soils_postprocessing.csv'

    # parameters
    LEAVES_MODEL = 'Soissons_byleafclass'
    update_parameters_all_models=None

    meteo, \
        g, \
        adel_wheat, \
        elongwheat_facade_, growthwheat_facade_, \
        farquharwheat_facade_, senescwheat_facade_, \
        cnwheat_facade_, fspmwheat_facade_ = initialisation_fspmwheat(ELONGWHEAT_TIMESTEP, 
                                                                        GROWTHWHEAT_TIMESTEP, 
                                                                        SENESCWHEAT_TIMESTEP, 
                                                                        CNWHEAT_TIMESTEP,
                                                                        INPUTS_FOLDER,
                                                                        METEO_FILENAME,
                                                                        AXES_INITIAL_STATE_FILENAME,
                                                                        ORGANS_INITIAL_STATE_FILENAME,
                                                                        HIDDENZONES_INITIAL_STATE_FILENAME,
                                                                        ELEMENTS_INITIAL_STATE_FILENAME,
                                                                        SOILS_INITIAL_STATE_FILENAME,
                                                                        PHYTOT_FILENAME,
                                                                        PLANT_DENSITY,
                                                                        LEAVES_MODEL,
                                                                        update_parameters_all_models,
                                                                        N_fertilizations=N_fertilizations)
    return PLANT_DENSITY, tillers_replications, N_fertilizations, \
            meteo, \
            g, \
            adel_wheat, \
            elongwheat_facade_, growthwheat_facade_, \
            farquharwheat_facade_, senescwheat_facade_, \
            cnwheat_facade_, fspmwheat_facade_ ,\
            OUTPUTS_DIRPATH,AXES_OUTPUTS_FILENAME, \
            ORGANS_OUTPUTS_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, \
            ELEMENTS_OUTPUTS_FILENAME, SOILS_OUTPUTS_FILENAME, \
            POSTPROCESSING_DIRPATH, AXES_POSTPROCESSING_FILENAME, \
            ORGANS_POSTPROCESSING_FILENAME, ORGANS_POSTPROCESSING_FILENAME, \
            HIDDENZONES_POSTPROCESSING_FILENAME, ELEMENTS_POSTPROCESSING_FILENAME, SOILS_POSTPROCESSING_FILENAME

def main_init_legume(foldin, foldout):
    foldout = os.path.join(foldout, "legume")

    fusms = "liste_usms_exemple.xls"
    ongletBatch = "exemple"

    # prépare une liste de l-system à lancer
    lsystem_simulations, names_simulations = initialisation_legume(foldin, foldout, fusms, ongletBatch)

    index_simulation = 0
    sim_id = names_simulations[index_simulation]
    lstring = lsystem_simulations[sim_id].axiom
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
    
    return lsystem_simulations, sim_id, lstring, dxyz, lsys_temp.cote

def run_light_legume(current_day, next_day_from_next_hour, meteo):
    '''Active le calcul de lumière pour l-egume

    à chaque fin de journée (condition si l'heure d'après est un nouveau jour)
    on fait le rayonnement global de la journée
    
    '''
    # on est toujours au même jour
    if current_day == next_day_from_next_hour :
        return 0, False
    
    else:
        filterday = meteo[meteo["DOY"] == current_day]
        mean_pari = filterday["PARi"].sum(axis=0)/len(filterday)

        # micromol.m-2.s-2 to J.cm²
        # conversion en global W.m-2
        mean_pari = mean_pari/2.02

        # W.m-2 = J.s-1.m-2
        # conversion en J.s-1.m-2
        mean_pari = mean_pari/10000

        # conversion en W.s-1.m-2 (8.64 tiré d'un fichier météo de l-egume)
        RG = mean_pari / 8.64

        return RG, True
    

def simulation(level_tesselation, SIMULATION_LENGTH, legumeinputs="", cnwheatinputs="", outfolderpath="", distrialgo="global", writegeo=False):
    # On va suivre la météo de cn-wheat
    # define the time step in hours for each simulator
    LIGHT_TIMESTEP = 4
    SENESCWHEAT_TIMESTEP = 1
    FARQUHARWHEAT_TIMESTEP = 1
    ELONGWHEAT_TIMESTEP = 1
    GROWTHWHEAT_TIMESTEP = 1
    CNWHEAT_TIMESTEP = 1

    # initialisation CN-Wheat
    PLANT_DENSITY, tillers_replications, N_fertilizations, \
    meteo, \
    g, \
    adel_wheat, \
    elongwheat_facade_, growthwheat_facade_, \
    farquharwheat_facade_, senescwheat_facade_, \
    cnwheat_facade_, fspmwheat_facade_ ,\
    OUTPUTS_DIRPATH,AXES_OUTPUTS_FILENAME, \
    ORGANS_OUTPUTS_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, \
    ELEMENTS_OUTPUTS_FILENAME, SOILS_OUTPUTS_FILENAME, \
    POSTPROCESSING_DIRPATH, AXES_POSTPROCESSING_FILENAME, \
    ORGANS_POSTPROCESSING_FILENAME, ORGANS_POSTPROCESSING_FILENAME, \
    HIDDENZONES_POSTPROCESSING_FILENAME, ELEMENTS_POSTPROCESSING_FILENAME,\
    SOILS_POSTPROCESSING_FILENAME = main_init_cnwheat(SENESCWHEAT_TIMESTEP, 
                                                        ELONGWHEAT_TIMESTEP, 
                                                        GROWTHWHEAT_TIMESTEP, 
                                                        CNWHEAT_TIMESTEP, outfolderpath, cnwheatinputs)
    
    # initialisation l-egume
    lsystem_simulations, sim_id, lstring,  dxyz, cote = main_init_legume(legumeinputs, outfolderpath)
    
    # -- SIMULATION PARAMETERS --
    START_TIME = 0    

    # paramètre communs
    environment = {}
    environment["coordinates"] = [48.85,0,0] # latitude, longitude, timezone
    environment["sky"] = "turtle46" # turtle à 46 directions par défaut
    environment["diffus"] = True
    environment["reflected"] = False
    environment["infinite"] = False
    environment["reflectance coefficients"] = [[0.1, 0.05], [0., 0.]]

    ratp_parameters = {}
    
    ratp_parameters["soil reflectance"] = [0., 0.]
    ratp_parameters["mu"] = [1.]
    ratp_parameters["tesselation level"] = level_tesselation
    if distrialgo=="voxel" : ratp_parameters["angle distrib algo"] = "compute voxel"
    elif distrialgo=="global" : ratp_parameters["angle distrib algo"] = "compute global"
    ratp_parameters["nb angle classes"] = 30

    ## Paramètres RATP ##
    # CN-Wheat : direct et diffus, taille de voxel dynamique
    environment["direct"] = True 
    ratp_parameters["voxel size"] = "dynamic"
    lghtratpcnwheat = LightVegeManager(environment=environment,
                                    lightmodel="ratp",
                                    lightmodel_parameters=ratp_parameters,
                                    main_unit="m")

    # l-egume : diffus, taille de voxel de la grille RIRI de l-egume
    environment["direct"] = False
    ratp_parameters["voxel size"] = [d*0.01 for d in dxyz] # on convertit en m
    ratp_parameters["xy max"] = cote
    lghtratplegume = LightVegeManager(environment=environment,
                                    lightmodel="ratp",
                                    lightmodel_parameters=ratp_parameters,
                                    main_unit="m")


    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------
        
    tot_light = 0.
    i_legume = 0

    # CN-wheat tables
    axes_all_data_list = []
    organs_all_data_list = []  # organs which belong to axes: roots, phloem, grains
    hiddenzones_all_data_list = []
    elements_all_data_list = []
    soils_all_data_list = []
    all_simulation_steps = []  # to store the steps of the simulation

    #for t_light in progressbar.progressbar(range(START_TIME, SIMULATION_LENGTH, LIGHT_TIMESTEP)):
    for t_light in progressbar.progressbar(range(START_TIME, SIMULATION_LENGTH, SENESCWHEAT_TIMESTEP)):
        print("--- ",t_light)

        light_start=time.time()
        # récupère les données météo
        PARi = meteo.loc[t_light, ['PARi']].iloc[0]
        DOY = meteo.loc[t_light, ['DOY']].iloc[0]
        hour = meteo.loc[t_light, ['hour']].iloc[0]
        PARi_next_hours = meteo.loc[range(t_light, t_light + LIGHT_TIMESTEP), ['PARi']].sum().values[0]
        next_day_next_hour = meteo.loc[t_light + LIGHT_TIMESTEP, ['DOY']].iloc[0]

        # vérifie si l'itération suivante est encore le jour? et lance le calcul de lumière
        light_run = light_run = (t_light % LIGHT_TIMESTEP == 0) and (PARi_next_hours > 0)

        if light_run:
            geometry = {}

            # scene cn-wheat : PlantGL
            cnwheat_scene = adel_wheat.scene(g)

            # scene l-egume : PlantGL
            legume_scene = lsystem_simulations[sim_id].sceneInterpretation(lstring)

            # copie de la scène
            geometry["scenes"] = [cnwheat_scene, legume_scene]
            geometry["transformations"] = {}
            geometry["transformations"]["scenes unit"] = ["m", "cm"]
            
            lghtratpcnwheat.init_scenes(geometry)
            lghtratpcnwheat.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
            lghtratpcnwheat.PAR_update_MTG(g)
        
        # sinon on copie Erel de l'itération précédente
        else:
            Erel = g.property('Erel')
            PARa_output = {k: v * PARi for k, v in Erel.items()}
            outputs = {}
            outputs.update({'PARa': PARa_output})
            for param in outputs.keys():
                if param not in g.properties():
                    g.add_property(param)
                # update the MTG
                g.property(param).update(outputs[param])
        
        RG, light_legume = run_light_legume(DOY, next_day_next_hour, meteo)

        if light_legume :
            lstring = lsystem_simulations[sim_id].derive(lstring, i_legume, 1)
        
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

            geometry = {}

            # scene cn-wheat : PlantGL
            cnwheat_scene = adel_wheat.scene(g)

            # scene l-egume : PlantGL
            legume_scene = lsystem_simulations[sim_id].sceneInterpretation(lstring)

            # copie de la scène
            geometry["scenes"] = [cnwheat_scene, legume_scene]
            
            lghtratplegume.init_scenes(geometry)
            lghtratplegume.run(PARi=RG, day=DOY, hour=hour, truesolartime=True, parunit="RG") 
            res_trans, res_abs_i = lghtratplegume.to_l_egume(m_lais, RG)

            iteration_legume_withoutlighting(lsystem_simulations[sim_id], res_trans, res_abs_i, tag_loop_inputs)
            i_legume += 1
        
        tot_light += (time.time() - light_start)

        iteration_fspmwheat_withoutlighting(meteo,
                                            g,
                                            adel_wheat,
                                            elongwheat_facade_,
                                            growthwheat_facade_,
                                            farquharwheat_facade_,
                                            senescwheat_facade_,
                                            cnwheat_facade_,
                                            fspmwheat_facade_,
                                            ELONGWHEAT_TIMESTEP, 
                                            GROWTHWHEAT_TIMESTEP, 
                                            SENESCWHEAT_TIMESTEP, 
                                            CNWHEAT_TIMESTEP,
                                            FARQUHARWHEAT_TIMESTEP,
                                            LIGHT_TIMESTEP,
                                            t_light,
                                            save_df=True,
                                            axes_all_data_list=axes_all_data_list,
                                            organs_all_data_list=organs_all_data_list,
                                            hiddenzones_all_data_list=hiddenzones_all_data_list,
                                            elements_all_data_list=elements_all_data_list,
                                            soils_all_data_list=soils_all_data_list,
                                            all_simulation_steps=all_simulation_steps,
                                            N_fertilizations=N_fertilizations,
                                            tillers_replications=tillers_replications)

        if writegeo :
            lghtratpcnwheat.VTKout(os.path.join(outfolderpath, "vtk/smallcouplage")+"_"+str(t_light)+"_", voxels=True)
    
    write_outputs_fspmwheat(OUTPUTS_DIRPATH,
                                POSTPROCESSING_DIRPATH,
                                AXES_OUTPUTS_FILENAME,
                                ORGANS_OUTPUTS_FILENAME,
                                HIDDENZONES_OUTPUTS_FILENAME,
                                ELEMENTS_OUTPUTS_FILENAME,
                                SOILS_OUTPUTS_FILENAME,
                                AXES_POSTPROCESSING_FILENAME,
                                ORGANS_POSTPROCESSING_FILENAME,
                                HIDDENZONES_POSTPROCESSING_FILENAME,
                                ELEMENTS_POSTPROCESSING_FILENAME,
                                SOILS_POSTPROCESSING_FILENAME,
                                axes_all_data_list,
                                organs_all_data_list,
                                hiddenzones_all_data_list,
                                elements_all_data_list,
                                soils_all_data_list,
                                all_simulation_steps)

if __name__ == "__main__":
    # valeur par défaut
    legumeinputs = os.path.normpath("C:/Users/mwoussen/cdd/codes/vegecouplelight/l-egume/legume/input/")
    cnwheatinputs = os.path.normpath("C:/Users/mwoussen/cdd/codes/vegecouplelight/WheatFspm/fspm-wheat/example/Vegetative_stages/inputs")
    outfolderpath = os.path.normpath("C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/couplage_smalltests/")
    writegeo = True
    SIMULATION_LENGTH = 60
    level_tesselation = 1
    distrialgo="global"

    # dossiers de sortie
    Create_Folders(outfolderpath, "cnwheat")
    Create_Folders(outfolderpath, "legume")
    Create_Folders(outfolderpath, "vtk")

    print('--- BEGIN ---')
    simulation(level_tesselation, SIMULATION_LENGTH, legumeinputs, cnwheatinputs, outfolderpath, distrialgo, writegeo)
    print('--- END ---')

        

    