from src.LightVegeManager import *
from src.FSPMWheat_facade import *

import os
import shutil
import sys
import time
import progressbar
import getopt

'''
Comparaison sur un temps avec un couvert dense entre CARIBU et RATP
Simu 1 : CARIBU + RATP sur la géo de CARIBU
Simu 2 : RATP seul

Compare :
    * PAR sortie (inutile de retenir les nuits)
    * temps de calcul
    * sorties CN-Wheat
'''

def Create_OutputsFolders(parentfolderpath):
    # dossier des données brutes
    dirName = parentfolderpath+"/brut"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
    
    # dossier du postprocessing
    dirName = parentfolderpath+"/postprocessing"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
    

def simulation_caribu(level_tesselation, inter_row, SIMULATION_LENGTH, outfolderpath):
    # -- SIMULATION PARAMETERS --
    START_TIME = 0
    PLANT_DENSITY = {1: 410}

    # define the time step in hours for each simulator
    LIGHT_TIMESTEP = 4
    SENESCWHEAT_TIMESTEP = 2
    FARQUHARWHEAT_TIMESTEP = 2
    ELONGWHEAT_TIMESTEP = 1
    GROWTHWHEAT_TIMESTEP = 1
    CNWHEAT_TIMESTEP = 1

    # fspm-wheat file name
    current_path = os.path.dirname(os.path.abspath(__file__))
    lvm_folder = "/".join(current_path.split("/")[:-1])
    if lvm_folder == "": lvm_folder = "/".join(current_path.split("\\")[:-1])
    Create_OutputsFolders(outfolderpath+"/caribu")
    
    # Path of the directory which contains the inputs of the model
    INPUTS_FOLDER = lvm_folder+'/WheatFspm/fspm-wheat/test/inputs'
    
    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    METEO_FILENAME = 'meteo_Ljutovac2002.csv'
    PHYTOT_FILENAME = 'phytoT.csv'

    # Path of the directory which contains the outputs of the model
    OUTPUTS_DIRPATH = outfolderpath+"/caribu/brut/"

    # Name of the CSV files which will contain the outputs of the model
    AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
    ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
    HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
    ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
    SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

    # Path of the directory which contains the postprocessing
    POSTPROCESSING_DIRPATH = outfolderpath+"/caribu/postprocessing/"
    
    # Name of the CSV files which will contain the postprocessing of the model
    AXES_POSTPROCESSING_FILENAME = 'axes_postprocessing.csv'
    ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
    HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
    ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
    SOILS_POSTPROCESSING_FILENAME = 'soils_postprocessing.csv'

    # parameters
    LEAVES_MODEL = 'Soissons_byleafclass'
    SL_RATIO_D = 0.25
    AGE_EFFECT_SENESCENCE = 10000
    VMAX_ROOTS_GROWTH_PREFLO = 0.02885625
    K_AMINO_ACIDS_EXPORT = 3E-5
    K_NITRATE_EXPORT = 1E-6

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
                                                                        SL_RATIO_D,
                                                                        AGE_EFFECT_SENESCENCE,
                                                                        VMAX_ROOTS_GROWTH_PREFLO,
                                                                        K_AMINO_ACIDS_EXPORT,
                                                                        K_NITRATE_EXPORT)
    
    # define lists of dataframes to store the inputs and the outputs of the models at each step.
    axes_all_data_list = []
    organs_all_data_list = []  # organs which belong to axes: roots, phloem, grains
    hiddenzones_all_data_list = []
    elements_all_data_list = []
    soils_all_data_list = []

    all_simulation_steps = []  # to store the steps of the simulation
    
    # Paramètres pré-simulation
    model_names = ["fspm-wheat"]
    coordinates = [43.,0,0] # latitude, longitude, timezone
    
    ## Paramètres CARIBU ##
    sun_sky_options = "mix"
    sun_algo="ratp" # soleil avec RATP sundirection dans Shortwave_Balance.f90 
    infinite=True # avec couvert inifini
    caribu_param = [sun_sky_options, sun_algo, infinite, None]
    caribu_sky = [] # ciel turtle à 46 directions
    caribu_rf = [[0.1, 0.05]] # réflectance, transmittance

    ## Paramètres RATP ##
    dv = 0.02 # m
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    ratp_mu = [1.]
    tesselate_level = level_tesselation
    distrib_algo = "compute voxel" # "file"
    distrib_option = 30
    infinite=True
    ratp_parameters = [dx, dy, dz, rs, ratp_mu, tesselate_level, distrib_algo, distrib_option, infinite]
    ratp_sky = [] # ciel turtle à 46 directions
    ratp_rf=[[0., 0]] # leaf reflectance PAR et NIR pour l'entité 0

    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    para_c=[]
    para_r=[]
    pari=[]
    parin=[]
    parsun=[]
    parsha=[]
    iter=[]
    shapes=[]
    caribu_times=[]
    ratp_times=[]
    tot_light = 0.
    for t_light in progressbar.progressbar(range(START_TIME, SIMULATION_LENGTH, LIGHT_TIMESTEP)):
        print("\n")
        light_start=time.time()
        # récupère les données météo
        PARi = meteo.loc[t_light, ['PARi_MA4']].iloc[0]
        DOY = meteo.loc[t_light, ['DOY']].iloc[0]
        hour = meteo.loc[t_light, ['hour']].iloc[0]
        PARi_next_hours = meteo.loc[range(t_light, t_light + LIGHT_TIMESTEP), ['PARi']].sum().values[0]

        # création d'un couvert hétérogène
        scene_etendu, domain = create_heterogeneous_canopy_copy(adel_wheat, g, nplants=100, var_plant_position=0.03, var_leaf_inclination=0.157, var_leaf_azimut=1.57, var_stem_azimut=0.157,
                                     plant_density=250, inter_row=inter_row)

        # vérifie si l'itération suivante est encore le jour? et lance le calcul de lumière
        if (t_light % LIGHT_TIMESTEP == 0) and (PARi_next_hours > 0):
            in_scenes = [scene_etendu]
            # recherche des tiges dans les différentes scènes
            id_entity = 0
            id_stems=whichstems_MTG(g, id_entity)

            # ajoute le pattern aux paramètres du modèle
            caribu_param[-1] = domain
            
            c_time = time.time()
            # Objet calcul de la lumière
            lghtcaribu = LightVegeManager(in_scenes, 
                                        in_names=model_names,
                                        id_stems=id_stems, 
                                        sky_parameters=caribu_sky,
                                        lightmodel="caribu", lightmodelparam=caribu_param, 
                                        rf=caribu_rf, 
                                        coordinates=coordinates)
            lghtcaribu.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
            c_time = time.time()-c_time
            
            lghtcaribu.PAR_update_MTG(g)
            
            r_time = time.time()
            lghtratp = LightVegeManager(in_scenes, 
                                        in_names=model_names,
                                        id_stems=id_stems, 
                                        sky_parameters=ratp_sky,
                                        lightmodel="ratp", lightmodelparam=ratp_parameters, 
                                        rf=ratp_rf, 
                                        coordinates=coordinates)
            lghtratp.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
            r_time = time.time()-r_time

            for key,items in g.property("PARa").items():
                para_c.append(items)
                para_r.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["PARa"].values[0])
                shapes.append(key)
                pari.append(lghtcaribu.shapes_outputs[lghtcaribu.shapes_outputs.ShapeId==key]["PARi"].values[0])
                parsun.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["SunlitPAR"].values[0])
                parsha.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["ShadedPAR"].values[0])
            iter.extend([t_light]*len(g.property("PARa")))
            parin.extend([PARi]*len(g.property("PARa")))
            caribu_times.extend([c_time]*len(g.property("PARa")))
            ratp_times.extend([r_time]*len(g.property("PARa")))
        
        # sinon on copie le PAR de l'itération précédente
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
                                            all_simulation_steps=all_simulation_steps)

    myoutputs = {"Shapes" : shapes, 
                    "Iteration" : iter,
                    "CARIBU time" : caribu_times,
                    "RATP time" : ratp_times,
                    "PAR input":parin, 
                    "PARa CARIBU" : para_c, 
                    "PARi CARIBU" : pari, 
                    "PARa RATP" : para_r, 
                    "SunlitPAR" : parsun,
                    "ShadedPAR" : parsha}
    df_myout = pd.DataFrame(myoutputs)
    df_myout.to_csv(outfolderpath+"/PAR_values_caribu_geom.csv")

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
    print("--- temps execution : ",tot_light)

def simulation_ratp(level_tesselation, inter_row, SIMULATION_LENGTH, outfolderpath):
    # -- SIMULATION PARAMETERS --
    START_TIME = 0
    PLANT_DENSITY = {1: 410}

    # define the time step in hours for each simulator
    LIGHT_TIMESTEP = 4
    SENESCWHEAT_TIMESTEP = 2
    FARQUHARWHEAT_TIMESTEP = 2
    ELONGWHEAT_TIMESTEP = 1
    GROWTHWHEAT_TIMESTEP = 1
    CNWHEAT_TIMESTEP = 1

    # fspm-wheat file name
    current_path = os.path.dirname(os.path.abspath(__file__))
    lvm_folder = "/".join(current_path.split("/")[:-1])
    if lvm_folder == "": lvm_folder = "/".join(current_path.split("\\")[:-1])
    Create_OutputsFolders(outfolderpath+"/ratp")
    
    # Path of the directory which contains the inputs of the model
    INPUTS_FOLDER = lvm_folder+'/WheatFspm/fspm-wheat/test/inputs'

    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    METEO_FILENAME = 'meteo_Ljutovac2002.csv'
    PHYTOT_FILENAME = 'phytoT.csv'

    # Path of the directory which contains the outputs of the model
    OUTPUTS_DIRPATH = outfolderpath+"/ratp/brut/"

    # Name of the CSV files which will contain the outputs of the model
    AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
    ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
    HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
    ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
    SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

    # Path of the directory which contains the postprocessing
    POSTPROCESSING_DIRPATH = outfolderpath+"/ratp/postprocessing/"
    
    # Name of the CSV files which will contain the postprocessing of the model
    AXES_POSTPROCESSING_FILENAME = 'axes_postprocessing.csv'
    ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
    HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
    ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
    SOILS_POSTPROCESSING_FILENAME = 'soils_postprocessing.csv'

    # parameters
    LEAVES_MODEL = 'Soissons_byleafclass'
    SL_RATIO_D = 0.25
    AGE_EFFECT_SENESCENCE = 10000
    VMAX_ROOTS_GROWTH_PREFLO = 0.02885625
    K_AMINO_ACIDS_EXPORT = 3E-5
    K_NITRATE_EXPORT = 1E-6

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
                                                                        SL_RATIO_D,
                                                                        AGE_EFFECT_SENESCENCE,
                                                                        VMAX_ROOTS_GROWTH_PREFLO,
                                                                        K_AMINO_ACIDS_EXPORT,
                                                                        K_NITRATE_EXPORT)

    # define lists of dataframes to store the inputs and the outputs of the models at each step.
    
    # Paramètres pré-simulation
    model_names = ["fspm-wheat"]
    coordinates = [43.,0,0] # latitude, longitude, timezone

    ## Paramètres RATP ##
    dv = 0.02 # m
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    ratp_mu = [1.]
    tesselate_level = level_tesselation
    distrib_algo = "compute voxel" # "file"
    distrib_option = 30
    infinite=True
    ratp_parameters = [dx, dy, dz, rs, ratp_mu, tesselate_level, distrib_algo, distrib_option, infinite]
    ratp_sky = [] # ciel turtle à 46 directions
    ratp_rf=[[0., 0]] # leaf reflectance PAR et NIR pour l'entité 0

    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    
    tot_light = 0.
    for t_light in progressbar.progressbar(range(START_TIME, SIMULATION_LENGTH, LIGHT_TIMESTEP)):
        axes_all_data_list = []
        organs_all_data_list = []  # organs which belong to axes: roots, phloem, grains
        hiddenzones_all_data_list = []
        elements_all_data_list = []
        soils_all_data_list = []
        all_simulation_steps = []  # to store the steps of the simulation

        para_r=[]
        parin=[]
        parsun=[]
        parsha=[]
        iter=[]
        shapes=[]
        ratp_times=[]

        print("\n")
        light_start=time.time()
        # récupère les données météo
        PARi = meteo.loc[t_light, ['PARi_MA4']].iloc[0]
        DOY = meteo.loc[t_light, ['DOY']].iloc[0]
        hour = meteo.loc[t_light, ['hour']].iloc[0]
        PARi_next_hours = meteo.loc[range(t_light, t_light + LIGHT_TIMESTEP), ['PARi']].sum().values[0]

        # création d'un couvert hétérogène
        scene_etendu, domain = create_heterogeneous_canopy_copy(adel_wheat, g, nplants=100, var_plant_position=0.03, var_leaf_inclination=0.157, var_leaf_azimut=1.57, var_stem_azimut=0.157,
                                     plant_density=250, inter_row=inter_row)

        # vérifie si l'itération suivante est encore le jour? et lance le calcul de lumière
        if (t_light % LIGHT_TIMESTEP == 0) and (PARi_next_hours > 0):
            in_scenes = [scene_etendu]
            # recherche des tiges dans les différentes scènes
            id_entity = 0
            id_stems=whichstems_MTG(g, id_entity)            

            r_time = time.time()
            lghtratp = LightVegeManager(in_scenes, 
                                        in_names=model_names,
                                        id_stems=id_stems, 
                                        sky_parameters=ratp_sky,
                                        lightmodel="ratp", lightmodelparam=ratp_parameters, 
                                        rf=ratp_rf, 
                                        coordinates=coordinates)
            lghtratp.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
            r_time = time.time()-r_time
            lghtratp.PAR_update_MTG(g)       

            # enregistre les valeurs de PAR
            for key,items in g.property("PARa").items():
                para_r.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["PARa"].values[0])
                shapes.append(key)
                parsun.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["SunlitPAR"].values[0])
                parsha.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["ShadedPAR"].values[0])
            iter.extend([t_light]*len(g.property("PARa")))
            parin.extend([PARi]*len(g.property("PARa")))
            ratp_times.extend([r_time]*len(g.property("PARa")))  

            myoutputs = {"Shapes" : shapes, 
                        "Iteration" : iter,
                        "RATP time" : ratp_times, 
                        "PAR input":parin, 
                        "PARa RATP" : para_r, 
                        "SunlitPAR" : parsun,
                        "ShadedPAR" : parsha}
            df_myout = pd.DataFrame(myoutputs)
            if not os.path.exists(outfolderpath+"/PAR_values_ratp_geom.csv"):
                df_myout.to_csv(outfolderpath+"/PAR_values_ratp_geom.csv", mode='w')
            else:
                df_myout.to_csv(outfolderpath+"/PAR_values_ratp_geom.csv", mode='a', header=False)   
        
        # sinon on copie le PAR de l'itération précédente
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
                                            all_simulation_steps=all_simulation_steps)

        append_outputs_fspmwheat(OUTPUTS_DIRPATH,
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
    
    print("--- temps execution : ",tot_light)

if __name__ == "__main__":
    #definition d'arguments avec getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s: l: i: n: o:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    # valeur par défaut
    level_tesselation=2
    inter_row=0.1
    nstep=16
    outfolderpath = "outputs/cn-wheat_dense"
    sim = 3

    # récupère les arguments en entrée
    for opt, arg in opts:
        if opt in ("-l"):
            level_tesselation = int(arg)
        elif opt in ("-o"):
            outfolderpath = str(arg)
        elif opt in ("-n"):
            nstep = int(arg)
        elif opt in ("-i"):
            inter_row = int(arg)
        elif opt in ("-s"):
            sim = int(arg)

    # supprime le contenant du dossier outputs si non vide
    for filename in os.listdir(outfolderpath):
        file_path = os.path.join(outfolderpath, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))
    
    # construit les dossiers s'ils n'existent pas
    dirName = outfolderpath+"/caribu"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
    dirName = outfolderpath+"/ratp"
    try:
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
    
    print("=== BEGIN ===")
    if sim==1:
        print("=== === CARIBU GEOM === ===")
        print("--- Simulation Cn-Wheat : niveau tesselation (RATP)=%i | inter row=%f | iterations=%i | outputs=%s"%(level_tesselation, inter_row, nstep, outfolderpath))
        simulation_caribu(level_tesselation, inter_row, nstep, outfolderpath)
    elif sim==2:
        print("=== === RATP GEOM === ===")
        print("--- Simulation Cn-Wheat : niveau tesselation (RATP)=%i | inter row=%f | iterations=%i | outputs=%s"%(level_tesselation, inter_row, nstep, outfolderpath))
        simulation_ratp(level_tesselation, inter_row, nstep, outfolderpath)
    elif sim==3:
        print("=== === ALL === ===")
        print("--- Simulation Cn-Wheat : niveau tesselation (RATP)=%i | inter row=%f | iterations=%i | outputs=%s"%(level_tesselation, inter_row, nstep, outfolderpath))
        simulation_caribu(level_tesselation, inter_row, nstep, outfolderpath)
        simulation_ratp(level_tesselation, inter_row, nstep, outfolderpath)
    print("=== END ===")
    