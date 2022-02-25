from src.LightVegeManager import *

def situation_1_cm():
    # les scènes en entrée
    s1 = pgl.Scene(r'scenes/scene_1ble.bgeom')
    s2 = pgl.Scene(r'scenes/scene_1luzerne_feuilles.geom')
    
    in_scenes = [s1, s2]
    
    # une liste [rescale, translation, unit] par scène, "none" si pas de transformation
    in_trans = [[3, "none", "m"], 
                ["none", Vector3(-30, -14, 0), "cm"]] 
    in_names = ["fspm-wheat", "l-egume"]
    
    # paramètres de RATP
    dx, dy, dz = 4, 4, 2 # cm
    latitude, longitude, timezone = 48.6, 7.76, 0.
    rs=[0.075, 0.2]
    rf=[0.051, 0.479]
    tesselate_level = 3
    distrib_algo = "file" # "compute"
    distrib_option = "inputs/distrib_leafinclination_exemples.dat" # 9
    ratp_parameters = [dx, dy, dz, latitude, longitude, timezone, rs, tesselate_level, distrib_algo, distrib_option]
    
    # Fichiers météo et données physiques des plantes
    meteopath = "inputs/meteo_exemple.mto"
    vegepath = "inputs/myveg_exemple.veg" 
    
    # Objet calcul de la lumière
    lght = LightVegeManager(in_scenes = in_scenes, 
                    in_transformations=in_trans, 
                    in_names=in_names,
                    lightmodel="ratp",
                    lightmodelparam=ratp_parameters,
                    rf=rf,
                    tesselation=True, 
                    scene_unit="cm")
    print(lght)
    # imprime la géométrie des scènes couplés (triangulation et voxels)
    #lght.VTKinit("outputs/static_simple/lvm_")
    
    # lance le bilan radiatif
    #lght.run(meteopath, vegepath)

    # imprime les résultats sur la triangulation
    #lght.VTKout("outputs/static_simple/lvm_")
    
    # outils calculs des distributions sur la scène
    lght.s5()
    #
    #lght.s2v()

def situation_2_m():
    import pandas as pd

    from alinea.adel.adel_dynamic import AdelDyn
    from alinea.adel.echap_leaf import echap_leaves
    from fspmwheat import elongwheat_facade

    # Path of the directory which contains the inputs of the model
    INPUTS_DIRPATH = r'C:\Users\mwoussen\cdd\codes\vegecouplelight\WheatFspm\fspm-wheat\test\inputs'

    # Number of seconds in 1 hour
    HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    METEO_FILENAME = 'meteo_Ljutovac2002.csv'

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                            ORGANS_INITIAL_STATE_FILENAME,
                            HIDDENZONES_INITIAL_STATE_FILENAME,
                            ELEMENTS_INITIAL_STATE_FILENAME,
                            SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframe = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
        inputs_dataframes[inputs_filename] = inputs_dataframe.where(inputs_dataframe.notnull(), None)

    ELONGWHEAT_TIMESTEP = 1
    
    # read adelwheat inputs at t0
    adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'))
    g = adel_wheat.load(dir=INPUTS_DIRPATH)

    # -- ELONGWHEAT (created first because it is the only facade to add new metamers) --
    # Initial states
    elongwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.HIDDENZONE_INPUTS if i in
                                                                   inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()
    elongwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.ELEMENT_INPUTS if i in
                                                                inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()
    elongwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.AXIS_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    phytoT = os.path.join(INPUTS_DIRPATH, 'phytoT.csv')

    # Update some parameters
    update_cnwheat_parameters = {'SL_ratio_d': 0.25}

    # create empty dataframes to shared data between the models
    shared_axes_inputs_outputs_df = pd.DataFrame()
    shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
    shared_elements_inputs_outputs_df = pd.DataFrame()

    # Facade initialisation
    elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(g,
                                                            ELONGWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                            elongwheat_axes_initial_state,
                                                            elongwheat_hiddenzones_initial_state,
                                                            elongwheat_elements_initial_state,
                                                            shared_axes_inputs_outputs_df,
                                                            shared_hiddenzones_inputs_outputs_df,
                                                            shared_elements_inputs_outputs_df,
                                                            adel_wheat, phytoT, update_cnwheat_parameters)
    
    
    # Update geometry
    adel_wheat.update_geometry(g)
    s1 = adel_wheat.scene(g)
    
    in_scenes = [s1]
    in_names = ["fspm-wheat"]

    # paramètres pour CARIBU
    # [diffuse_model(string), skyazimuths(float), skyzeniths(float), latitude(float), sun_sky_option(string)]
    latitude=48.85
    diffuse_model="soc"
    azimuts=4
    zeniths=5
    sun_sky_options="mix"
    lightmodelparam = [diffuse_model, azimuts, zeniths, latitude, sun_sky_options]
    stems=[(19,0),(34,0)] 

    # Objet calcul de la lumière
    lght = LightVegeManager(in_scenes, in_names=in_names, id_stems=stems, lightmodel="caribu", lightmodelparam=lightmodelparam)
    
    meteo = pd.read_csv(os.path.join(INPUTS_DIRPATH, METEO_FILENAME), index_col='t')
    PARi = meteo.loc[0, ['PARi_MA4']].iloc[0]
    DOY = meteo.loc[0, ['DOY']].iloc[0]
    hour = meteo.loc[0, ['hour']].iloc[0]
    lght.run(PARi=PARi, day=DOY, hour=hour)
    print(lght.shapes_outputs)

    lght.VTKout("outputs\dynamic_cnwheat\mycnwheat",0)


if __name__ == "__main__":
    situation_2_m()