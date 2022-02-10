import elongwheat
from src.LightVegeManager import *

import pandas as pd
from alinea.adel.adel_dynamic import AdelDyn
from alinea.adel.echap_leaf import echap_leaves

from fspmwheat import caribu_facade
from fspmwheat import cnwheat_facade
from fspmwheat import elongwheat_facade
from fspmwheat import farquharwheat_facade
from fspmwheat import growthwheat_facade
from fspmwheat import senescwheat_facade
from fspmwheat import fspmwheat_facade

import progressbar

# Number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

def initialisation_fspmwheat(ELONGWHEAT_TIMESTEP, 
                                GROWTHWHEAT_TIMESTEP, 
                                SENESCWHEAT_TIMESTEP, 
                                CNWHEAT_TIMESTEP,
                                PLANT_DENSITY):
    # Path of the directory which contains the inputs of the model
    INPUTS_DIRPATH = r'C:\Users\mwoussen\cdd\codes\vegecouplelight\WheatFspm\fspm-wheat\test\inputs'

    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                            ORGANS_INITIAL_STATE_FILENAME,
                            HIDDENZONES_INITIAL_STATE_FILENAME,
                            ELEMENTS_INITIAL_STATE_FILENAME,
                            SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframe = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
        inputs_dataframes[inputs_filename] = inputs_dataframe.where(inputs_dataframe.notnull(), None)


    # create empty dataframes to shared data between the models
    shared_axes_inputs_outputs_df = pd.DataFrame()
    shared_organs_inputs_outputs_df = pd.DataFrame()
    shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
    shared_elements_inputs_outputs_df = pd.DataFrame()
    shared_soils_inputs_outputs_df = pd.DataFrame()

    # -- ADEL and MTG CONFIGURATION --

    # read adelwheat inputs at t0
    adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'))
    g = adel_wheat.load(dir=INPUTS_DIRPATH)

    # ---------------------------------------------
    # ----- CONFIGURATION OF THE FACADES -------
    # ---------------------------------------------

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

    # -- SENESCWHEAT --
    # Initial states
    senescwheat_roots_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
        senescwheat_facade.converter.ROOTS_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_ROOTS_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    senescwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        senescwheat_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    senescwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        senescwheat_facade.converter.AXES_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Update some parameters
    update_cnwheat_parameters = {'AGE_EFFECT_SENESCENCE': 10000}

    # Facade initialisation
    senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(g,
                                                               SENESCWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               senescwheat_roots_initial_state,
                                                               senescwheat_axes_initial_state,
                                                               senescwheat_elements_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df, update_cnwheat_parameters)

    # -- FARQUHARWHEAT --
    # Initial states
    farquharwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    farquharwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                     farquharwheat_elements_initial_state,
                                                                     farquharwheat_axes_initial_state,
                                                                     shared_elements_inputs_outputs_df)

    # -- GROWTHWHEAT --
    # Initial states
    growthwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.HIDDENZONE_INPUTS if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.ELEMENT_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_root_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
        growthwheat_facade.converter.ROOT_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.ROOT_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.AXIS_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Update some parameters
    update_cnwheat_parameters = {'VMAX_ROOTS_GROWTH_PREFLO': 0.02885625}

    # Facade initialisation
    growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(g,
                                                               GROWTHWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               growthwheat_hiddenzones_initial_state,
                                                               growthwheat_elements_initial_state,
                                                               growthwheat_root_initial_state,
                                                               growthwheat_axes_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_hiddenzones_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df, update_cnwheat_parameters)

    # -- CNWHEAT --
    # Initial states
    cnwheat_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ORGANS_VARIABLES if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.HIDDENZONE_VARIABLES if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ELEMENTS_VARIABLES if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.SOILS_VARIABLES if i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

    # Update some parameters
    update_cnwheat_parameters = {'roots': {'K_AMINO_ACIDS_EXPORT': 3E-5,
                                           'K_NITRATE_EXPORT': 1E-6}}

    # Facade initialisation
    cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                   CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                   PLANT_DENSITY,
                                                   update_cnwheat_parameters,
                                                   cnwheat_organs_initial_state,
                                                   cnwheat_hiddenzones_initial_state,
                                                   cnwheat_elements_initial_state,
                                                   cnwheat_soils_initial_state,
                                                   shared_axes_inputs_outputs_df,
                                                   shared_organs_inputs_outputs_df,
                                                   shared_hiddenzones_inputs_outputs_df,
                                                   shared_elements_inputs_outputs_df,
                                                   shared_soils_inputs_outputs_df)
    # -- FSPMWHEAT --
    # Facade initialisation
    fspmwheat_facade_ = fspmwheat_facade.FSPMWheatFacade(g)

    # Update geometry
    adel_wheat.update_geometry(g)

    return g, \
            adel_wheat, \
            elongwheat_facade_, \
            growthwheat_facade_, \
            farquharwheat_facade_,   \
            senescwheat_facade_, \
            cnwheat_facade_,\
            fspmwheat_facade_

if __name__ == "__main__":
    # -- SIMULATION PARAMETERS --
    START_TIME = 0
    SIMULATION_LENGTH = 1200
    PLANT_DENSITY = {1: 1}

    # define the time step in hours for each simulator
    CARIBU_TIMESTEP = 4
    SENESCWHEAT_TIMESTEP = 2
    FARQUHARWHEAT_TIMESTEP = 2
    ELONGWHEAT_TIMESTEP = 1
    GROWTHWHEAT_TIMESTEP = 1
    CNWHEAT_TIMESTEP = 1
    
    g, adel_wheat, \
        elongwheat_facade_, growthwheat_facade_, \
        farquharwheat_facade_, senescwheat_facade_, \
        cnwheat_facade_, fspmwheat_facade_ = initialisation_fspmwheat(ELONGWHEAT_TIMESTEP, 
                                                                        GROWTHWHEAT_TIMESTEP, 
                                                                        SENESCWHEAT_TIMESTEP, 
                                                                        CNWHEAT_TIMESTEP,
                                                                        PLANT_DENSITY)


    # paramètres pour LVM
    # [diffuse_model(string), skyazimuths(float), skyzeniths(float), latitude(float), sun_sky_option(string)]
    latitude=48.85
    diffuse_model="soc"
    azimuts=4
    zeniths=5
    sun_sky_options="mix"
    lightmodelparam = [diffuse_model, azimuts, zeniths, latitude, sun_sky_options]
    in_names = ["fspm-wheat"]

    # Path of the directory which contains the inputs of the model
    INPUTS_DIRPATH = r'C:\Users\mwoussen\cdd\codes\vegecouplelight\WheatFspm\fspm-wheat\test\inputs'
    METEO_FILENAME = 'meteo_Ljutovac2002.csv'
    # Name of the CSV files which contains the meteo data
    meteo = pd.read_csv(os.path.join(INPUTS_DIRPATH, METEO_FILENAME), index_col='t')
    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------
    
    para=[]
    pari=[]
    erel=[]
    iter=[]
    shapes=[]

    for t_caribu in progressbar.progressbar(range(START_TIME, SIMULATION_LENGTH, CARIBU_TIMESTEP)):
        # run Caribu
        PARi = meteo.loc[t_caribu, ['PARi_MA4']].iloc[0]
        DOY = meteo.loc[t_caribu, ['DOY']].iloc[0]
        hour = meteo.loc[t_caribu, ['hour']].iloc[0]
        PARi_next_hours = meteo.loc[range(t_caribu, t_caribu + CARIBU_TIMESTEP), ['PARi']].sum().values[0]

        # vérifie si l'itération suivante est encore le jour? et lance le calcul de lumière
        if (t_caribu % CARIBU_TIMESTEP == 0) and (PARi_next_hours > 0):
            in_scenes = [adel_wheat.scene(g)]
            # recherche des tiges dans les différentes scènes
            id_stems=[]
            for i in range(len(in_scenes)):
                id_stems.extend(whichstems_MTG(g, i))
            
            # Objet calcul de la lumière
            lght = LightVegeManager(in_scenes, in_names=in_names, id_stems=id_stems, lightmodel="caribu", lightmodelparam=lightmodelparam)
            lght.run(PARi=PARi, day=DOY, hour=hour)
            lght.PAR_update_MTG(g)
            # write VTK
            lght.VTKout("outputs/dynamic_cnwheat/cnwheat_dyn_PAR_",t_caribu)
        
        # sinon on copie le PAR de l'itération d'avant
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

        # enregistre les valeurs de PAR
        for key,items in g.property("PARa").items():
            para.append(items)
            shapes.append(key)
        for key,items in g.property("Erel").items():
            erel.append(items)
        iter.extend([t_caribu]*len(g.property("PARa")))
        pari.extend([PARi]*len(g.property("PARa")))

        # suite de la simu
        for t_senescwheat in range(t_caribu, t_caribu + CARIBU_TIMESTEP, SENESCWHEAT_TIMESTEP):
            # run SenescWheat
            senescwheat_facade_.run()

            # Run the rest of the model if the plant is alive
            for t_farquharwheat in range(t_senescwheat, t_senescwheat + SENESCWHEAT_TIMESTEP, FARQUHARWHEAT_TIMESTEP):
                # get the meteo of the current step
                Ta, ambient_CO2, RH, Ur = meteo.loc[t_farquharwheat, ['air_temperature_MA2', 'ambient_CO2_MA2', 'humidity_MA2', 'Wind_MA2']]

                # run FarquharWheat
                farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)

                for t_elongwheat in range(t_farquharwheat, t_farquharwheat + FARQUHARWHEAT_TIMESTEP, ELONGWHEAT_TIMESTEP):
                    # run ElongWheat
                    Tair, Tsoil = meteo.loc[t_elongwheat, ['air_temperature', 'soil_temperature']]
                    elongwheat_facade_.run(Tair, Tsoil, option_static=False)

                    # Update geometry
                    adel_wheat.update_geometry(g)
                    adel_wheat.save(g, t_elongwheat, "outputs/dynamic_cnwheat/")

                    for t_growthwheat in range(t_elongwheat, t_elongwheat + ELONGWHEAT_TIMESTEP, GROWTHWHEAT_TIMESTEP):
                        # run GrowthWheat
                        growthwheat_facade_.run()

                        for t_cnwheat in range(t_growthwheat, t_growthwheat + GROWTHWHEAT_TIMESTEP, CNWHEAT_TIMESTEP):
                            if t_cnwheat > 0:

                                # run CNWheat
                                Tair = meteo.loc[t_elongwheat, 'air_temperature']
                                Tsoil = meteo.loc[t_elongwheat, 'soil_temperature']
                                cnwheat_facade_.run(Tair, Tsoil)

    myoutputs = {"Shapes" : shapes, "PARa" : para, "Erel" : erel, "Iteration" : iter, "PARi":pari}
    df_myout = pd.DataFrame(myoutputs)
    df_myout.to_csv("outputs/dynamic_cnwheat/my_wheat_resuts.csv")