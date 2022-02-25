import pandas as pd
import os

from alinea.adel.adel_dynamic import AdelDyn
from alinea.adel.echap_leaf import echap_leaves

from fspmwheat import cnwheat_facade
from fspmwheat import elongwheat_facade
from fspmwheat import farquharwheat_facade
from fspmwheat import growthwheat_facade
from fspmwheat import senescwheat_facade
from fspmwheat import fspmwheat_facade

from cnwheat import tools as cnwheat_tools
from cnwheat import simulation as cnwheat_simulation

# Number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

def initialisation_fspmwheat(ELONGWHEAT_TIMESTEP, 
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
                                K_NITRATE_EXPORT):

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                            ORGANS_INITIAL_STATE_FILENAME,
                            HIDDENZONES_INITIAL_STATE_FILENAME,
                            ELEMENTS_INITIAL_STATE_FILENAME,
                            SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframe = pd.read_csv(os.path.join(INPUTS_FOLDER, inputs_filename))
        inputs_dataframes[inputs_filename] = inputs_dataframe.where(inputs_dataframe.notnull(), None)


    # create empty dataframes to shared data between the models
    shared_axes_inputs_outputs_df = pd.DataFrame()
    shared_organs_inputs_outputs_df = pd.DataFrame()
    shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
    shared_elements_inputs_outputs_df = pd.DataFrame()
    shared_soils_inputs_outputs_df = pd.DataFrame()

    # -- ADEL and MTG CONFIGURATION --

    # read adelwheat inputs at t0
    adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model=LEAVES_MODEL))
    g = adel_wheat.load(dir=INPUTS_FOLDER)

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

    phytoT = os.path.join(INPUTS_FOLDER, PHYTOT_FILENAME)

    # Update some parameters
    update_cnwheat_parameters = {'SL_ratio_d': SL_RATIO_D}

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
    update_cnwheat_parameters = {'AGE_EFFECT_SENESCENCE': AGE_EFFECT_SENESCENCE}

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
    update_cnwheat_parameters = {'VMAX_ROOTS_GROWTH_PREFLO': VMAX_ROOTS_GROWTH_PREFLO}

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
    update_cnwheat_parameters = {'roots': {'K_AMINO_ACIDS_EXPORT': K_AMINO_ACIDS_EXPORT,
                                           'K_NITRATE_EXPORT': K_NITRATE_EXPORT}}

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

    meteo = pd.read_csv(os.path.join(INPUTS_FOLDER, METEO_FILENAME), index_col='t')

    return meteo,\
            g, \
            adel_wheat, \
            elongwheat_facade_, \
            growthwheat_facade_, \
            farquharwheat_facade_,   \
            senescwheat_facade_, \
            cnwheat_facade_,\
            fspmwheat_facade_

def iteration_fspmwheat_withoutlighting(meteo,
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
                                        save_df=False,
                                        axes_all_data_list = [],
                                        organs_all_data_list = [],
                                        hiddenzones_all_data_list = [],
                                        elements_all_data_list = [],
                                        soils_all_data_list = [],
                                        all_simulation_steps = [],
                                        save_geom=False,
                                        GEOM_FOLDER=""):
    # suite de la simu
    for t_senescwheat in range(t_light, t_light + LIGHT_TIMESTEP, SENESCWHEAT_TIMESTEP):
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
                
                if save_geom : adel_wheat.save(g, t_elongwheat, GEOM_FOLDER)

                for t_growthwheat in range(t_elongwheat, t_elongwheat + ELONGWHEAT_TIMESTEP, GROWTHWHEAT_TIMESTEP):
                    # run GrowthWheat
                    growthwheat_facade_.run()

                    for t_cnwheat in range(t_growthwheat, t_growthwheat + GROWTHWHEAT_TIMESTEP, CNWHEAT_TIMESTEP):
                        if t_cnwheat > 0:

                            # run CNWheat
                            Tair = meteo.loc[t_elongwheat, 'air_temperature']
                            Tsoil = meteo.loc[t_elongwheat, 'soil_temperature']
                            cnwheat_facade_.run(Tair, Tsoil)

                        if save_df:
                            # append outputs at current step to global lists
                            axes_outputs, elements_outputs, hiddenzones_outputs, organs_outputs, soils_outputs = fspmwheat_facade_.build_outputs_df_from_MTG()

                            all_simulation_steps.append(t_cnwheat)
                            axes_all_data_list.append(axes_outputs)
                            organs_all_data_list.append(organs_outputs)
                            hiddenzones_all_data_list.append(hiddenzones_outputs)
                            elements_all_data_list.append(elements_outputs)
                            soils_all_data_list.append(soils_outputs)

def write_outputs_fspmwheat(OUTPUTS_DIRPATH,
                            DESIRED_AXES_OUTPUTS_FILENAME,
                            DESIRED_ORGANS_OUTPUTS_FILENAME,
                            DESIRED_HIDDENZONES_OUTPUTS_FILENAME,
                            DESIRED_ELEMENTS_OUTPUTS_FILENAME,
                            DESIRED_SOILS_OUTPUTS_FILENAME,
                            ACTUAL_AXES_OUTPUTS_FILENAME,
                            ACTUAL_ORGANS_OUTPUTS_FILENAME,
                            ACTUAL_HIDDENZONES_OUTPUTS_FILENAME,
                            ACTUAL_ELEMENTS_OUTPUTS_FILENAME,
                            ACTUAL_SOILS_OUTPUTS_FILENAME,
                            axes_all_data_list,
                            organs_all_data_list,
                            hiddenzones_all_data_list,
                            elements_all_data_list,
                            soils_all_data_list,
                            all_simulation_steps,
                            PRECISION = 4,
                            overwrite_desired_data=False):
    # compare actual to desired outputs at each scale level (an exception is raised if the test failed)
    for (outputs_df_list,
         desired_outputs_filename,
         actual_outputs_filename,
         index_columns,
         state_variables_names) \
            in ((axes_all_data_list, DESIRED_AXES_OUTPUTS_FILENAME,
                 ACTUAL_AXES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.AXES_T_INDEXES, cnwheat_simulation.Simulation.AXES_STATE),
                (organs_all_data_list, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ORGANS_T_INDEXES, cnwheat_simulation.Simulation.ORGANS_STATE),
                (hiddenzones_all_data_list, DESIRED_HIDDENZONES_OUTPUTS_FILENAME, ACTUAL_HIDDENZONES_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES, cnwheat_simulation.Simulation.HIDDENZONE_STATE),
                (elements_all_data_list, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES, cnwheat_simulation.Simulation.ELEMENTS_STATE),
                (soils_all_data_list, DESIRED_SOILS_OUTPUTS_FILENAME, ACTUAL_SOILS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.SOILS_T_INDEXES, cnwheat_simulation.Simulation.SOILS_STATE)):
        outputs_df = pd.concat(outputs_df_list, keys=all_simulation_steps, sort=False)
        outputs_df.reset_index(0, inplace=True)
        outputs_df.rename({'level_0': 't'}, axis=1, inplace=True)
        outputs_df = outputs_df.reindex(index_columns + outputs_df.columns.difference(index_columns).tolist(), axis=1, copy=False)
        outputs_df = outputs_df.loc[:, state_variables_names]  # compare only the values of the compartments
        cnwheat_tools.compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, desired_outputs_filename,
                                                actual_outputs_filename, precision=PRECISION, overwrite_desired_data=overwrite_desired_data)