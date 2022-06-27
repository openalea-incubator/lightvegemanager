import pandas as pd
import numpy as np
import os
import warnings

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
                                update_parameters_all_models=None,
                                N_fertilizations=None):

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                            ORGANS_INITIAL_STATE_FILENAME,
                            HIDDENZONES_INITIAL_STATE_FILENAME,
                            ELEMENTS_INITIAL_STATE_FILENAME,
                            SOILS_INITIAL_STATE_FILENAME):
        #inputs_dataframe = pd.read_csv(os.path.join(INPUTS_FOLDER, inputs_filename))
        inputs_dataframe = pd.read_csv(INPUTS_FOLDER+"/"+inputs_filename)
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

    # Update parameters if specified
    if update_parameters_all_models and 'elongwheat' in update_parameters_all_models:
        update_parameters_elongwheat = update_parameters_all_models['elongwheat']
    else:
        update_parameters_elongwheat = None

    # # Update some parameters
    # update_cnwheat_parameters = {'SL_ratio_d': SL_RATIO_D}

    # Facade initialisation
    elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(g,
                                                            ELONGWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                            elongwheat_axes_initial_state,
                                                            elongwheat_hiddenzones_initial_state,
                                                            elongwheat_elements_initial_state,
                                                            shared_axes_inputs_outputs_df,
                                                            shared_hiddenzones_inputs_outputs_df,
                                                            shared_elements_inputs_outputs_df,
                                                            adel_wheat, phytoT, 
                                                            update_parameters_elongwheat)

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
    # update_cnwheat_parameters = {'AGE_EFFECT_SENESCENCE': AGE_EFFECT_SENESCENCE}

    # Update parameters if specified
    if update_parameters_all_models and 'senescwheat' in update_parameters_all_models:
        update_parameters_senescwheat = update_parameters_all_models['senescwheat']
    else:
        update_parameters_senescwheat = None

    # Facade initialisation
    senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(g,
                                                               SENESCWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               senescwheat_roots_initial_state,
                                                               senescwheat_axes_initial_state,
                                                               senescwheat_elements_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df, 
                                                               update_parameters_senescwheat)

    # -- FARQUHARWHEAT --
    # Initial states
    farquharwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    farquharwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
    update_parameters_farquharwheat = {'SurfacicProteins': False, 'NSC_Retroinhibition': False}
    
    # Facade initialisation
    farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                     farquharwheat_elements_initial_state,
                                                                     farquharwheat_axes_initial_state,
                                                                     shared_elements_inputs_outputs_df,
                                                                     update_parameters_farquharwheat)

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
    # update_cnwheat_parameters = {'VMAX_ROOTS_GROWTH_PREFLO': VMAX_ROOTS_GROWTH_PREFLO}

    # Update parameters if specified
    if update_parameters_all_models and 'growthwheat' in update_parameters_all_models:
        update_parameters_growthwheat = update_parameters_all_models['growthwheat']
    else:
        update_parameters_growthwheat = None

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
                                                               shared_axes_inputs_outputs_df, 
                                                               update_parameters_growthwheat)

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
    # update_cnwheat_parameters = {'roots': {'K_AMINO_ACIDS_EXPORT': K_AMINO_ACIDS_EXPORT,
    #                                        'K_NITRATE_EXPORT': K_NITRATE_EXPORT}}

    # Update parameters if specified
    if update_parameters_all_models and 'cnwheat' in update_parameters_all_models:
        update_parameters_cnwheat = update_parameters_all_models['cnwheat']
    else:
        update_parameters_cnwheat = {}

    # Facade initialisation
    cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                   CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                   PLANT_DENSITY,
                                                   update_parameters_cnwheat,
                                                   cnwheat_organs_initial_state,
                                                   cnwheat_hiddenzones_initial_state,
                                                   cnwheat_elements_initial_state,
                                                   cnwheat_soils_initial_state,
                                                   shared_axes_inputs_outputs_df,
                                                   shared_organs_inputs_outputs_df,
                                                   shared_hiddenzones_inputs_outputs_df,
                                                   shared_elements_inputs_outputs_df,
                                                   shared_soils_inputs_outputs_df)
    
    # Run cnwheat with constant nitrates concentration in the soil if specified
    if N_fertilizations is not None and 'constant_Conc_Nitrates' in N_fertilizations.keys():
        cnwheat_facade_.soils[(1, 'MS')].constant_Conc_Nitrates = True
        cnwheat_facade_.soils[(1, 'MS')].nitrates = N_fertilizations['constant_Conc_Nitrates'] * cnwheat_facade_.soils[(1, 'MS')].volume
    
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
                                        N_fertilizations=None,
                                        tillers_replications=None,
                                        save_geom=False,
                                        GEOM_FOLDER=""):
    # suite de la simu
    for t_senescwheat in range(t_light, t_light + SENESCWHEAT_TIMESTEP, SENESCWHEAT_TIMESTEP):
        # run SenescWheat
        senescwheat_facade_.run()

        # Run the rest of the model if the plant is alive
        for t_farquharwheat in range(t_senescwheat, t_senescwheat + SENESCWHEAT_TIMESTEP, FARQUHARWHEAT_TIMESTEP):
            # get the meteo of the current step
            Ta, ambient_CO2, RH, Ur = meteo.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind']]

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
                        # N fertilization if any
                        if N_fertilizations is not None and len(N_fertilizations) > 0:
                            if t_cnwheat in N_fertilizations.keys():
                                cnwheat_facade_.soils[(1, 'MS')].nitrates += N_fertilizations[t_cnwheat]
                        
                        if t_cnwheat > 0:

                            # run CNWheat
                            Tair = meteo.loc[t_elongwheat, 'air_temperature']
                            Tsoil = meteo.loc[t_elongwheat, 'soil_temperature']
                            cnwheat_facade_.run(Tair, Tsoil, tillers_replications)

                        if save_df:
                            # append outputs at current step to global lists
                            axes_outputs, elements_outputs, hiddenzones_outputs, organs_outputs, soils_outputs = fspmwheat_facade_.build_outputs_df_from_MTG()

                            all_simulation_steps.append(t_cnwheat)
                            axes_all_data_list.append(axes_outputs)
                            organs_all_data_list.append(organs_outputs)
                            hiddenzones_all_data_list.append(hiddenzones_outputs)
                            elements_all_data_list.append(elements_outputs)
                            soils_all_data_list.append(soils_outputs)

def save_df_to_csv(df, outputs_filepath, precision):
    """
    Save pandas dataframes to csv
    :param pandas.DataFrame df: a pandas dataframe to be saved
    :param str outputs_filepath: the name of the CSV file to be saved
    :param int precision: number of decimals in CSV file
    """
    try:
        df.to_csv(outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(precision))
    except IOError as err:
        path, filename = os.path.split(outputs_filepath)
        filename = os.path.splitext(filename)[0]
        newfilename = 'ACTUAL_{}.csv'.format(filename)
        newpath = os.path.join(path, newfilename)
        df.to_csv(newpath, na_rep='NA', index=False, float_format='%.{}f'.format(precision))
        warnings.warn('[{}] {}'.format(err.errno, err.strerror))
        warnings.warn('File will be saved at {}'.format(newpath))

def write_outputs_fspmwheat(OUTPUTS_DIRPATH,
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
                            all_simulation_steps,
                            PRECISION = 4,
                            run_postprocessing=True):
       
    # save des données brutes
    outputs_df_dict = {}
    for (outputs_df_list,
         outputs_filename,
         index_columns) \
            in ((axes_all_data_list, AXES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.AXES_T_INDEXES),
                (organs_all_data_list, ORGANS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ORGANS_T_INDEXES),
                (hiddenzones_all_data_list, HIDDENZONES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES),
                (elements_all_data_list, ELEMENTS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES),
                (soils_all_data_list, SOILS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.SOILS_T_INDEXES)):
        data_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
        outputs_df = pd.concat(outputs_df_list, keys=all_simulation_steps, sort=False)
        outputs_df.reset_index(0, inplace=True)
        outputs_df.rename({'level_0': 't'}, axis=1, inplace=True)
        outputs_df = outputs_df.reindex(index_columns + outputs_df.columns.difference(index_columns).tolist(), axis=1, copy=False)
        outputs_df.fillna(value=np.nan, inplace=True)  # Convert back None to NaN
        save_df_to_csv(outputs_df, data_filepath, PRECISION)
        outputs_file_basename = outputs_filename.split('.')[0]
        outputs_df_dict[outputs_file_basename] = outputs_df.reset_index()
    
    # construit un premier postprocessing
    if run_postprocessing:
        time_grid = list(outputs_df_dict.values())[0].t
        delta_t = (time_grid.loc[1] - time_grid.loc[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        axes_postprocessing_file_basename = AXES_POSTPROCESSING_FILENAME.split('.')[0]
        hiddenzones_postprocessing_file_basename = HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]
        organs_postprocessing_file_basename = ORGANS_POSTPROCESSING_FILENAME.split('.')[0]
        elements_postprocessing_file_basename = ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]
        soils_postprocessing_file_basename = SOILS_POSTPROCESSING_FILENAME.split('.')[0]

        postprocessing_df_dict = {}
        (postprocessing_df_dict[axes_postprocessing_file_basename],
         postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
         postprocessing_df_dict[organs_postprocessing_file_basename],
         postprocessing_df_dict[elements_postprocessing_file_basename],
         postprocessing_df_dict[soils_postprocessing_file_basename]) \
            = cnwheat_facade.CNWheatFacade.postprocessing(axes_outputs_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                          hiddenzone_outputs_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                          organs_outputs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                          elements_outputs_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                          soils_outputs_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]],
                                                          delta_t=delta_t)

        for postprocessing_file_basename, postprocessing_filename in ((axes_postprocessing_file_basename, AXES_POSTPROCESSING_FILENAME),
                                                                      (hiddenzones_postprocessing_file_basename, HIDDENZONES_POSTPROCESSING_FILENAME),
                                                                      (organs_postprocessing_file_basename, ORGANS_POSTPROCESSING_FILENAME),
                                                                      (elements_postprocessing_file_basename, ELEMENTS_POSTPROCESSING_FILENAME),
                                                                      (soils_postprocessing_file_basename, SOILS_POSTPROCESSING_FILENAME)):
            postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
            postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))

def append_outputs_fspmwheat(OUTPUTS_DIRPATH,
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
                            axes_data_list,
                            organs_data_list,
                            hiddenzones_data_list,
                            elements_data_list,
                            soils_data_list,
                            all_simulation_steps,
                            PRECISION = 4,
                            run_postprocessing=True,
                            delta_t="none"):
       
    # save des données brutes
    outputs_df_dict = {}
    for (outputs_df_list,
         outputs_filename,
         index_columns) \
            in ((axes_data_list, AXES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.AXES_T_INDEXES),
                (organs_data_list, ORGANS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ORGANS_T_INDEXES),
                (hiddenzones_data_list, HIDDENZONES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES),
                (elements_data_list, ELEMENTS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES),
                (soils_data_list, SOILS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.SOILS_T_INDEXES)):
        #data_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
        data_filepath = OUTPUTS_DIRPATH+"/"+outputs_filename
        outputs_df = pd.concat(outputs_df_list, keys=all_simulation_steps, sort=False)
        outputs_df.reset_index(0, inplace=True)
        outputs_df.rename({'level_0': 't'}, axis=1, inplace=True)
        outputs_df = outputs_df.reindex(index_columns + outputs_df.columns.difference(index_columns).tolist(), axis=1, copy=False)
        outputs_df.fillna(value=np.nan, inplace=True)  # Convert back None to NaN
        if not os.path.exists(data_filepath):
            outputs_df.to_csv(data_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION), mode='w')
        else:
            outputs_df.to_csv(data_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION), mode='a', header=False)
        
        outputs_file_basename = outputs_filename.split('.')[0]
        outputs_df_dict[outputs_file_basename] = outputs_df.reset_index()
    
    # construit un premier postprocessing
    if run_postprocessing:
        if delta_t=="none":
            time_grid = list(outputs_df_dict.values())[0].t
            delta_t = (time_grid.loc[1] - time_grid.loc[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        axes_postprocessing_file_basename = AXES_POSTPROCESSING_FILENAME.split('.')[0]
        hiddenzones_postprocessing_file_basename = HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]
        organs_postprocessing_file_basename = ORGANS_POSTPROCESSING_FILENAME.split('.')[0]
        elements_postprocessing_file_basename = ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]
        soils_postprocessing_file_basename = SOILS_POSTPROCESSING_FILENAME.split('.')[0]

        postprocessing_df_dict = {}
        (postprocessing_df_dict[axes_postprocessing_file_basename],
         postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
         postprocessing_df_dict[organs_postprocessing_file_basename],
         postprocessing_df_dict[elements_postprocessing_file_basename],
         postprocessing_df_dict[soils_postprocessing_file_basename]) \
            = cnwheat_facade.CNWheatFacade.postprocessing(axes_outputs_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                          hiddenzone_outputs_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                          organs_outputs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                          elements_outputs_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                          soils_outputs_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]],
                                                          delta_t=delta_t)

        for postprocessing_file_basename, postprocessing_filename in ((axes_postprocessing_file_basename, AXES_POSTPROCESSING_FILENAME),
                                                                      (hiddenzones_postprocessing_file_basename, HIDDENZONES_POSTPROCESSING_FILENAME),
                                                                      (organs_postprocessing_file_basename, ORGANS_POSTPROCESSING_FILENAME),
                                                                      (elements_postprocessing_file_basename, ELEMENTS_POSTPROCESSING_FILENAME),
                                                                      (soils_postprocessing_file_basename, SOILS_POSTPROCESSING_FILENAME)):
            #postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
            postprocessing_filepath = POSTPROCESSING_DIRPATH+"/"+postprocessing_filename
            if not os.path.exists(postprocessing_filepath) :
                postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION), mode='w')
            else:
                postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION), mode='a', header=False)


# copie de _create_heterogeneous_canopy dans fspmwheat/caribu_facade pour l'utiliser en dehors de la classe
def create_heterogeneous_canopy_copy(geometrical_model, mtg, nplants=50, var_plant_position=0.03, var_leaf_inclination=0.157, var_leaf_azimut=1.57, var_stem_azimut=0.157,
                                     plant_density=250, inter_row=0.15):
        """
        Duplicate a plant in order to obtain a heterogeneous canopy.

        :param int nplants: the desired number of duplicated plants
        :param float var_plant_position: variability for plant position (m)
        :param float var_leaf_inclination: variability for leaf inclination (rad)
        :param float var_leaf_azimut: variability for leaf azimut (rad)
        :param float var_stem_azimut: variability for stem azimut (rad)

        :return: duplicated heterogenous scene and its domain
        :rtype: openalea.plantgl.all.Scene, (float)
        """
        from alinea.adel.Stand import AgronomicStand
        import openalea.plantgl.all as plantgl
        import random

        # Load scene
        initial_scene = geometrical_model.scene(mtg)

        alea_canopy = pd.DataFrame()

        # Planter
        stand = AgronomicStand(sowing_density=plant_density, plant_density=plant_density, inter_row=inter_row, noise=var_plant_position)
        _, domain, positions, _ = stand.smart_stand(nplants=nplants, at=inter_row, convunit=1)

        random.seed(1234)

        # Built alea table if does not exist yet
        if alea_canopy.empty:
            elements_vid_list = []
            for mtg_plant_vid in mtg.components_iter(mtg.root):
                for mtg_axis_vid in mtg.components_iter(mtg_plant_vid):
                    for mtg_metamer_vid in mtg.components_iter(mtg_axis_vid):
                        for mtg_organ_vid in mtg.components_iter(mtg_metamer_vid):
                            for mtg_element_vid in mtg.components_iter(mtg_organ_vid):
                                if mtg.label(mtg_element_vid) == 'LeafElement1':
                                    elements_vid_list.append(mtg_element_vid)
            elements_vid_df = pd.DataFrame({'vid': elements_vid_list, 'tmp': 1})
            positions_df = pd.DataFrame({'pos': range(len(positions)),
                                         'tmp': 1,
                                         'azimut_leaf': 0,
                                         'inclination_leaf': 0})
            alea = pd.merge(elements_vid_df, positions_df, on=['tmp'])
            alea = alea.drop('tmp', axis=1)
            for vid in elements_vid_list:
                np.random.seed(vid)
                alea.loc[alea['vid'] == vid, 'azimut_leaf'] = np.random.uniform(-var_leaf_azimut, var_leaf_azimut, size=len(positions))
                alea.loc[alea['vid'] == vid, 'inclination_leaf'] = np.random.uniform(-var_leaf_inclination, var_leaf_inclination, size=len(positions))
            alea_canopy = alea

        # Duplication and heterogeneity
        duplicated_scene = plantgl.Scene()
        position_number = 0
        for pos in positions:
            azimut_stem = random.uniform(-var_stem_azimut, var_stem_azimut)
            for shp in initial_scene:
                if mtg.label(shp.id) == 'StemElement':
                    rotated_geometry = plantgl.EulerRotated(azimut_stem, 0, 0, shp.geometry)
                    translated_geometry = plantgl.Translated(plantgl.Vector3(pos), rotated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape
                elif mtg.label(shp.id) == 'LeafElement1':
                    # Add shp.id in alea_canopy if not in yet:
                    if shp.id not in list(alea_canopy['vid']):
                        new_vid_df = pd.DataFrame({'vid': shp.id, 'pos': range(len(positions))})
                        np.random.seed(shp.id)
                        new_vid_df['azimut_leaf'] = np.random.uniform(-var_leaf_azimut, var_leaf_azimut, size=len(positions))
                        new_vid_df['inclination_leaf'] = np.random.uniform(-var_leaf_inclination, var_leaf_inclination, size=len(positions))
                        alea_canopy = alea_canopy.copy().append(new_vid_df, sort=False)
                    # Translation to origin
                    anchor_point = mtg.get_vertex_property(shp.id)['anchor_point']
                    trans_to_origin = plantgl.Translated(-anchor_point, shp.geometry)
                    # Rotation variability
                    azimut = alea_canopy.loc[(alea_canopy.pos == position_number) & (alea_canopy.vid == shp.id), 'azimut_leaf'].values[0]
                    inclination = alea_canopy.loc[(alea_canopy.pos == position_number) & (alea_canopy.vid == shp.id), 'inclination_leaf'].values[0]
                    rotated_geometry = plantgl.EulerRotated(azimut, inclination, 0, trans_to_origin)
                    # Restore leaf base at initial anchor point
                    translated_geometry = plantgl.Translated(anchor_point, rotated_geometry)
                    # Translate leaf to new plant position
                    translated_geometry = plantgl.Translated(pos, translated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape
            position_number += 1

        return duplicated_scene, domain