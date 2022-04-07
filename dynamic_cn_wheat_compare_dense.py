from src.LightVegeManager import *
from src.FSPMWheat_facade import *

from fspmwheat import caribu_facade

import time
import progressbar

def simulation(SIMULATION_LENGTH, write, outpath):
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
    INPUTS_FOLDER = r'C:\Users\mwoussen\cdd\codes\vegecouplelight\WheatFspm\fspm-wheat\test\inputs'
    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    METEO_FILENAME = 'meteo_Ljutovac2002.csv'
    PHYTOT_FILENAME = 'phytoT.csv'

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
    
    # Paramètres pré-simulation
    model_names = ["fspm-wheat"]
    coordinates = [43.,0,0] # latitude, longitude, timezone
    
    ## Paramètres CARIBU ##
    sun_sky_options = "mix"
    sun_algo="ratp" # soleil avec RATP sundirection dans Shortwave_Balance.f90 
    infinite=True # pas de couvert infini
    caribu_param = [sun_sky_options, sun_algo, infinite, None]
    caribu_sky = [] # ciel turtle à 46 directions
    caribu_rf = [[0.1, 0.05]] # réflectance, transmittance

    ## Paramètres RATP ##
    dv = 0.02 # m
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    ratp_mu = [1.]
    tesselate_level = 5
    distrib_algo = "compute voxel" # "file"
    distrib_option = 45
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
                                     plant_density=250, inter_row=0.15)

        # vérifie si l'itération suivante est encore le jour? et lance le calcul de lumière
        if (t_light % LIGHT_TIMESTEP == 0) and (PARi_next_hours > 0):
            in_scenes = [scene_etendu]
            # recherche des tiges dans les différentes scènes
            id_entity = 0
            id_stems=whichstems_MTG(g, id_entity)

            # ajoute le pattern aux paramètres du modèle
            caribu_param[-1] = domain
            
            # Objet calcul de la lumière
            lghtcaribu = LightVegeManager(in_scenes, 
                                        in_names=model_names,
                                        id_stems=id_stems, 
                                        sky_parameters=caribu_sky,
                                        lightmodel="caribu", lightmodelparam=caribu_param, 
                                        rf=caribu_rf, 
                                        coordinates=coordinates)
            lghtcaribu.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
            lghtcaribu.PAR_update_MTG(g)

            lghtratp = LightVegeManager(in_scenes, 
                                        in_names=model_names,
                                        id_stems=id_stems, 
                                        sky_parameters=ratp_sky,
                                        lightmodel="ratp", lightmodelparam=ratp_parameters, 
                                        rf=ratp_rf, 
                                        coordinates=coordinates)
            lghtratp.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)

            print(lghtcaribu.shapes_outputs)
            print(lghtratp.shapes_outputs)
            
            # write VTK
            if write:
                lghtcaribu.VTKout("outputs/"+outpath+"/cnwheat_caribu_PAR_",t_light)
                lghtratp.VTKout("outputs/"+outpath+"/cnwheat_ratp_PAR_",t_light)
        
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

        # enregistre les valeurs de PAR
        for key,items in g.property("PARa").items():
            para_c.append(items)
            para_r.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["PARa"].values[0])
            shapes.append(key)
            pari.append(lghtcaribu.shapes_outputs[lghtcaribu.shapes_outputs.ShapeId==key]["PARi"].values[0])
            parsun.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["SunlitPAR"].values[0])
            parsha.append(lghtratp.shapes_outputs[lghtratp.shapes_outputs.ShapeId==key]["ShadedPAR"].values[0])
        iter.extend([t_light]*len(g.property("PARa")))
        parin.extend([PARi]*len(g.property("PARa")))

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
                                            t_light)

    myoutputs = {"Shapes" : shapes, 
                    "Iteration" : iter, 
                    "PAR input":parin, 
                    "PARa CARIBU" : para_c, 
                    "PARi CARIBU" : pari, 
                    "PARa RATP" : para_r, 
                    "SunlitPAR" : parsun,
                    "ShadedPAR" : parsha}
    df_myout = pd.DataFrame(myoutputs)
    
    if write: df_myout.to_csv("outputs/dynamic_cn_wheat_csv/"+outpath+".csv")

    print("--- temps execution : ",tot_light)

if __name__ == "__main__":
    nstep=4
    write=True
    outpath = "dynamic_cnwheat_dense_2"
    simulation(nstep,write, outpath)
    print("=== END ===")
    