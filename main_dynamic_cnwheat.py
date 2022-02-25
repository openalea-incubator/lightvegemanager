from src.LightVegeManager import *
from src.FSPMWheat_facade import *

import time
import progressbar

def simu_caribu(SIMULATION_LENGTH, outpath, write):
    # -- SIMULATION PARAMETERS --
    START_TIME = 0
    PLANT_DENSITY = {1: 1}

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
    
    diffuse_model="soc"
    azimuts=9
    zeniths=5

    sun_sky_options="mix"
    sun_algo="ratp"
    lightmodelparam = [sun_sky_options, sun_algo]
    in_names = ["fspm-wheat"]

    rf = [[0.1, 0.05]]

    latitude, longitude, timezone = 46.58, 0.38, 0.
    coordinates = [latitude, longitude, timezone]


    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    para=[]
    pari=[]
    erel=[]
    iter=[]
    shapes=[]
    tot_light = 0.
    for t_light in progressbar.progressbar(range(START_TIME, SIMULATION_LENGTH, LIGHT_TIMESTEP)):
        
        light_start=time.time()
        # récupère les données météo
        PARi = meteo.loc[t_light, ['PARi_MA4']].iloc[0]
        DOY = meteo.loc[t_light, ['DOY']].iloc[0]
        hour = meteo.loc[t_light, ['hour']].iloc[0]
        PARi_next_hours = meteo.loc[range(t_light, t_light + LIGHT_TIMESTEP), ['PARi']].sum().values[0]

        # vérifie si l'itération suivante est encore le jour? et lance le calcul de lumière
        if (t_light % LIGHT_TIMESTEP == 0) and (PARi_next_hours > 0):
            in_scenes = [adel_wheat.scene(g)]
            # recherche des tiges dans les différentes scènes
            id_stems=[]
            for i in range(len(in_scenes)):
                id_stems.extend(whichstems_MTG(g, i))
            
            # Objet calcul de la lumière
            lght = LightVegeManager(in_scenes, in_names=in_names, id_stems=id_stems, 
                                        sky_parameters=[azimuts, zeniths, diffuse_model],
                                        lightmodel="caribu", lightmodelparam=lightmodelparam, 
                                        rf=rf, coordinates=coordinates)
            
            lght.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1")
            lght.PAR_update_MTG(g)

            print(lght.shapes_outputs)
            
            # write VTK
            if write:
                lght.VTKout("outputs/"+outpath+"/cnwheat_dyn_PAR_",t_light)
        
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
            para.append(items)
            shapes.append(key)
        for key,items in g.property("Erel").items():
            erel.append(items)
        iter.extend([t_light]*len(g.property("PARa")))
        pari.extend([PARi]*len(g.property("PARa")))

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

    myoutputs = {"Shapes" : shapes, "PARa" : para, "Erel" : erel, "Iteration" : iter, "PARi":pari}
    df_myout = pd.DataFrame(myoutputs)
    
    if write: df_myout.to_csv("outputs/dynamic_cn_wheat_csv/"+outpath+".csv")

    print("--- CARIBU lumiere temps execution : ",tot_light)
                                
def simu_ratp(SIMULATION_LENGTH, dv, outpath, write):
    # -- SIMULATION PARAMETERS --
    START_TIME = 0
    PLANT_DENSITY = {1: 1}

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

    # paramètres pour le modèle de lumière
    in_scenes=[adel_wheat.scene(g)]
    in_names = ["fspm-wheat"]
    
    # paramètres de RATP
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    rf=[[0., 0]] # leaf reflectance PAR et NIR pour l'entité 0
    tesselate_level = 5
    distrib_algo = "compute" # "file"
    distrib_option = 9
    ratp_parameters = [dx, dy, dz, rs, tesselate_level, distrib_algo, distrib_option]

    sky_parameters=[] # turtle 46 directions par défaut

    latitude, longitude, timezone = 46.58, 0.38, 0.
    coordinates = [latitude, longitude, timezone]

    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    para=[]
    pari=[]
    erel=[]
    iter=[]
    shapes=[]

    tot_light = 0.
    for t_light in progressbar.progressbar(range(START_TIME, SIMULATION_LENGTH, LIGHT_TIMESTEP)):
        # récupère les données météo
        PARi = meteo.loc[t_light, ['PARi_MA4']].iloc[0]
        DOY = meteo.loc[t_light, ['DOY']].iloc[0]
        hour = meteo.loc[t_light, ['hour']].iloc[0]
        PARi_next_hours = meteo.loc[range(t_light, t_light + LIGHT_TIMESTEP), ['PARi']].sum().values[0]

        light_start = time.time()

        # vérifie si l'itération suivante est encore le jour? et lance le calcul de lumière
        if (t_light % LIGHT_TIMESTEP == 0) and (PARi_next_hours > 0):
            in_scenes = [adel_wheat.scene(g)]
            
            # recherche des tiges dans les différentes scènes
            id_stems=[]
            for i in range(len(in_scenes)):
                id_stems.extend(whichstems_MTG(g, i))
            
            # Objet calcul de la lumière
            lght = LightVegeManager(in_scenes, 
                                in_names=in_names,
                                id_stems=id_stems, 
                                lightmodel="ratp", 
                                lightmodelparam=ratp_parameters,
                                sky_parameters=sky_parameters, 
                                rf=rf,
                                coordinates=coordinates)
            
            lght.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1")
            lght.PAR_update_MTG(g)


            print(lght.shapes_outputs)

            # write VTK
            if write:
                #lght.VTKinit("outputs/"+outpath+"/cnwheat_dyn_init_"+str(t_light)+"_")
                lght.VTKout("outputs/"+outpath+"/cnwheat_dyn_PAR_",t_light)
        
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
            para.append(items)
            shapes.append(key)
        for key,items in g.property("Erel").items():
            erel.append(items)
        iter.extend([t_light]*len(g.property("PARa")))
        pari.extend([PARi]*len(g.property("PARa")))

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
    
    print("--- RATP lumiere temps execution : ",tot_light)
    myoutputs = {"Shapes" : shapes, "PARa" : para, "Erel" : erel, "Iteration" : iter, "PARi":pari}
    df_myout = pd.DataFrame(myoutputs)
    if write: df_myout.to_csv("outputs/dynamic_cn_wheat_csv/"+outpath+".csv")

if __name__ == "__main__":
    nstep=1200
    write=True

    out = "dynamic_cnwheat_lvmcaribu"
    
    simu_caribu(nstep,out,write)
    print("=== END CARIBU===")
    
    dv = 0.01 #m
    out = "dynamic_cnwheat_lvmratp_1cm"
    #simu_ratp(nstep, dv, out, write)
    print("=== END RATP 1cm===")
    
    dv = 0.005 #m
    out = "dynamic_cnwheat_lvmratp_5mm"
    #simu_ratp(nstep, dv, out, write)

    print("=== END ===")
    