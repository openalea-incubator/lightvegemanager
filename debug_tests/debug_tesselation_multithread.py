from src.LightVegeManager import *

from alinea.adel.adel_dynamic import AdelDyn
from alinea.adel.echap_leaf import echap_leaves
import openalea.plantgl.all as pgl_all

import time


def simulation_wheat(multithread):
    # fspm-wheat file name
    INPUTS_FOLDER = r'C:\Users\mwoussen\cdd\codes\vegecouplelight\WheatFspm\fspm-wheat\test\inputs'

    # CN-Wheat parameters
    LEAVES_MODEL = 'Soissons_byleafclass'

    # read adelwheat inputs at t0
    adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model=LEAVES_MODEL))
    g = adel_wheat.load(dir=INPUTS_FOLDER)
    
    # Paramètres pré-simulation
    model_names = ["fspm-wheat"]
    coordinates = [46.,0,0] # latitude, longitude, timezone

    ## Paramètres RATP ##
    dv = 0.02 # m
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    ratp_mu = [1.]
    tesselate_level = 7
    distrib_algo = "compute" # "file"
    distrib_option = 45
    infinite=False
    ratp_parameters = [dx, dy, dz, rs, ratp_mu,tesselate_level, distrib_algo, distrib_option,infinite]
    ratp_sky = [] # ciel turtle à 46 directions
    ratp_rf=[[0., 0]] # leaf reflectance PAR et NIR pour l'entité 0

    in_scenes = [adel_wheat.scene(g)]
    lghtratp = LightVegeManager(in_scenes, 
                                in_names=model_names,
                                sky_parameters=ratp_sky,
                                lightmodel="ratp", lightmodelparam=ratp_parameters, 
                                rf=ratp_rf, 
                                coordinates=coordinates)

def simulation_plaque(multithread):
    # Paramètres pré-simulation
    model_names = ["fspm-wheat"]
    coordinates = [46.,0,0] # latitude, longitude, timezone

    ## Paramètres RATP ##
    dv = 1.1 # m
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    ratp_mu = [1.]
    tesselate_level = 10
    distrib_algo = "compute" # "file"
    distrib_option = 45
    infinite=False
    ratp_parameters = [dx, dy, dz, rs, ratp_mu,tesselate_level, distrib_algo, distrib_option,infinite]
    ratp_sky = [] # ciel turtle à 46 directions
    ratp_rf=[[0., 0]] # leaf reflectance PAR et NIR pour l'entité 0

    geom1 = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,0),(0,2,0)],[range(4)]) # plaque dessous
    geom2 = pgl_all.FaceSet([(0,1,1),(2,1,1), (2,3,1),(0,3,1)],[range(4)]) # plaque dessus, décalée en y (axe est-ouest) vers l'ouest
    s = pgl_all.Scene([pgl_all.Shape(geom1, pgl_all.Material((250,0,0),1), 888),
                        pgl_all.Shape(geom2, pgl_all.Material((0,250,0),1), 999)])
    lghtratp = LightVegeManager([s], 
                                in_names=model_names,
                                sky_parameters=ratp_sky,
                                lightmodel="ratp", lightmodelparam=ratp_parameters, 
                                rf=ratp_rf, 
                                coordinates=coordinates,
                                global_scene_tesselate_level=0,
                                multithread=multithread)

    # VTK de la scène avec x+ = North
    path_out = "outputs/debug_multithreading_tesselation/plaque"
    lghtratp.VTKinit(path_out)

if __name__ == "__main__":
    multithread = 0
    start=time.time()
    simulation_plaque(multithread)
    print("temps : ",time.time()-start)

    multithread = 3
    start=time.time()
    simulation_plaque(multithread)
    print("temps : ",time.time()-start)
    print("=== END ===")
    