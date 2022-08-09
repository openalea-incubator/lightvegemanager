import os
import sys

# si le package est déjà installé
try :
    from src.LightVegeManager import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
    from src.LightVegeManager import *
    
import openalea.plantgl.all as pgl_all


def simulation(geom, hour, situation, direct, diffus, ratp_mu, dv):
    print("--- "+situation)

    geometry = {}
    environment = {}
    ratp_parameters = {}
    caribu_parameters = {}

    # Paramètres pré-simulation
    s = pgl_all.Scene([pgl_all.Shape(geom, pgl_all.Material((250,0,0),1), 888)])
    plaque = [s]
    geometry["scenes"] = plaque

    environment["coordinates"] = [0. ,0. ,0.] # latitude, longitude, timezone
    environment["sky"] = "turtle46" # turtle à 46 directions par défaut
    environment["diffus"] = diffus
    environment["direct"] = direct
    environment["reflected"] = False
    environment["reflectance coefficients"] = [[0.1, 0.05]]
    environment["infinite"] = False

    ## Paramètres CARIBU ##
    caribu_parameters["sun algo"] = "ratp"

    ## Paramètres RATP ##
    dx, dy, dz = dv, dv, dv # m
    ratp_parameters["voxel size"] = [dx, dy, dz]
    ratp_parameters["soil reflectance"] = [0., 0.]
    ratp_parameters["mu"] = ratp_mu
    ratp_parameters["tesselation level"] = 7
    ratp_parameters["angle distrib algo"] = "compute global"
    ratp_parameters["nb angle classes"] = 45


    PARi = 500 # forçage arbitraire en µmol.m-2.s-1
    day = 264 # jour arbitraire (21 septembre)
    print("--- day %i, hour %i"%(day, hour))
    print("\n\t C A R I B U")
    #  CARIBU
    lghtcaribu = LightVegeManager(environment=environment,
                                    lightmodel="caribu", 
                                    lightmodel_parameters=caribu_parameters, 
                                    global_scene_tesselate_level=5)
    lghtcaribu.init_scenes(geometry)                                    
    lghtcaribu.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
    print(lghtcaribu.shapes_outputs)

    # VTK de la scène avec x+ = North
    path_out = "outputs/debug_plaque_ratp-caribu/caribu_"+str(day)+"_"+str(hour)+"h"+"_"+situation
    lghtcaribu.VTKout(path_out, 1)

    #  RATP
    print("\n\t R A T P")
    lghtratp = LightVegeManager(geometry=geometry,
                                environment=environment,
                                lightmodel="ratp", 
                                lightmodel_parameters=ratp_parameters)
    lghtratp.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
    print(lghtratp.shapes_outputs)

    # VTK de la scène avec x+ = North
    path_out = "outputs/debug_plaque_ratp-caribu/ratp_"+str(day)+"_"+str(hour)+"h"+"_"+situation
    lghtratp.VTKout(path_out, 1, voxels=True)


    # ligne du soleil (rotation de 180° autour de z pour se mettre dans l'espace x+ = North)
    VTKline(Vector3(float(-lghtcaribu.sun[0][1][0])*2, float(-lghtcaribu.sun[0][1][1])*2, float(lghtcaribu.sun[0][1][2])*2),
            Vector3(float(lghtcaribu.sun[0][1][0])*2, float(lghtcaribu.sun[0][1][1])*2, -float(lghtcaribu.sun[0][1][2])*2),
            "outputs/debug_plaque_ratp-caribu/sun_day"+str(day)+"_"+str(hour)+"h.vtk")
    print("\n")

if __name__ == "__main__":
    # plaque horizontale
    geom = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,0),(0,2,0)],[range(4)])
    # soleil au zenith + uniquement le direct
    hour = 12
    situation = "plaque_horiz_zenith_direct"
    dv = 10. # m
    ratp_mu = [1.]
    simulation(geom, hour, situation, True, False, ratp_mu, dv)

    # soleil au zenith + diffus + direct
    situation = "plaque_horiz_zenith_directdiffus"
    sun_sky_options = "mix"
    dv = 0.5 # m
    simulation(geom, hour, situation, True, True, ratp_mu, dv)

    # plaque inclinée
    geom = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,1),(0,2,1)],[range(4)])
    # soleil au zenith + direct
    situation = "plaque_incli_zenith_direct_noclump"
    simulation(geom, hour, situation, True, False, ratp_mu, dv)
    # ajoute du clumping
    ratp_mu = [3.]
    situation = "plaque_incli_zenith_direct_withclump"
    simulation(geom, hour, situation, True, False, ratp_mu, dv)

    # direct+diffus
    situation = "plaque_incli_zenith_directdiffus"
    ratp_mu = [1.]
    simulation(geom, hour, situation, True, True, ratp_mu, dv)