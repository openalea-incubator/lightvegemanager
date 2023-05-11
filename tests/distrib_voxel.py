import os
import sys

# si le package est déjà installé
try :
    from  LightVegeManager import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
    from  LightVegeManager import *
    
import openalea.plantgl.all as pgl_all


def simulation(dv):
    geom = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,0),(0,2,0)],[range(4)])
    s = pgl_all.Scene([pgl_all.Shape(geom, pgl_all.Material((250,0,0),1), 888)])

    geometry = {}
    environment = {}
    ratp_parameters = {}

    # Paramètres pré-simulation
    geometry["scenes"] = [s]

    environment["coordinates"] = [40. ,0. ,0.] # latitude, longitude, timezone
    environment["sky"] = "turtle46" # turtle à 46 directions par défaut
    environment["diffus"] = False
    environment["direct"] = True
    environment["reflected"] = False
    environment["infinite"] = False

    ## Paramètres RATP ##
    dx, dy, dz = dv, dv, dv # m
    ratp_parameters["voxel size"] = [dx, dy, dz]
    ratp_parameters["soil reflectance"] = [0., 0.]
    ratp_parameters["reflectance coefficients"] = [[0.1, 0.05]]
    ratp_parameters["mu"] = [1.]
    ratp_parameters["tesselation level"] = 7
    ratp_parameters["angle distrib algo"] = "compute voxel"
    ratp_parameters["nb angle classes"] = 45

    hour = 12   
    PARi = 500 # forçage arbitraire en µmol.m-2.s-1
    day = 264 # jour arbitraire (21 septembre)
    print("--- day %i, hour %i"%(day, hour))

    #  RATP
    lghtratp = LightVegeManager(environment=environment,
                                    lightmodel="ratp",
                                    lightmodel_parameters=ratp_parameters)
    lghtratp.build(geometry)
    lghtratp.run(energy=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
    print(lghtratp.elements_outputs)

    # VTK de la scène avec x+ = North
    path_out = "outputs/debug_distrib_voxel/ratp_"+str(day)+"_"+str(hour)+"h"
    lghtratp.VTK_light(path_out)


if __name__ == "__main__":
    # plaque horizontale
    dv = 1. # m
    simulation(dv)
