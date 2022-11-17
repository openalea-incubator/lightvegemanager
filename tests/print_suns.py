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

def run_simu(lightvegemanager, day, hour, latitude, longitude, timezone, truesolartime):
    PARi = 500
    environment["coordinates"] = [latitude, longitude, timezone]                               
    lightvegemanager.run(energy=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=truesolartime, printsun=True)

if __name__ == "__main__":
    # Paramètres scène en entrée
    geom = pgl_all.FaceSet([(0,0,2),(2,0,2), (2,2,0),(0,2,0)],[range(4)]) # plaque horizontale
    s = pgl_all.Scene([pgl_all.Shape(geom, pgl_all.Material((250,0,0),1), 888)])

    geometry = {}
    environment = {}
    caribu_parameters = {}

    # Paramètres pré-simulation
    geometry["scenes"] = [s]
    geometry["transformations"] = {}
    geometry["transformations"]["scenes unit"] = ["m"]
    geometry["transformations"]["xyz orientation"] = ["x+ = S"]

    environment["sky"] = [2,2,"soc"]
    environment["diffus"] = True
    environment["direct"] = True
    environment["reflected"] = False
    environment["caribu opt"] = {} 
    environment["caribu opt"]["par"] = (0.10, 0.05)
    environment["infinite"] = False

    ## Paramètres CARIBU ##
    caribu_parameters["sun algo"] = "ratp"

    # Objet calcul de la lumière
    lghtcaribu = LightVegeManager(environment=environment,
                                    lightmodel="caribu",
                                    lightmodel_parameters=caribu_parameters)
    lghtcaribu.init_scenes(geometry)

    run_simu(lghtcaribu, 201, 12., 0., 0., 0., True)
    run_simu(lghtcaribu, 201, 12., 40., 0., 0., True)
    run_simu(lghtcaribu, 201, 16, 40, 12, 2., True)
    run_simu(lghtcaribu, 201, 16, 40, 12, 2., False)
    run_simu(lghtcaribu, 347, 16, 46, 0, 0., True)
