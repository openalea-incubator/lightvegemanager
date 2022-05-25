from src.LightVegeManager import *
import openalea.plantgl.all as pgl_all

# Paramètres scène en entrée
geom1 = pgl_all.FaceSet([(0,0,2),(2,0,2), (2,2,0),(0,2,0)],[range(4)]) # plaque horizontale
s = pgl_all.Scene([pgl_all.Shape(geom1, pgl_all.Material((250,0,0),1), 888)])

geometry = {}
environment = {}
ratp_parameters = {}
caribu_parameters = {}

# Paramètres pré-simulation
geometry["scenes"] = [s]
geometry["transformations"] = {}
geometry["transformations"]["scenes unit"] = ["m"]
geometry["transformations"]["xyz orientation"] = ["x+ = S"]

environment["sky"] = [2,2,"soc"] # turtle à 46 directions par défaut
environment["diffus"] = True
environment["direct"] = True
environment["reflected"] = False
environment["reflectance coefficients"] = [[0.1, 0.05]]
environment["infinite"] = False

## Paramètres CARIBU ##
caribu_parameters["sun algo"] = "ratp"

def run_simu(day, hour, latitude, longitude, timezone, truesolartime):
    # récupère les données météo
    PARi = 500
    environment["coordinates"] = [latitude, longitude, timezone]

    # Objet calcul de la lumière
    lghtcaribu = LightVegeManager(geometry=geometry,
                                    environment=environment,
                                    lightmodel="caribu",
                                    lightmodel_parameters=caribu_parameters)
    lghtcaribu.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=truesolartime, printsun=True)

run_simu(201, 12., 0., 0., 0., True)
run_simu(201, 12., 40., 0., 0., True)
run_simu(201, 16, 40, 12, 2., True)
run_simu(201, 16, 40, 12, 2., False)
run_simu(347, 16, 46, 0, 0., True)
