from src.LightVegeManager import *
import openalea.plantgl.all as pgl_all

"""
vérification de la prise en charge des éléments tiges
plusieurs plaques horizontales dans un voxel comme suit (censé correspondre à l'hypothèse turbid medium) : 
    -------------
    |   ----    |
    |     ____  |
    | ----      |
    -------------

On compare la moyenne du PAR sur les plaques entre CARIBU et RATP

"""

# scène géométrique
# NE
geom1 = pgl_all.FaceSet([(0.1,0.1,0.96),(0.1,0.6,0.96), (0.6,0.6,0.96),(0.6,0.1,0.96)],[range(4)])

geom2 = pgl_all.FaceSet([(0.1,0.5,0.806),(0.1,0.9,0.806), (0.6,0.9,0.806),(0.6,0.5,0.806)],[range(4)])

geom3 = pgl_all.FaceSet([(0.5,0.1,0.646),(0.5,0.7,0.646), (0.9,0.7,0.646),(0.9,0.1,0.646)],[range(4)])

geom4 = pgl_all.FaceSet([(0.3,0.5,0.486),(0.3,0.9,0.486), (0.8,0.9,0.486),(0.8,0.5,0.486)],[range(4)])

geom5 = pgl_all.FaceSet([(0.3,0.1,0.32),(0.3,0.5,0.32), (0.8,0.5,0.32),(0.8,0.1,0.32)],[range(4)])

geom6 = pgl_all.FaceSet([(0.2,0.4,0.16),(0.2,0.9,0.16), (0.6,0.9,0.16),(0.6,0.4,0.16)],[range(4)])

s = pgl_all.Scene([pgl_all.Shape(geom1, pgl_all.Material((250,0,0),1), 888),
                    pgl_all.Shape(geom2, pgl_all.Material((250,0,0),1), 888),
                    pgl_all.Shape(geom3, pgl_all.Material((250,0,0),1), 888),
                    pgl_all.Shape(geom4, pgl_all.Material((250,0,0),1), 888),
                    pgl_all.Shape(geom5, pgl_all.Material((250,0,0),1), 888),
                    pgl_all.Shape(geom6, pgl_all.Material((250,0,0),1), 888)])

geometry = {}
environment = {}
ratp_parameters = {}
caribu_parameters = {}

# Paramètres pré-simulation
geometry["scenes"] = [s]
geometry["stems id"] = [(888, 0)]

environment["coordinates"] = [0., 0., 0.] # latitude, longitude, timezone
environment["sky"] = "turtle46" # turtle à 46 directions par défaut
environment["diffus"] = False
environment["direct"] = True
environment["reflected"] = False
environment["reflectance coefficients"] = [[0.1, 0.05]]
environment["infinite"] = False

## Paramètres CARIBU ##
caribu_parameters["sun algo"] = "caribu"
lghtcaribu = LightVegeManager(environment=environment,
                                lightmodel="caribu",
                                lightmodel_parameters=caribu_parameters)
lghtcaribu.init_scenes(geometry)

## Paramètres RATP ##
dv = 1. # m
dx, dy, dz = dv, dv, dv # m
ratp_parameters["voxel size"] = [dx, dy, dz]
ratp_parameters["soil reflectance"] = [0., 0.]
ratp_parameters["mu"] = [1.]
ratp_parameters["tesselation level"] = 0
ratp_parameters["angle distrib algo"] = "compute global"
ratp_parameters["nb angle classes"] = 30
lghtratp = LightVegeManager(environment=environment,
                            lightmodel="ratp",
                            lightmodel_parameters=ratp_parameters)
lghtratp.init_scenes(geometry)

# imprime la scène
lghtratp.VTKinit("outputs/stems/")

# situation
PARi=500
DOY=201
hour=12.

# calcul
lghtcaribu.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True, printsun=True)
lghtratp.run(PARi=PARi, day=DOY, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)

# résultats
print("=== CARIBU ===")
print(lghtcaribu.shapes_outputs)
print("=== RATP ===")
print(lghtratp.shapes_outputs)

lghtcaribu.VTKout("outputs/stems/", iteration=0)
