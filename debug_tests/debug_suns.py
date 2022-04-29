from src.LightVegeManager import *
import openalea.plantgl.all as pgl_all

# Paramètres scène en entrée
geom1 = pgl_all.FaceSet([(0,0,2),(2,0,2), (2,2,0),(0,2,0)],[range(4)]) # plaque horizontale
s = pgl_all.Scene([pgl_all.Shape(geom1, pgl_all.Material((250,0,0),1), 888)])
scenes = [s]
model_names = ["fspm-wheat"]
transformations = [["none","none","m","x+ = S"]]
rf = [[0.1, 0.05]] # réflectance, transmittance

# Paramètres CARIBU
sun_sky_options="mix"
sun_algo="ratp"
infinite=False
caribu_param = [sun_sky_options, sun_algo, None]
sky = [2,2,"soc"]

def run_simu(day, hour, latitude, longitude, timezone, truesolartime):
    # récupère les données météo
    PARi = 500
    coordinates = [latitude, longitude, timezone]

    # Objet calcul de la lumière
    lghtcaribu = LightVegeManager(scenes, 
                                    in_names=model_names,
                                    in_transformations=transformations,
                                    sky_parameters=sky,
                                    lightmodel="caribu", lightmodelparam=caribu_param, 
                                    rf=rf, 
                                    coordinates=coordinates)
    lghtcaribu.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=truesolartime, printsun=True)

run_simu(201, 12., 0., 0., 0., True)
run_simu(201, 12., 40., 0., 0., True)
run_simu(201, 16, 40, 12, 2., True)
run_simu(201, 16, 40, 12, 2., False)
run_simu(347, 16, 46, 0, 0., True)
