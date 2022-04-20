from src.LightVegeManager import *
import openalea.plantgl.all as pgl_all

## Paramètres scène en entrée ##
# plaque l'une au dessus de l'autre
#       
#       
#                           -------
# East ---> West (==y)  -------
#
geom1 = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,0),(0,2,0)],[range(4)]) # plaque dessous
geom2 = pgl_all.FaceSet([(0,1,1),(2,1,1), (2,3,1),(0,3,1)],[range(4)]) # plaque dessus, décalée en y (axe est-ouest) vers l'ouest
s = pgl_all.Scene([pgl_all.Shape(geom1, pgl_all.Material((250,0,0),1), 888),
                    pgl_all.Shape(geom2, pgl_all.Material((0,250,0),1), 999)])
in_scenes = [s]
model_names = ["fspm-wheat"]
rf = [[0.1, 0.05]] # réflectance, transmittance

## Paramètres CARIBU ##
sun_sky_options="sun" # calcule que le direct
sun_algo="caribu" # soleil avec GenSun 
infinite=False # pas de couvert infini
caribu_param = [sun_sky_options, sun_algo, None]
sky = [] # ciel turtle à 46 directions

coordinates = [0,0,0] # latitude, longitude, timezone
PARi = 500 # forçage arbitraire en µmol.m-2.s-1
day = 264 # jour arbitraire (21 septembre)

# Déclaration de l'objet
lghtcaribu = LightVegeManager(in_scenes=in_scenes,
                                # in_transformations=None, pas de transformations
                                in_names=model_names,
                                sky_parameters=sky,
                                lightmodel="caribu", lightmodelparam=caribu_param, 
                                rf=rf, 
                                coordinates=coordinates,
                                global_scene_tesselate_level=5)

# heure matin, soleil arrive de l'est, les 2 plaques reçoivent tout le rayonnement
hour = 8
print("--- day %i, hour %i"%(day, hour))
lghtcaribu.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True, printsun=True)

# VTK de la scène avec x+ = North
path_out = "outputs/debug_caribu_xy-NorthEast/scene_day"+str(day)+"_"+str(hour)+"h"
lghtcaribu.VTKout(path_out, 1)
# ligne du soleil (rotation de 180° autour de z pour se mettre dans l'espace x+ = North)
VTKline(Vector3(float(lghtcaribu.sun[0][1][0])*2, float(lghtcaribu.sun[0][1][1])*2, float(lghtcaribu.sun[0][1][2])*2),
        Vector3(-float(lghtcaribu.sun[0][1][0])*2, -float(lghtcaribu.sun[0][1][1])*2, -float(lghtcaribu.sun[0][1][2])*2),
        "outputs/debug_caribu_xy-NorthEast/sun_day"+str(day)+"_"+str(hour)+"h.vtk")

# heure soir, soleil arrive de l'ouest, la plaque du bas est à l'ombre
day = 264
hour = 16
print("--- day %i, hour %i"%(day, hour))
lghtcaribu.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True, printsun=True)

# VTK de la scène avec x+ = North
path_out = "outputs/debug_caribu_xy-NorthEast/scene_day"+str(day)+"_"+str(hour)+"h"
lghtcaribu.VTKout(path_out, 1)
# ligne du soleil (rotation de 180° autour de z pour se mettre dans l'espace x+ = North)
VTKline(Vector3(float(lghtcaribu.sun[0][1][0])*2, float(lghtcaribu.sun[0][1][1])*2, float(lghtcaribu.sun[0][1][2])*2),
        Vector3(-float(lghtcaribu.sun[0][1][0])*2, -float(lghtcaribu.sun[0][1][1])*2, -float(lghtcaribu.sun[0][1][2])*2),
        "outputs/debug_caribu_xy-NorthEast/sun_day"+str(day)+"_"+str(hour)+"h.vtk")
print("--- done")