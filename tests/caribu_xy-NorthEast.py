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


def run_and_vtk(lightvegemanager, folderout, PARi, day, hour) :
    '''run puis ecrit les resultats au VTK avec LightVegeManager
    Args:
        lightvegemanager: après l'appel de init_scene
        folderout: string, chemin du dossier de sortie
        PARi : float, rayonnement en entrée
        day : int
        hour : int
    '''
    print("--- day %i, hour %i"%(day, hour))

    # Calcul du rayonnement
    lightvegemanager.run(energy=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True, printsun=True)

    # VTK de la scène avec x+ = North
    path_out = folderout + "scene_day"+str(day)+"_"+str(hour)+"h"
    lightvegemanager.VTKout(path_out, 1)
    
    # ligne du soleil (rotation de 180° autour de z pour se mettre dans l'espace x+ = North)
    VTKline(Vector3(float(lightvegemanager.sun[0][1][0])*2, float(lightvegemanager.sun[0][1][1])*2, float(lightvegemanager.sun[0][1][2])*2),
            Vector3(-float(lightvegemanager.sun[0][1][0])*2, -float(lightvegemanager.sun[0][1][1])*2, -float(lightvegemanager.sun[0][1][2])*2),
           folderout + "sun_day"+str(day)+"_"+str(hour)+"h.vtk")


if __name__ == "__main__":
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

    # input dict
    geometry = {}
    environment = {}
    caribu_parameters = {}

    # Paramètres pré-simulation
    geometry["scenes"] = [s]

    environment["coordinates"] = [0. ,0. ,0.] # latitude, longitude, timezone
    environment["direct"] = True
    environment["diffus"] = False
    environment["reflected"] = False
    environment["caribu opt"] = {} 
    environment["caribu opt"]["par"] = (0.10, ) # plaques opaques
    environment["infinite"] = False

    ## Paramètres CARIBU ##
    caribu_parameters["sun algo"] = "caribu"

    # Déclaration de l'objet
    lghtcaribu = LightVegeManager(environment=environment,
                                    lightmodel="caribu",
                                    lightmodel_parameters=caribu_parameters)
    # création de la scène dans l'objet
    lghtcaribu.init_scenes(geometry, global_scene_tesselate_level=5)

    # dossier de sortie
    folderout = "outputs/debug_caribu_xy-NorthEast/"
    # forçage arbitraire en µmol.m-2.s-1
    PARi = 500 
    # jour arbitraire (21 septembre)
    day = 264 
    
    # heure matin, soleil arrive de l'est, les 2 plaques reçoivent tout le rayonnement
    hour = 8
    run_and_vtk(lghtcaribu, folderout, PARi, day, hour)

    # heure soir, soleil arrive de l'ouest, la plaque du bas est à l'ombre
    hour = 16
    run_and_vtk(lghtcaribu, folderout, PARi, day, hour)

    print("--- done")