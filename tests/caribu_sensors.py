import sys
import os

# si le package est déjà installé
try :
    from src.LightVegeManager import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
    from src.LightVegeManager import *



# Test des sensors dans caribu sur une grille vide
def testsensors():
    print("=== START ===")

    ## INITIALISATION LIGHTVEGEMANAGER
    environment = {}
    caribu_parameters = {}

    # Paramètres pré-simulation
    environment["sky"] = [4, 5, "soc"]
    environment["diffus"] = True
    environment["direct"] = False
    environment["reflected"] = False

    # note : si True -> réception, si False -> atténuation
    environment["infinite"] = True

    environment["caribu opt"] = {} 
    environment["caribu opt"]["par"] = (0.10, 0.07)

    caribu_parameters["sun algo"] = "caribu"
    dxyz = [0.04, 0.04, 0.02]
    nxyz = [10, 10, 1]
    orig=[0.,0.,0.]
    caribu_parameters["sensors"] = ["grid", dxyz, nxyz, orig, "outputs/caribu_sensors/", "vtk"]

    lghtcaribu = LightVegeManager(environment=environment,
                                lightmodel="caribu",
                                lightmodel_parameters=caribu_parameters, 
                                main_unit="m")

    # crée un petit carré
    points = [(0, 0, 0), (0.02, 0, 0), (0.02, 0, 0.01), (0, 0, 0.01)]
    normals = [(0, 1, 0) for i in range(4)]
    indices = [(0, 1, 2, 3)]

    carre = pgl.QuadSet(points, indices, normals, indices)
    carre = pgl.Translated(geometry=carre, translation=(0.2, 0.2, 0.005))

    geometry = {}
    geometry["scenes"] = [pgl.Scene([pgl.Shape(geometry=carre, id=1)])]
    geometry["domain"] = ((0,0), (nxyz[0] * dxyz[0], nxyz[1] * dxyz[1]))

    # Calcul du rayonnement sur un jour arbitraire
    lghtcaribu.init_scenes(geometry)
    lghtcaribu.VTKinit("outputs/caribu_sensors/")
    lghtcaribu.run(energy=1, day=100, hour=12, truesolartime=True, parunit="RG")

    print("=== END ===")

if __name__ == "__main__":
    testsensors()