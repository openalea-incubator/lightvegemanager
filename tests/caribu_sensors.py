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
def testsensors(pgl_scene, folder_vtk_outputs):
    print("=== START ===")

    ## INITIALISATION LIGHTVEGEMANAGER
    environment = {}
    caribu_parameters = {}

    # Paramètres pré-simulation
    environment["sky"] = [4, 5, "soc"]  # ["file", "runscripts/legume/sky_5.data"] [4, 5, "soc"] "turtle46"
    environment["coordinates"] = [44, 12, 1]
    environment["diffus"] = True
    environment["direct"] = False
    environment["reflected"] = False
    environment["infinite"] = False



    # grille de capteurs
    dxyz = [0.04, 0.04, 0.02] # en m
    nxyz = [10, 10, 1] # nombre de capteurs dans chaque axe
    orig = [0.,0.,0.]
    caribu_parameters["sensors"] = ["grid", dxyz, nxyz, orig, folder_vtk_outputs, "vtk"]

    # autres paramètres de CARIBU
    caribu_parameters["debug"] = False
    caribu_parameters["sun algo"] = "caribu"
    caribu_parameters["caribu opt"] = {} 
    caribu_parameters["caribu opt"]["par"] = (0.10, ) # plaque opaque 

    lghtcaribu = LightVegeManager(environment=environment,
                                lightmodel="caribu",
                                lightmodel_parameters=caribu_parameters, 
                                main_unit="m")

    geometry = {}
    geometry["scenes"] = [pgl_scene]
    reduction_domain = 0
    geometry["domain"] = ((orig[0] + reduction_domain, orig[1] + reduction_domain), 
                            (nxyz[0] * dxyz[0] - reduction_domain, nxyz[1] * dxyz[1] - reduction_domain))

    # Calcul du rayonnement sur un jour arbitraire
    day = 100
    hour = 16
    truesolartime = True
    lghtcaribu.build(geometry)
    lghtcaribu.VTK_nolight(folder_vtk_outputs+"init_")
    lghtcaribu.run(energy=500, day=day, hour=hour, truesolartime=truesolartime, parunit="RG")
    lghtcaribu.VTK_light(folder_vtk_outputs)
    # lghtcaribu.VTKsun(folder_vtk_outputs, day, hour, truesolartime)

    print("=== END ===")

if __name__ == "__main__":
    # # Empty
    # points = [(0, 0, -1), (1e-6, 0, -1), (1e-6, 1e-6, -1), (0, 1e-6, -1)]
    # normals = [(0, 0, 1) for i in range(4)]
    # indices = [(0, 1, 2, 3)]
    # carre = pgl.QuadSet(points, indices, normals, indices)
    # s = pgl.Scene([pgl.Shape(geometry=carre, id=1)])
    

    # folder_vtk_outputs = "outputs/caribu_sensors/empty_"
    # testsensors(s, folder_vtk_outputs)

    # Plaque horizontale
    points = [(0, 0, 0), (0.2, 0, 0), (0.2, 0.1, 0), (0, 0.1, 0)]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]
    carre = pgl.QuadSet(points, indices, normals, indices)
    carre = pgl.Translated(geometry=carre, translation=(0.1, 0.2, 0.019))
    s = pgl.Scene([pgl.Shape(geometry=carre, id=1)])

    folder_vtk_outputs = "outputs/caribu_sensors/horiz_"
    testsensors(s, folder_vtk_outputs)

    # # plaque verticale
    # points = [(0, 0, 0), (0.2, 0, 0), (0.2, 0, 0.1), (0, 0, 0.1)]
    # normals = [(0, 1, 0) for i in range(4)]
    # indices = [(0, 1, 2, 3)]
    # carre = pgl.QuadSet(points, indices, normals, indices)
    # carre = pgl.Translated(geometry=carre, translation=(0.1, 0.1, 0.01))
    # s = pgl.Scene([pgl.Shape(geometry=carre, id=1)])

    # folder_vtk_outputs = "outputs/caribu_sensors/vert_"
    # testsensors(s, folder_vtk_outputs)