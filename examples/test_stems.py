from lightvegemanager.LVM import LightVegeManager

import openalea.plantgl.all as pgl_all

"""
    Stems processing checking
    -------------------------

    Scene for verification is a set of horizontal planes in one voxel as so:

    -------------
    |   ----    |
    |     ____  |
    | ----      |
    -------------

    We compare the PAR mean on all the planes between CARIBU and RATP
    
    Stems processing
        * CARIBU
            Computing of incident PAR is done with no transmittance, triangles are opaques

        * RATP
            triangles area are divided by 2
"""

def run_print(lghtcaribu, lghtratp, PARi, day, hour, i):
    # calcul
    lghtcaribu.run(energy=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
    lghtratp.run(energy=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)

    # print résultats
    print("=== CARIBU ===")
    print(lghtcaribu.elements_outputs)
    print("=== RATP ===")
    print(lghtratp.elements_outputs)
    print("\n")

    # ecriture des résultats en VTK
    lghtratp.VTK_light("outputs/stems/teststemsratp", i=i)
    lghtcaribu.VTK_light("outputs/stems/teststemscaribu", i=i)

if __name__ == "__main__":
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
    environment["infinite"] = False

    ## Paramètres CARIBU ##
    caribu_parameters["sun algo"] = "caribu"
    caribu_parameters["caribu opt"] = {}
    caribu_parameters["caribu opt"]["par"] = (0.10, 0.05)
    lghtcaribu = LightVegeManager(environment=environment,
                                    lightmodel="caribu",
                                    lightmodel_parameters=caribu_parameters)
    lghtcaribu.build(geometry)

    ## Paramètres RATP ##
    dv = 1. # m
    dx, dy, dz = dv, dv, dv # m
    ratp_parameters["voxel size"] = [dx, dy, dz]
    ratp_parameters["soil reflectance"] = [0., 0.]
    ratp_parameters["mu"] = [1.]
    ratp_parameters["tesselation level"] = 0
    ratp_parameters["angle distrib algo"] = "compute global"
    ratp_parameters["nb angle classes"] = 30
    ratp_parameters["reflectance coefficients"] = [[0.1, 0.05]]
    lghtratp = LightVegeManager(environment=environment,
                                lightmodel="ratp",
                                lightmodel_parameters=ratp_parameters)
    lghtratp.build(geometry)

    # imprime la grille de voxels de RATP
    lghtratp.VTK_nolight("outputs/stems/teststemsratp")

    ## situation 1
    PARi=500
    day=100.
    hour=12
    run_print(lghtcaribu, lghtratp, PARi, day, hour, 1)

    ## situation 2
    hour=17.5
    run_print(lghtcaribu, lghtratp, PARi, day, hour, 2)

    print("=== END ===")
