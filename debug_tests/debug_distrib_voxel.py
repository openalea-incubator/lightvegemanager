from src.LightVegeManager import *
import openalea.plantgl.all as pgl_all


def simulation(dv):
    geom = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,0),(0,2,0)],[range(4)])
    s = pgl_all.Scene([pgl_all.Shape(geom, pgl_all.Material((250,0,0),1), 888)])
    plaque = [s]

    ## Paramètres RATP ##
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    ratp_rf=[[0., 0]] # leaf reflectance PAR et NIR pour l'entité 0
    ratp_mu = [1.]
    tesselate_level = 7
    distrib_algo = "compute voxel" # "file"
    distrib_option = 45
    infinite=False
    ratp_parameters = [dx, dy, dz, rs, ratp_mu, tesselate_level, distrib_algo, distrib_option, infinite]
    ratp_sky = [0]

    hour = 12   
    coordinates = [40,0,0] # latitude, longitude, timezone
    PARi = 500 # forçage arbitraire en µmol.m-2.s-1
    day = 264 # jour arbitraire (21 septembre)
    print("--- day %i, hour %i"%(day, hour))

    #  RATP
    lghtratp = LightVegeManager(in_scenes=plaque,
                                    sky_parameters=ratp_sky,
                                    lightmodel="ratp", lightmodelparam=ratp_parameters, 
                                    rf=ratp_rf, coordinates=coordinates)
    lghtratp.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
    print(lghtratp.shapes_outputs)

    # VTK de la scène avec x+ = North
    path_out = "outputs/debug_distrib_voxel/ratp_"+str(day)+"_"+str(hour)+"h"
    lghtratp.VTKout(path_out, 1, voxels=True)


if __name__ == "__main__":
    # plaque horizontale
    dv = 1. # m
    simulation(dv)
