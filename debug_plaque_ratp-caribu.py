from src.LightVegeManager import *
import openalea.plantgl.all as pgl_all


def simulation(geom, hour, situation, sun_sky_options, ratp_sky, ratp_mu, dv):
    print("--- "+situation)
    s = pgl_all.Scene([pgl_all.Shape(geom, pgl_all.Material((250,0,0),1), 888)])
    plaque = [s]
    model_names = ["fspm-wheat"]
    
    ## Paramètres CARIBU ##
    sun_algo="ratp" # soleil avec RATP sundirection dans Shortwave_Balance.f90 
    infinite=False # pas de couvert infini
    caribu_param = [sun_sky_options, sun_algo, infinite]
    caribu_sky = [] # ciel turtle à 46 directions
    caribu_rf = [[0.1, 0.05]] # réflectance, transmittance

    ## Paramètres RATP ##
    dx, dy, dz = dv, dv, dv # m
    rs=[0., 0] # Soil reflectance in PAR and NIR bands
    ratp_rf=[[0., 0]] # leaf reflectance PAR et NIR pour l'entité 0
    tesselate_level = 7
    distrib_algo = "compute global" # "file"
    distrib_option = 45
    infinite=False
    ratp_parameters = [dx, dy, dz, rs, ratp_mu, tesselate_level, distrib_algo, distrib_option, infinite]

    coordinates = [0,0,0] # latitude, longitude, timezone
    PARi = 500 # forçage arbitraire en µmol.m-2.s-1
    day = 264 # jour arbitraire (21 septembre)
    print("--- day %i, hour %i"%(day, hour))
    print("\t C A R I B U")
    #  CARIBU
    lghtcaribu = LightVegeManager(in_scenes=plaque,
                                    sky_parameters=caribu_sky,
                                    lightmodel="caribu", lightmodelparam=caribu_param, 
                                    rf=caribu_rf, 
                                    coordinates=coordinates,
                                    global_scene_tesselate_level=5)
    lghtcaribu.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
    print(lghtcaribu.shapes_outputs)

    # VTK de la scène avec x+ = North
    path_out = "outputs/debug_plaque_ratp-caribu/caribu_"+str(day)+"_"+str(hour)+"h"+"_"+situation
    lghtcaribu.VTKout(path_out, 1)

    #  RATP
    lghtratp = LightVegeManager(in_scenes=plaque,
                                    sky_parameters=ratp_sky,
                                    lightmodel="ratp", lightmodelparam=ratp_parameters, 
                                    rf=ratp_rf, coordinates=coordinates)
    lghtratp.run(PARi=PARi, day=day, hour=hour, parunit="micromol.m-2.s-1", truesolartime=True)
    print(lghtratp.shapes_outputs)

    # VTK de la scène avec x+ = North
    path_out = "outputs/debug_plaque_ratp-caribu/ratp_"+str(day)+"_"+str(hour)+"h"+"_"+situation
    lghtratp.VTKout(path_out, 1, voxels=True)


    # ligne du soleil (rotation de 180° autour de z pour se mettre dans l'espace x+ = North)
    VTKline(Vector3(float(-lghtcaribu.sun[0][1][0])*2, float(-lghtcaribu.sun[0][1][1])*2, float(lghtcaribu.sun[0][1][2])*2),
            Vector3(float(lghtcaribu.sun[0][1][0])*2, float(lghtcaribu.sun[0][1][1])*2, -float(lghtcaribu.sun[0][1][2])*2),
            "outputs/debug_plaque_ratp-caribu/sun_day"+str(day)+"_"+str(hour)+"h.vtk")
    print("\n")

if __name__ == "__main__":
    # plaque horizontale
    geom = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,0),(0,2,0)],[range(4)])
    # soleil au zenith + uniquement le direct
    hour = 12
    situation = "plaque_horiz_zenith_direct"
    sun_sky_options = "sun"
    dv = 10. # m
    ratp_sky = [0]
    ratp_mu = [1.]
    simulation(geom, hour, situation, sun_sky_options, ratp_sky, ratp_mu, dv)

    # soleil au zenith + diffus + direct
    situation = "plaque_horiz_zenith_directdiffus"
    sun_sky_options = "mix"
    dv = 0.5 # m
    ratp_sky = []
    simulation(geom, hour, situation, sun_sky_options, ratp_sky, ratp_mu, dv)

    # plaque inclinée
    geom = pgl_all.FaceSet([(0,0,0),(2,0,0), (2,2,1),(0,2,1)],[range(4)])
    # soleil au zenith + direct
    situation = "plaque_incli_zenith_direct_noclump"
    sun_sky_options = "sun"
    ratp_sky = [0]
    simulation(geom, hour, situation, sun_sky_options, ratp_sky, ratp_mu, dv)
    # ajoute du clumping
    ratp_mu = [3.]
    situation = "plaque_incli_zenith_direct_withclump"
    simulation(geom, hour, situation, sun_sky_options, ratp_sky, ratp_mu, dv)

    # direct+diffus
    situation = "plaque_incli_zenith_directdiffus"
    sun_sky_options = "mix"
    ratp_sky = []
    ratp_mu = [1.]
    simulation(geom, hour, situation, sun_sky_options, ratp_sky, ratp_mu, dv)