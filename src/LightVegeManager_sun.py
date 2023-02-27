'''
construction du soleil
transfert à CARIBU
soit algo de calcul de CARIBU
soit algo de calcul 
'''
import math

from PyRATP.pyratp import pyratp
from alinea.caribu.sky_tools import Gensun
from alinea.caribu.sky_tools import GetLightsSun
from alinea.caribu.sky_tools import Sun

def ratp_sun(day, hour, coordinates, truesolartime) :                   
    az, ele = 5.,9. # variables fantômes (non récupérées)
    pyratp.shortwave_balance.sundirection(ele, az, 
                                            coordinates[0], 
                                            coordinates[1], 
                                            coordinates[2], 
                                            day, hour, truesolartime)            
    degtorad = math.pi/180

    azrad = pyratp.shortwave_balance.azdeg
    # peut avoir un nan à 12h 
    if math.isnan(azrad) and hour==12. : 
        # critère empirique basé sur l'algo de CARIBU (GenSun)
        if coordinates[0] >= 21.11:
            azrad = 0.
        else:
            azrad = 180.
    # passage de l'azimuth en South=0° clockwise (RATP) à North=0° 
    # clockwise (CARIBU.Gensun)
    azrad = -azrad * degtorad

    zenirad = pyratp.shortwave_balance.hdeg * degtorad

    # le vecteur pointe du ciel vers le sol
    sunx = suny = math.cos(zenirad)
    sunx *= math.cos(azrad)
    suny *= math.sin(azrad)
    sunz = -math.sin(zenirad)
    
    return (sunx, suny, sunz)

def caribu_sun(day, hour, coordinates, truesolartime) :
    # si l'heure est locale, calcul de l'heure solaire 
    # (algorithme pris dans RATP, sundirection: shortwave_balance.f90)
    if not truesolartime:
        om =0.017202*(day-3.244)
        teta=om+0.03344*math.sin(om)*(1+0.021*math.cos(om))-1.3526
        tphi=0.91747*math.sin(teta)/math.cos(teta)
        dphi=math.atan(tphi)-om+1.3526
        if dphi+1. <= 0: dphi=(dphi+1+1000.*math.pi % math.pi)-1.
        eqntime=dphi*229.2
        hour =hour+coordinates[2]+coordinates[1]/15.-eqntime/60. % 24.

    # Conversion de la latitude en radian!
    sun = Gensun.Gensun()(1., day, hour, coordinates[0]*math.pi/180)
    sun = GetLightsSun.GetLightsSun(sun)
    sun_str_split = sun.split(' ')
    
    return (float(sun_str_split[1]), 
            float(sun_str_split[2]), 
            float(sun_str_split[3]))

def print_sun(day, hour, coordinates, truesolartime) :
    """Méthode qui imprime les coordonnées du soleil (pour le débuggage)
        séparée pour la lisibilité dans run(...)
    """
    print("---\t SUN COORDONATES\t ---")
    print("--- Convention x+ = North, vector from sky to floor")
    print("--- azimut: south clockwise E = -90° W = 90°\t zenith: \
                                zenith = 0° horizon = 90°")
    print("--- day: %i \t hour: %i \t latitude: %0.2f °" % 
                                                (day, hour, coordinates[0]))
    
    if truesolartime : 
        caribuhour=hour
        print("--- true solar time")
    
    else :
        om =0.017202*(day-3.244)
        teta=om+0.03344*math.sin(om)*(1+0.021*math.cos(om))-1.3526
        tphi=0.91747*math.sin(teta)/math.cos(teta)
        dphi=math.atan(tphi)-om+1.3526
        if dphi+1. <= 0: dphi=(dphi+1+1000.*math.pi % math.pi)-1.
        eqntime=dphi*229.2
        caribuhour =hour+coordinates[2]+coordinates[1]/15.-eqntime/60. % 24. 
        print("--- local time, true solar time is %0.2f"%(caribuhour))
    
    # CARIBU et RATP sont dans la même convention de coordonnées, x+ = North
    print("--- RATP ---")
    az, ele = 5.,9. # variables fantômes (non récupérées)
    pyratp.shortwave_balance.sundirection(ele, az, 
                                            coordinates[0], 
                                            coordinates[1], 
                                            coordinates[2], 
                                            day, hour, 
                                            truesolartime)        
    degtorad = math.pi/180
    azrad = pyratp.shortwave_balance.azdeg
    if math.isnan(azrad) and hour==12. : 
        if coordinates[0] >= 21.11:
            azrad = 0.
        else:
            azrad = 180.
    azrad = -azrad * degtorad
    elerad = pyratp.shortwave_balance.hdeg * degtorad

    sunx = suny = math.cos(elerad)
    sunx *= math.cos(azrad)
    suny *= math.sin(azrad)
    sunz = -math.sin(elerad)

    print("\t azimut: %0.3f \t zenith: %0.3f" % \
                        (-azrad*180/math.pi, pyratp.shortwave_balance.hdeg))
    print("\t x: %0.3f \t y: %0.3f \t z: %0.3f" % (sunx, suny, sunz))
    
    print("--- CARIBU ---")
    suncaribu = Sun.Sun()
    suncaribu._set_pos_astro(day, caribuhour, coordinates[0]*math.pi/180)
    sun_str_split = suncaribu.toLight().split(' ')

    print("\t azimut: %0.3f \t zenith: %0.3f" % 
                (-suncaribu.azim*180/math.pi, 90-(suncaribu.elev*180/math.pi)))
    print("\t x: %0.3f \t y: %0.3f \t z: %0.3f" % 
                                                    (float(sun_str_split[1]),
                                                    float(sun_str_split[2]),
                                                    float(sun_str_split[3])))
    print("\n")

