'''
création du ciel dans le format du lightmodel

entrée :
    - "turtle46" : ciels oc à 46 directions en forme de tortue
    - filepath : string indique un chemin
    - [nb_azimut, nb_zenith, "soc" ou "uoc"]

sortie :
    - RATP : pyratp.skyvault
    - CARIBU : [tuple((poids, tuple((dir[0], dir[1], dir[2])))), ... ]
'''
from PyRATP.pyratp.skyvault import Skyvault
from alinea.caribu.sky_tools import turtle
from alinea.caribu.sky_tools import GenSky
from alinea.caribu.sky_tools import GetLight

import math
import numpy as np



def RATPsky(skytype) :
    '''
        entrée : environnement["sky"]

        sortie : pyratp.skyvault
    '''
    input_sky = skytype
    output_sky = []

    if input_sky == "turtle46" : output_sky = Skyvault.initialise()
    elif isinstance(input_sky, str) and \
        input_sky != "turtle46" :
        output_sky = Skyvault.read(input_sky)

    elif isinstance(input_sky, list) and \
        len(input_sky) == 3 :
        ele, \
        azi,  \
        omega, \
            pc = discrete_sky(input_sky[0],  \
                                input_sky[1], \
                                input_sky[2])
        output_sky = Skyvault.initialise(ele, azi, omega, pc)
    
    else :
        raise ValueError("Unknown sky parameters : can be either \
                            'turtle46' or a string filepath or    \
                            [nb_azimut, nb_zenith, 'soc' or 'uoc'] ")
    
    return output_sky


def CARIBUsky(skytype) :
    '''
        entrée : environnement["sky"]

        sortie : [tuple((poids, tuple((dir[0], dir[1], dir[2])))), ... ]
    '''
    input_sky = skytype
    output_sky = []

    # ciel tortue à 46 directions
    if input_sky == "turtle46":
        turtle_list = turtle.turtle()
        output_sky = []
        for i,e in enumerate(turtle_list[0]):
            dir = turtle_list[2][i]
            t = tuple((e, tuple((dir[0], dir[1], dir[2]))))
            output_sky.append(t)

    # chemin d'un fichier
    elif isinstance(input_sky, str) and \
        input_sky != "turtle46" :
        
        # lecture du fichier
        listGene = []
        f = open(input_sky)
        ndir=int(f.readline().strip().split('\t')[0])
        hmoy=np.zeros(ndir)
        azmoy=np.zeros(ndir)
        omega=np.zeros(ndir)
        pc=np.zeros(ndir)
        for n in range(ndir):
            listGene.append(f.readline().strip().split('\t'))
        tabGene=np.array(listGene)
        tabGene = np.cast['float64'](tabGene)
        hmoy=np.transpose(tabGene)[0]*math.pi / 180
        azmoy=np.transpose(tabGene)[1]*math.pi / 180
        omega=np.transpose(tabGene)[2]
        pc=np.transpose(tabGene)[3]
        f.close()
        
        # convertit au format CARIBU
        output_sky = []
        for i, p in enumerate(pc) :
            dir=[0,0,0]
            dir[0] = dir[1] = math.sin(hmoy[i])
            dir[0] *= math.cos(azmoy[i])
            dir[1] *= math.sin(azmoy[i])
            dir[2] = -math.cos(hmoy[i])
            t = tuple((float(p), tuple((dir[0], dir[1], dir[2]))))
            output_sky.append(t)

    elif isinstance(input_sky, list) and \
        len(input_sky) == 3 :
        #: (Energy, soc/uoc, azimuts, zenits)
        sky_string = GetLight.GetLight(GenSky.GenSky()(1., input_sky[2], 
                                                            input_sky[0], 
                                                            input_sky[1]))  

        for string in sky_string.split('\n'):
            if len(string) != 0:
                string_split = string.split(' ')
                t = tuple((float(string_split[0]), 
                            tuple((float(string_split[1]), 
                                    float(string_split[2]), 
                                    float(string_split[3])))))
                output_sky.append(t)

    else :
        raise ValueError("Unknown sky parameters : can be either \
                            'turtle46' or a string filepath or    \
                            [nb_azimut, nb_zenith, 'soc' or 'uoc'] ")
                            
    return output_sky


def writeskyfile(h, az, omega, weights, filepath) :
    '''
    Parameters:
        - h: Elevation angle (degrees) pointing to the center of a sky fraction
        - az: Azimuth angle (degrees) pointing to the center of a sky fraction
        - omega: Solid angle associated to the direction of incidence
        - weights: Relative contribution of the sky fraction to the sky illumination)
        - filepath : name of the file to write
    '''
    return


def discrete_sky(n_azimuts, n_zeniths, sky_type) :
    ele=[]
    azi=[]
    da = 2 * math.pi / n_azimuts
    dz = math.pi / 2 / n_zeniths    
    todeg = 180/math.pi          
    for j in range(n_azimuts):
        for k in range(n_zeniths):
            azi.append((j * da + da / 2)*todeg)
            ele.append((k * dz + dz / 2)*todeg)
    n = n_azimuts*n_zeniths
    
    omega=[2*math.pi/n]*n
    
    def uoc (teta, dt, dp):
        """ teta: angle zenithal; phi: angle azimutal du soleil """
        dt /= 2.
        x = math.cos(teta-dt)
        y = math.cos(teta+dt)
        E = (x*x-y*y)*dp/2./math.pi
        return E
    def soc (teta, dt, dp):
        """ teta: angle zenithal; phi: angle azimutal du soleil """
        dt /= 2.
        x = math.cos(teta-dt)
        y = math.cos(teta+dt)
        E = (3/14.*(x*x-y*y) + 6/21.*(x*x*x-y*y*y))*dp/math.pi
        return E
    pc=[]
    for j in range(n_azimuts):
        for k in range(n_zeniths):
            azim, elv = j * da + da / 2, k * dz + dz / 2
            I=0
            if(sky_type=='soc'):
                I = soc(elv, dz, da)

    return ele, azi, omega, pc


