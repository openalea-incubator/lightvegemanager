'''
gestion des entrées des modèles de lumières

'''
from alinea.caribu.sky_tools import spitters_horaire

from PyRATP.pyratp.vegetation import Vegetation
from PyRATP.pyratp.micrometeo import MicroMeteo

def RATP_vegetation(parameters, angle_distrib, reflected) :
    entities_param = []
    if parameters["angle distrib algo"] != "compute voxel":
        for id, mu_ent in enumerate(parameters["mu"]):
            if reflected : 
                reflectance_coef = parameters["reflectance coefficients"][id]
            else :
                reflectance_coef = [0., 0.]

            entities_param.append({
                                    'mu' : mu_ent,
                                    'distinc' : angle_distrib["global"][id],
                                    'rf' : reflectance_coef
                                    })

        return Vegetation.initialise(entities_param)
    else :
        for id, mu_ent in enumerate(parameters["mu"]):
            if reflected : 
                reflectance_coef = parameters["reflectance coefficients"][id]
            else :
                reflectance_coef = [0., 0.]
            entities_param.append({
                                    'mu' : mu_ent,
                                    'rf' : reflectance_coef
                                    })

        return Vegetation.initialise(entities_param, 
                                        pervoxel=True, 
                                        distribvox=angle_distrib["voxel"])

def RATP_meteo(energy,    
                day, 
                hour, 
                coordinates,
                parunit, 
                truesolartime,
                direct,
                diffus) :

    if parunit == "micromol.m-2.s-1" :
        #: Spitters's model estimating for the diffuse:direct ratio
        # coefficient 2.02 : 4.6 (conversion en W.m-2) x 0.439 (PAR -> global)
        RdRs = spitters_horaire.RdRsH(Rg=energy/2.02, 
                                        DOY=day, 
                                        heureTU=hour, 
                                        latitude=coordinates[0])
        
        # coeff 4.6 : https://www.researchgate.net/post/Can-I-convert-PAR-photo-active-radiation-value-of-micro-mole-M2-S-to-Solar-radiation-in-Watt-m2
        energy = energy / 4.6  # W.m-2
    else:
        RdRs = spitters_horaire.RdRsH(Rg=energy/0.439, 
                                        DOY=day, 
                                        heureTU=hour, 
                                        latitude=coordinates[0])
    
    # PAR et Dif en W.m^-2
    if diffus :
        # Direct et diffus
        if direct: 
            rdif = energy*RdRs
        # uniquement diffus
        else:
            rdif = energy
    # uniquement direct
    else :  rdif = 0.
        
    return MicroMeteo.initialise(doy=day, 
                                    hour=hour, 
                                    Rglob=energy, 
                                    Rdif=rdif, 
                                    truesolartime=truesolartime)