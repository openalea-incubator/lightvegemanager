'''
    
    LightVegeManager_RATPinputs
    ****************************

    Manage vegetation and meteo input informations for RATP

    The argument `parameters` refers to one the three inputs dict of LightVegeManager. It is 
    structured as so:

    .. code:: python

        ratp_args = {
                # Grid specifications
                "voxel size" : [dx, dy, dz],
                "voxel size" : "dynamic",
                
                "origin" : [xorigin, yorigin, zorigin],
                "origin" : [xorigin, yorigin],

                "number voxels" : [nx, ny, nz],
                "grid slicing" : "ground = 0."
                "tesselation level" : int

                # Leaf angle distribution
                "angle distrib algo" : "compute global",
                "angle distrib algo" : "compute voxel",
                "angle distrib algo" : "file",

                "nb angle classes" : int,
                "angle distrib file" : filepath,

                # Vegetation type
                "soil reflectance" : [reflectance_band0, reflectance_band1, ...],
                "reflectance coefficients" : [reflectance_band0, reflectance_band1, ...],
                "mu" : [mu_scene0, mu_scene1, ...]
            }

    .. seealso:: inputs.rst

'''
from alinea.caribu.sky_tools import spitters_horaire

from alinea.pyratp.vegetation import Vegetation
from alinea.pyratp.micrometeo import MicroMeteo


def RATP_vegetation(parameters, angle_distrib, reflected) :
    """Initialise a RATP Vegetation object from LightVegeManager input datas

    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param angle_distrib: leaf angle distribution
    :type angle_distrib: dict
    :param reflected: if the user wishes to activate reflected radiations
    :type reflected: bool
    :return: Vegetation types contains clumoing effect ratio, leaf angle distribution and reflectance/transmittance of leaves for each specy
    :rtype: PyRATP.pyratp.vegetation.Vegetation
    """    
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
    """Initialise a RATP MicroMeteo object from LightVegeManager input parameters

    :param energy: input ray energy
    :type energy: float
    :param day: input day
    :type day: int
    :param hour: input hour
    :type hour: int
    :param coordinates: [latitude, longitude, timezone]
    :type coordinates: list
    :param parunit: unit of energy argument, "micromol.m-2.s-1" or else
    :type parunit: string
    :param truesolartime: activates true solar time or local time to compute sun position
    :type truesolartime: bool
    :param direct: if direct rays are activated
    :type direct: bool
    :param diffus: if diffuse rays are activated
    :type diffus: bool
    :return: input meteorological data at current time step
    :rtype: PyRATP.pyratp.micrometeo.MicroMeteo
    """
    # RATP expects W.m-2 for input energy
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