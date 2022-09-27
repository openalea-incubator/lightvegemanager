from ftplib import error_reply
from random import uniform
import sys
import os
import getopt
import time

import numpy as np
import pandas as pd

try :
    from src.LightVegeManager import *
    from src.l_egume_template import *

except ModuleNotFoundError:
    # ajoute le dossier lightvegemanager dans le sys.path
    sys.path.insert(1, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
    from src.LightVegeManager import *
    from src.l_egume_template import *


## INITIALISATION LIGHTVEGEMANAGER
environment = {}
ratp_parameters_legume = {}
ratp_parameters_plantgl = {}

# Paramètres pré-simulation
environment["names"] = ["l-egume"]
environment["coordinates"] = [46.43,0,0] # latitude, longitude, timezone
environment["sky"] = "turtle46" # ["file", "runscripts/legume/sky_5.data"] # "turtle46" # turtle à 46 directions par défaut
environment["diffus"] = True
environment["direct"] = False
environment["reflected"] = False
environment["reflectance coefficients"] = [[0., 0.]]
environment["infinite"] = True


## PARAMETRES RATP scene grille l-egume ##
ratp_parameters_legume["voxel size"] = [0.04,0.04,0.02]
ratp_parameters_legume["soil reflectance"] = [0., 0.]
ratp_parameters_legume["mu"] = [1.]
ratp_parameters_legume["origin"] = [0., 0., 0.]

lghtratp = LightVegeManager(environment=environment,
                                lightmodel="ratp",
                                lightmodel_parameters=ratp_parameters_legume,
                                main_unit="m")

m_lais = np.zeros((1, 10, 4, 4))
for iz in range(6,10):
    for iy in range(4):
        for ix in range(4):
            m_lais[0][iz][iy][ix] = uniform(0.5, 2)

m_lais[0][7][int(uniform(1,4))][int(uniform(1,4))] = 0
m_lais[0][7][int(uniform(1,4))][int(uniform(1,4))] = 0
m_lais[0][7][int(uniform(1,4))][int(uniform(1,4))] = 0
m_lais[0][8][int(uniform(1,4))][int(uniform(1,4))] = 0
m_lais[0][8][int(uniform(1,4))][int(uniform(1,4))] = 0
m_lais[0][9][int(uniform(1,4))][int(uniform(1,4))] = 0

## Paramètres météo ## 
doy = 100
hour = 12
energy = 0.48*1*10000/(3600*24)#flux PAR journalier moyen en W.m-2 / RG en j.cm-2

# Surface foliaire et distribution des angles de l-egume
scene_legume = {}
scene_legume["LA"] = m_lais
scene_legume["distrib"] = [[0.1382,0.1664,0.1972,0.1925,0.1507,0.0903,0.0425,0.0172,0.005]]

geometry = {}
geometry["scenes"] = [scene_legume]

lghtratp.init_scenes(geometry)
lghtratp.run(PARi=energy, day=doy, hour=hour, truesolartime=True, parunit="RG")   
res_trans_2, res_abs_i_2 = lghtratp.to_l_egume(m_lais, energy)

# vérification par couche
for iz in range(m_lais.shape[1]) : 
    if np.sum(m_lais[0][iz][:][:]) > 0 :
        vox_data = lghtratp.voxels_outputs[lghtratp.voxels_outputs.Nz==(iz-lghtratp.legume_empty_layers)+1]
        print("couche %i -> area : l-egume : %.8f ratp : %.8f | sum res_trans ratp : %.8f  | sum res_abs ratp : %.8f  " % 
                                                                                                                (iz, 
                                                                                                                    np.sum(m_lais[0][iz][:][:]), 
                                                                                                                    np.sum(vox_data['Area']),
                                                                                                                    np.sum(res_trans_2[iz][:][:]), 
                                                                                                                    np.sum(res_abs_i_2[0][iz][:][:])))
print("end")