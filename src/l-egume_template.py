import numpy as np
import os

from openalea.lpy import *
import legume
import IOxls
import IOtable
import run_legume_usm as runl
import ShootMorpho as sh
import daily_loop as loop
import RIRI5 as riri

def initialisation_legume(pathin, foldin, foldout, fusms, ongletBatch):
    
    #lecture de la liste des usm
    mn_path = os.path.join(pathin, foldin, fusms)
    usms = IOxls.xlrd.open_workbook(mn_path)
    ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))

    # cree la liste de L-systems et liste des noms
    lsystem_simulations = {}
    names_simulations = []
    for i in range(len(ls_usms['ID_usm'])):
        if int(ls_usms['torun'][i]) == 1:  # si 1 dans la colonne 'torun' l'ajoute a la liste
            #mylsys = runl.lsystemInputOutput_usm(path_, fxls, i, foldin=foldin, ongletBatch=ongletBatch)
            mylsys = runl.lsystemInputOutput_usm(fusms, foldin=foldin, ongletBatch=ongletBatch, i=i, path_OUT=foldout)
            name = list(mylsys)[0]
            names_simulations.append(name)
            lsystem_simulations[name] = mylsys[name]

    return lsystem_simulations, names_simulations



def iteration_legume_withoutlighting():
    return