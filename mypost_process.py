import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

## comparaison RATP et CARIBU
caribu = pd.read_csv("outputs/dynamic_cn_wheat_csv/dynamic_cnwheat_dense_lvmcaribu.csv")
ratp1cm = pd.read_csv("outputs/dynamic_cn_wheat_csv/dynamic_cnwheat_dense_lvmratp_1cm.csv")

iter=[]
iter_tige=[]
iter_feuille=[]
carpar_feuilles=[]
carpar_tiges=[]
rapar1cm_feuilles=[]
rapar1cm_tiges=[]
ratio_1cm_tige=[]
ratio_1cm_feuille=[]

iter_detail = [[],[],[],[],[]]
feuilles_detail = [[],[],[],[],[]]

for ite in caribu["Iteration"]:    
    for shape in caribu[caribu.Iteration==ite]["Shapes"]:
        car = caribu[(caribu.Iteration==ite) & (caribu.Shapes==shape)]["PARa"].values[0]
        if shape in ratp1cm[(ratp1cm.Iteration==ite)]["Shapes"].values:
            ra1 = ratp1cm[(ratp1cm.Iteration==ite) & (ratp1cm.Shapes==shape)]["PARa"].values[0]
            if car!=0:
                if shape == 19 or shape == 34:
                    ratio_1cm_tige.append(ra1/car)
                    iter_tige.append(ite)
                else:
                    ratio_1cm_feuille.append(ra1/car)
                    
                    if shape==813:
                        iter_detail[0].append(ite)
                        feuilles_detail[0].append(ra1-car)
                    
                    elif shape==814:
                        iter_detail[1].append(ite)
                        feuilles_detail[1].append(ra1-car)

                    elif shape==51:
                        iter_detail[2].append(ite)
                        feuilles_detail[2].append(ra1-car)

                    elif shape==815:
                        iter_detail[3].append(ite)
                        feuilles_detail[3].append(ra1-car)

                    elif shape==816:
                        iter_detail[4].append(ite)
                        feuilles_detail[4].append(ra1-car)

for l in feuilles_detail:
    print(stats.describe(np.array(l)))
    print(np.median(np.array(l)))

plt.plot(iter_detail[0], feuilles_detail[0],label='shape 813', marker="+")
plt.plot(iter_detail[1], feuilles_detail[1],label='shape 814', marker="+")
plt.plot(iter_detail[2], feuilles_detail[2],label='shape 51', marker="+")
plt.plot(iter_detail[3], feuilles_detail[3],label='shape 815', marker="+")
plt.plot(iter_detail[4], feuilles_detail[4],label='shape 816', marker="+")
plt.xlabel("Itération temps")
plt.ylabel("Différence ratp.PARa - caribu.PARa en µmol.m-2.s-1")
plt.title("Analyse sur chaque feuille pour un couvert dense")
plt.legend()

plt.show()
