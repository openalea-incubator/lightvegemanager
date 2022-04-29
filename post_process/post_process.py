from matplotlib import markers
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

## comparaison RATP et CARIBU
all = pd.read_csv("outputs/dynamic_cn_wheat_csv/saved/dynamic_cnwheat_dense_2.csv")

iter_tige=[[],[]]
tige_detail=[[],[]]
iter_detail = [[],[],[],[],[]]
feuilles_detail = [[],[],[],[],[]]
diff_tiges=[]

memite=-1
for ite in all["Iteration"]:    
    if ite != memite:
        for shape in all[all.Iteration==ite]["Shapes"]:
            car = all[(all.Iteration==ite) & (all.Shapes==shape)]["PARi CARIBU"].values[0]
            ra1 = all[(all.Iteration==ite) & (all.Shapes==shape)]["PARa RATP"].values[0]
            val = 100*abs(car-ra1)/all[(all.Iteration==ite) & (all.Shapes==shape)]["PAR input"].values[0]
            if car!=0 :# and car < all[(all.Iteration==ite) & (all.Shapes==shape)]["PAR input"].values[0]:
                if shape == 19 :
                    tige_detail[0].append(val)
                    iter_tige[0].append(ite)
                elif shape==34:
                    tige_detail[1].append(val)
                    iter_tige[1].append(ite)
                else:
                    if shape==813:
                        iter_detail[0].append(ite)
                        feuilles_detail[0].append(val)
                        
                    elif shape==814:
                        iter_detail[1].append(ite)
                        feuilles_detail[1].append(val)

                    elif shape==51:
                        iter_detail[2].append(ite)
                        feuilles_detail[2].append(val)

                    elif shape==815:
                        iter_detail[3].append(ite)
                        feuilles_detail[3].append(val)

                    elif shape==816:
                        iter_detail[4].append(ite)
                        feuilles_detail[4].append(val)
    memite=ite

'''
print("STATS TIGES :")
print(stats.describe(np.array(diff_tiges)))
print(np.median(np.array(diff_tiges)))

print("\n-----------------------------------\n")
print("STATS FEUILLES :")
for i, s in enumerate([813, 814, 51, 815, 816]):
    print("\t SHAPE ",s)
    print(stats.describe(np.array(feuilles_detail[i])))
    print(np.median(np.array(feuilles_detail[i])))
'''
fig, ((ax1), (ax2)) = plt.subplots(2, 1)
fig.subplots_adjust(hspace=0.5)
ax1.plot(iter_tige[0], tige_detail[0],label='shape 19 (tube)')
ax1.plot(iter_tige[1], tige_detail[1],label='shape 34 (capuchon)')

ax2.plot(iter_detail[0], feuilles_detail[0],label='shape 813')
ax2.plot(iter_detail[1], feuilles_detail[1],label='shape 814')
ax2.plot(iter_detail[2], feuilles_detail[2],label='shape 51')
ax2.plot(iter_detail[3], feuilles_detail[3],label='shape 815')
ax2.plot(iter_detail[4], feuilles_detail[4],label='shape 816',marker='+')

#ax1.set_title("Différence sur les tiges")
ax1.set_title("Erreur relative sur les tiges")
ax1.set_xlabel("Itération temps (# ligne fichier météo)")
#ax1.set_ylabel("Différence CARIBU-PAR (µmol.m-2.s-1)")
ax1.set_ylabel("Erreur relative à CARIBU (%)")
ax1.legend()

#ax2.set_title("Différence sur les feuilles")
ax2.set_title("Erreur relative sur les feuilles")
ax2.set_xlabel("Itération temps (# ligne fichier météo)")
#ax2.set_ylabel("Différence CARIBU-PAR (µmol.m-2.s-1)")
ax2.set_ylabel("Erreur relative à CARIBU (%)")
ax2.legend()

plt.show()