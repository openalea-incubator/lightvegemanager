import pandas as pd
import matplotlib.pyplot as plt

fspm = pd.read_csv("WheatFspm/mau/one_wheat_resuts.csv")
lvm = pd.read_csv("outputs/dynamic_cnwheat/my_wheat_resuts.csv")

diff=[]
iter=[]
for ite in fspm["Iteration"]:
    for shape in fspm[fspm.Iteration==ite]["Shapes"]:
        par_true = fspm[(fspm.Iteration==ite) & (fspm.Shapes==shape)]["PARa"].values[0]
        par_mine = lvm[(lvm.Iteration==ite) & (lvm.Shapes==shape)]["PARa"].values[0]
        diff.append(par_true-par_mine)
        iter.append(ite)

plt.scatter(iter, diff)
plt.show()
