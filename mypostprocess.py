from matplotlib import markers
import pandas as pd
import matplotlib.pyplot as plt


## Validation de CARIBU dans lightvegemanager
fspm = pd.read_csv("WheatFspm/mau/one_wheat_resuts.csv")
lvm = pd.read_csv("outputs/dynamic_cn_wheat_csv/wheat_resuts_caribu.csv")
'''
diff=[]
iter=[]
for ite in fspm["Iteration"]:
    for shape in fspm[fspm.Iteration==ite]["Shapes"]:
        par_true = fspm[(fspm.Iteration==ite) & (fspm.Shapes==shape)]["PARa"].values[0]
        par_mine = lvm[(lvm.Iteration==ite) & (lvm.Shapes==shape)]["PARa"].values[0]
        diff.append(par_true-par_mine)
        iter.append(ite)

fi1=plt.figure()
plt.scatter(iter, diff, marker='+')
plt.title("Implémentation de CARIBU dans lightvegemanager")
plt.xlabel("Itération temps")
plt.ylabel("Différence du PAR\n (cnwheat - lightvegemanager(caribu)) en µmol.m-2.s-1")
fi1.show()
'''

## comparaison RATP et CARIBU
caribu = pd.read_csv("outputs/dynamic_cn_wheat_csv/dynamic_cnwheat_lvmcaribu.csv")
ratp1cm = pd.read_csv("outputs/dynamic_cn_wheat_csv/dynamic_cnwheat_lvmratp_1cm.csv")
ratp5mm = pd.read_csv("outputs/dynamic_cn_wheat_csv/dynamic_cnwheat_lvmratp_5mm.csv")

iter=[]
iter_tige=[]
iter_feuille=[]
carpar_feuilles=[]
carpar_tiges=[]
rapar1cm_feuilles=[]
rapar5mm_feuilles=[]
rapar1cm_tiges=[]
rapar5mm_tiges=[]
ratio_1cm_tige=[]
ratio_5mm_feuille=[]
ratio_1cm_feuille=[]
ratio_5mm_tige=[]
for ite in caribu["Iteration"]:
    if ite < 600:
        sum_carib_tige=0
        sum_carib_feui=0
        sum_ra1_tige=0
        sum_ra1_feui=0
        sum_ra5_tige=0
        sum_ra5_feui=0
        for shape in caribu[caribu.Iteration==ite]["Shapes"]:
            if shape == 19 or shape == 34:
                sum_carib_tige += caribu[(caribu.Iteration==ite) & (caribu.Shapes==shape)]["PARa"].values[0]
            else:
                sum_carib_feui += caribu[(caribu.Iteration==ite) & (caribu.Shapes==shape)]["PARa"].values[0]

        for shape in ratp1cm[ratp1cm.Iteration==ite]["Shapes"]:
            if shape == 19 or shape == 34:
                sum_ra1_tige += ratp1cm[(ratp1cm.Iteration==ite) & (ratp1cm.Shapes==shape)]["PARa"].values[0]
            else:
                sum_ra1_feui += ratp1cm[(ratp1cm.Iteration==ite) & (ratp1cm.Shapes==shape)]["PARa"].values[0]
        
        for shape in ratp5mm[ratp5mm.Iteration==ite]["Shapes"]:
            if shape == 19 or shape == 34:
                sum_ra5_tige += ratp5mm[(ratp5mm.Iteration==ite) & (ratp5mm.Shapes==shape)]["PARa"].values[0]
            else:
                sum_ra5_feui += ratp5mm[(ratp5mm.Iteration==ite) & (ratp5mm.Shapes==shape)]["PARa"].values[0]
    
        iter.append(ite)
        carpar_feuilles.append(sum_carib_feui)
        carpar_tiges.append(sum_carib_tige)
        rapar1cm_feuilles.append(sum_ra1_feui)
        rapar5mm_feuilles.append(sum_ra5_feui)
        rapar1cm_tiges.append(sum_ra1_tige)
        rapar5mm_tiges.append(sum_ra5_tige)
    
    for shape in caribu[caribu.Iteration==ite]["Shapes"]:
        car = caribu[(caribu.Iteration==ite) & (caribu.Shapes==shape)]["PARa"].values[0]
        if shape in ratp1cm[(ratp1cm.Iteration==ite)]["Shapes"].values:
            ra1 = ratp1cm[(ratp1cm.Iteration==ite) & (ratp1cm.Shapes==shape)]["PARa"].values[0]
            ra5 = ratp5mm[(ratp5mm.Iteration==ite) & (ratp5mm.Shapes==shape)]["PARa"].values[0]
            if car!=0:
                if shape == 19 or shape == 34:
                    ratio_1cm_tige.append(ra1/car)
                    ratio_5mm_tige.append(ra5/car)
                    iter_tige.append(ite)
                else:
                    ratio_1cm_feuille.append(ra1/car)
                    ratio_5mm_feuille.append(ra5/car)
                    iter_feuille.append(ite)

'''
fig, ((ax1), (ax2)) = plt.subplots(2, 1)
fig.subplots_adjust(hspace=0.5)
ax1.plot(iter, carpar_feuilles,label='CARIBU')
ax1.plot(iter, rapar1cm_feuilles,label='RATP vox 1cm')
ax1.plot(iter, rapar5mm_feuilles,label='RATP vox 5mm')
ax1.set_title("PAR absorbé cumulé sur les feuilles")
ax1.set_xlabel("Itération temps")
ax1.set_ylabel("PARa cumulé en µmol.m-2.s-1")
ax1.legend()

ax2.plot(iter, carpar_tiges,label='CARIBU')
ax2.plot(iter, rapar1cm_tiges,label='RATP vox 1cm')
ax2.plot(iter, rapar5mm_tiges,label='RATP vox 5mm')
ax2.set_title("PAR absorbé cumulé sur les tiges")
ax2.set_xlabel("Itération temps")
ax2.set_ylabel("PARa cumulé en µmol.m-2.s-1")
ax2.legend()
'''

fig, ((ax1), (ax2)) = plt.subplots(2, 1)
fig.subplots_adjust(hspace=0.4)
ax1.scatter(iter_tige, ratio_1cm_tige,label='RATP vox 1cm', marker="+")
ax1.scatter(iter_tige, ratio_5mm_tige,label='RATP vox 5mm', marker="+")
ax1.plot(iter_tige,[sum(ratio_1cm_tige)/len(ratio_1cm_tige)]*len(iter_tige),label='moyenne RATP vox 1cm = %.3f'%(sum(ratio_1cm_tige)/len(ratio_1cm_tige)))
ax1.plot(iter_tige,[sum(ratio_5mm_tige)/len(ratio_5mm_tige)]*len(iter_tige),label='moyenne RATP vox 5mm = %.3f'%(sum(ratio_5mm_tige)/len(ratio_5mm_tige)))
ax1.set_title("Comparaison des PAR absorbés sur les tiges")
ax1.set_xlabel("Itération temps")
ax1.set_ylabel("Ratio ratp.PARa/caribu.PARa")
ax1.legend()

ax2.scatter(iter_feuille, ratio_1cm_feuille,label='RATP vox 1cm', marker="+")
ax2.scatter(iter_feuille, ratio_5mm_feuille,label='RATP vox 5mm', marker="+")
ax2.plot(iter_feuille,[sum(ratio_1cm_feuille)/len(ratio_1cm_feuille)]*len(iter_feuille),label='moyenne RATP vox 1cm = %.3f'%(sum(ratio_1cm_feuille)/len(ratio_1cm_feuille)))
ax2.plot(iter_feuille,[sum(ratio_5mm_feuille)/len(ratio_5mm_feuille)]*len(iter_feuille),label='moyenne RATP vox 5mm = %.3f'%(sum(ratio_5mm_feuille)/len(ratio_5mm_feuille)))
ax2.set_title("Comparaison des PAR absorbés sur les feuilles")
ax2.set_xlabel("Itération temps")
ax2.set_ylabel("Ratio ratp.PARa/caribu.PARa")
ax2.legend()

plt.show()


print("tige")
print("moyenne %f %f" % (sum(ratio_1cm_tige)/len(ratio_1cm_tige),sum(ratio_5mm_tige)/len(ratio_5mm_tige)))
print("max %f %f"%(max(ratio_1cm_tige),max(ratio_5mm_tige)))
print("min %f %f"%(min(ratio_1cm_tige),min(ratio_5mm_tige)))
print("\n")
print("feuille")
print("moyenne %f %f"%(sum(ratio_1cm_feuille)/len(ratio_1cm_feuille),sum(ratio_5mm_feuille)/len(ratio_5mm_feuille)))
print("max %f %f"%(max(ratio_1cm_feuille),max(ratio_5mm_feuille)))
print("min %f %f"%(min(ratio_1cm_feuille),min(ratio_5mm_feuille)))