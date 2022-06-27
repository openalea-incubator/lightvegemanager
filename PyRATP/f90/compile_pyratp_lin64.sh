#!/bin/bash

# on est dans le dossier PyRATP/f90
# fichier encodés en UTF-8
# création de l'entete
f2py -m pyratp --debug-capi --overwrite-signature -h pyratp.pyf mod_Cocnstant_ValuesF2PY.f90 mod_Grid3DF2PY_64bit.f90 mod_SkyvaultF2PY.f90 mod_Vegetation_TypesF2PY.f90 mod_Dir_InterceptionF2PY.f90 mod_Hemi_InterceptionF2PY.f90 mod_MicrometeoF2PY.f90 mod_Shortwave_BalanceF2PY.f90 mod_Energy_BalanceF2PY.f90 mod_PhotosynthesisF2PY.f90 mod_MinerPhenoF2PY.f90 prog_RATP.f90
# compilation de la librairie
f2py -c --debug  --fcompiler=gnu95 pyratp.pyf mod_Cocnstant_ValuesF2PY.f90 mod_Grid3DF2PY_64bit.f90 mod_SkyvaultF2PY.f90 mod_Vegetation_TypesF2PY.f90 mod_Dir_InterceptionF2PY.f90 mod_Hemi_InterceptionF2PY.f90 mod_MicrometeoF2PY.f90 mod_Shortwave_BalanceF2PY.f90 mod_Energy_BalanceF2PY.f90 mod_PhotosynthesisF2PY.f90 mod_MinerPhenoF2PY.f90 prog_RATP.f90
# copie dans les source python et simplifie le nom de la librairie
cp pyratp.*.so ../pyratp/pyratp.so

