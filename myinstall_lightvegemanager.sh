#!/bin/bash

# convertit les retours chariots
dos2unix PyRATP/f90/compile_pyratp_lin64.sh

# rend exécutable les scripts du dossier runscripts
chmod +x runscripts/cnwheat/*
chmod +x runscripts/legume/*

# compile la partie fortran de PyRATP (pour linux 64bits)
# Sur le méso@LR : à lancer dans un conda avec numpy!
cd PyRATP/f90
bash compile_pyratp_lin64.sh
cd ../..

# installe Wheat-Fspm pour avoir les données d'entrée
git clone --recurse-submodules https://github.com/openalea-incubator/WheatFspm.git
cd WheatFspm
git submodule update --remote

# conversion retour chariot windows to unix
find . -type f -print0 | xargs -0 dos2unix