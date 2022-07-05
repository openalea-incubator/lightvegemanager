#!/bin/bash

# compile la partie fortran de PyRATP (pour linux 64bits)
# Sur le méso@LR : à lancer dans un conda avec numpy!
cd PyRATP/f90
bash compile_pyratp_lin64.sh

# installe Wheat-Fspm pour avoir les données d'entrée
git clone --recurse-submodules https://github.com/openalea-incubator/WheatFspm.git
cd WheatFspm
git submodule update --remote