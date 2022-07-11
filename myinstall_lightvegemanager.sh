#!/bin/bash

# installe Wheat-Fspm pour avoir les données d'entrée
git clone --recurse-submodules https://github.com/openalea-incubator/WheatFspm.git
cd WheatFspm
git submodule update --remote

# conversion retour chariot windows to unix
find . -type f -print0 | xargs -0 dos2unix