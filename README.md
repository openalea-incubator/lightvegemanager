# Gestion de la lumière entre des modèles de plantes

## Notes sur le fonctionnement
Le package PyRATP est appelé en local et fait parti du package de lightvegemanager. Son développement est en parallèle de l'outil de gestion.

## Installation de l'outil
1) Création d'un environnement conda avec miniconda3
    ```bash
    conda create -n myenvname python=3.7 \
                    openalea.mtg=2.0.5 openalea.plantgl=3.14.1 openalea.lpy=3.9.2 alinea.caribu=8.0.7 alinea.astk=2.1.0 xlrd=2.0.1\
                     coverage=6.3 nose=1.3.7 sphinx=4.4.0 statsmodels=0.13.1 numpy=1.20.3 scipy=1.7.3 pandas=1.3.4 progressbar2=3.37.1\
                      -c conda-forge -c fredboudon
    ```

2) Se placer dans l'environnement créé : `conda activate myenvname`

3) Installation de Adel-Wheat (requis par WheatFspm)
   - Linux
        ```bash
        wget https://github.com/rbarillot/adel/archive/python3.zip
        unzip python3.zip && cd adel-python3
        python setup.py develop
        ```

    - Windows
        - télécharger et déziper l'archive : https://github.com/rbarillot/adel/archive/python3.zip
        - lancer la commande suivante dans le dossier `adel-python3`
            ```bash
            python setup.py develop
            ```
4) Installation de Wheat-fspm
    1) dans une consol Git :
        ```bash
        git clone --recurse-submodules https://github.com/openalea-incubator/WheatFspm.git
        cd WheatFspm
        git submodule update --remote
        ```
    2) installation dans l'environnement conda (dans le dossier `WheatFspm`)
        ```bash
        python setup.py develop
        ```

5) Installation de l-egume
    1) dans une consol Git :
        ```bash
        git clone -b Develop https://github.com/glouarn/l-egume.git
        ```
    2) installation dans l'environnement conda (dans le dossier `l-egume`)
        ```bash
        python setup.py develop
        ```

6) installation de LightVegeManager
    1) dans une consol Git :
        ```bash
        git clone https://github.com/mwoussen/lightvegemanager
        ```
    2) installation dans l'environnement conda (dans le dossier `lightvegemanager`)
        ```bash
        python setup.py develop
        ```

7) Installation de PyRATP en local, compilation avec numpy
   Dans le dossier `lightvegemanager/PyRATP/f90`
   - Linux
        ```bash
        bash compile_pyratp_lin64.sh
        ```
    - Windows
        ```bash
        .\compile_pyratp_win64.bat
        ```
        /!\ attention : les chemins relatifs de la ligne 5 (appel de `gfortran.exe`) peuvent être erronés, il faudra les modifier manuellement.

## Utilisation sur le meso@LR avec un conteneur singularity comme environnement d'exécution
* avoir un pyratp.so compilé par le conteneur pour être utilisable
* les scripts python à exécuter doivent avoir la permission d'exécution : `chmod +x script.py`
* tous les fichiers compatibles pour unix 

### compilation de la partie fortran de PyRATP
Il y a un script pour linux 64 bits et Windows 64 bits pour compiler la partie fortran avec f2py.
Ces scripts sont à lancer dans un environnement avec numpy d'installé (normalement installé dans le conteneur), avant de run les simulations.

1) d'abord rentrer dans le conteneur
    ```bash
    singularity shell mycontainer.sif
    ```

2) puis lancer la compilation fortran -> python 
    ```bash
    cd PyRATP/f90
    bash compile_pyratp_lin64.sh
    ```

3) vérifier que la nouvelle librairie a bien été installée dans PyRATP/pyratp sous le nom pyratp.so
    ```bash
    ls -lt PyRATP/pyratp/*.so
    ```

### mise à jour du dépôt 
1) on revient au dernier commit (nécessaire car dos2unix a modifié tous les fichiers et entraine des conflits)
    ```bash
    git reset --hard HEAD
    ```

2) on télécharge le nouveau commit du dépot
    ```bash
    git pull origin master
    ```

3) on reconvertit les fichiers au format unix
    ```bash
    # lance la derniere ligne du fichier
    tail -1 myinstall_lightvegemanager.sh | bash
    ```
    ou
    ```bash
    find . -type f -print0 | xargs -0 dos2unix
    ```

4) Regarde quels fichiers ont été modifiés au dernier commit
    ```bash
    git diff --name-only HEAD HEAD~1
    ```

5) Si la partie fortran de RATP a été modifiée on la recompile avec le conteneur et f2py