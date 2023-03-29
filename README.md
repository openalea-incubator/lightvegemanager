# Light management between plant models

## Dependencies
- python 3
- numpy
- pandas
- Openalea MTG
- Openalea Lpy
- Openalea PlantGL
- CARIBU
- PyRATP_Mobidiv (https://github.com/mwoussen/PyRATP_MobiDiv)

## Installing the tool
1) Create a conda environment with miniconda3
    ```bash
    conda create -n myenvname python=3.7 openalea.mtg=2.0.5 openalea.plantgl=3.14.1 openalea.lpy=3.9.2 alinea.caribu=8.0.7 alinea.astk=2.1.0 xlrd=2.0.1 coverage=6.3 nose=1.3.7 sphinx=4.4.0 statsmodels=0.13.1 numpy=1.20.3 scipy=1.7.3 pandas=1.3.4 progressbar2=3.37.1 -c conda-forge -c fredboudon
    ```

2) Place yourself in the created environment  : `conda activate myenvname`

3) Temporary step: add modifications to CARIBU in order to use virtual sensors and soilmesh
   1) download the file https://raw.githubusercontent.com/mwoussen/caribu/master/src/alinea/caribu/CaribuScene.py and https://raw.githubusercontent.com/mwoussen/caribu/master/src/alinea/caribu/caribu_shell.py
   2) replace the two python files in your conda environment folder:
      -   in Windows: `User/Anaconda3/envs/myenvname/Lib/site-packages/alinea.caribu-*/alinea/caribu/`
      -   in Linux: `/opt/miniconda/lib/python3.7/site-packages/alinea.caribu-*/alinea/caribu/`

4) Install Adel-Wheat (required by WheatFspm)
   - Linux:
        ```bash
        wget https://github.com/rbarillot/adel/archive/python3.zip
        unzip python3.zip && cd adel-python3
        python setup.py develop
        ```

    - Windows:
        - download et unzip the archive : https://github.com/rbarillot/adel/archive/python3.zip
        - run the following command in the folder `adel-python3`
            ```bash
            python setup.py develop
            ```
5) Install Wheat-fspm
    1) Git console :
        ```bash
        git clone --recurse-submodules https://github.com/openalea-incubator/WheatFspm.git
        cd WheatFspm
        git submodule update --remote
        ```
    2) installation in the conda environment (in folder `WheatFspm`)
        ```bash
        python setup.py develop
        ```

6) Install l-egume
    1) Git console :
        ```bash
        git clone -b Develop https://github.com/glouarn/l-egume.git
        ```
    2) installation in the conda environment (in folder `l-egume`)
        ```bash
        python setup.py develop
        ```

7) Install LightVegeManager
    1) Git console :
        ```bash
        git clone https://github.com/mwoussen/lightvegemanager
        ```
    2) installation in the conda environment (in folder `lightvegemanager`)
        ```bash
        python setup.py develop
        ```

8) Install PyRATP_Mobidiv
    1) Git console :
        ```bash
        git clone https://github.com/mwoussen/PyRATP_MobiDiv
        ```
    2) installation in the conda environment (in folder `PyRATP_MobiDiv`)
        ```bash
        make mode=develop
        make clean
        ```