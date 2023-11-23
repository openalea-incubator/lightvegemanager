# Light management between plant models

![](doc/_img/illus_lightvegemanager.png)

**Authors** : Maurane Woussen, Romain Barillot, Didier Combes, Gaëtan Louarn

**Institutes** : INRAE

**Status** : Python package 

**License** : [Cecill-C](https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)

**URL** : https://github.com/mwoussen/lightvegemanager

**Documentation** : https://lightvegemanager.readthedocs.io/en/latest/

## Overview

LightVegeManager is a Python package made for plant modeling and managing enlightment. It serves as an interface to merge several plant models and a light model. Applications involved by this tools are part of OpenAlea platform.

For pratical uses of this tool with mixed crops, see [Plant Fusion](https://github.com/mwoussen/plantfusion)

## Required dependencies

- python 3
- numpy
- pandas

## Geometric dependencies

- [Openalea PlantGL](https://github.com/openalea/plantgl)
- [Openalea MTG](https://github.com/openalea/mtg)

## Light Model dependencies

- [CARIBU](https://github.com/openalea-incubator/caribu)
- [PyRATP (fork)](https://github.com/mwoussen/PyRATP)
- [RiRi5](https://github.com/glouarn/riri5)

## Installing the tool

### Targeted installation

1) Create a conda environment with miniconda3 and required dependencies
    ```bash
        conda create -n myenvname python=3 numpy=1.20.3 pandas -c conda-forge
    ```

2) Place yourself in your conda environment
    ```bash
        conda activate myenvname
    ```

2) Then install your preferred packages

    PlantGL
    ```bash
        conda install openalea.plantgl -c openalea3
    ```
    
    MTG
    ```bash
        conda install openalea.mtg -c openalea3
    ```

    CARIBU
    ```bash
        conda install alinea.caribu alinea.astk -c openalea3
    ```

    PyRATP
    1) Git console :
        ```bash
        git clone -b update_mobidiv https://github.com/mwoussen/PyRATP
        ```
    2) installation in the conda environment (in folder `PyRATP`)
        ```bash
        make mode=develop
        make clean
        ```
    RiRi5
    1) Git console :
        ```bash
        git clone https://github.com/glouarn/riri5
        ```
    2) installation in the conda environment (in folder `riri5`)
        ```bash
        python setup.py develop
        ```
    
    Development tools
    ```bash
        conda install pytest sphinx sphinx-rtd-theme -c conda-forge
        pip install nbsphinx pgljupyter
    ```
    
### Complete installation

1) Create a conda environment with miniconda3
    ```bash
        conda create -n myenvname python=3 openalea.mtg openalea.plantgl alinea.caribu alinea.astk numpy=1.20.3 pandas pytest sphinx sphinx-rtd-theme -c conda-forge -c openalea3
    ```

2) Place yourself in your conda environment
    ```bash
        conda activate myenvname
    ```

3) Install PyRATP
    1) Git console :
        ```bash
        git clone -b update_mobidiv https://github.com/mwoussen/PyRATP
        ```
    2) installation in the conda environment (in folder `PyRATP`)
        ```bash
        make mode=develop
        make clean
        ```

4) Install RiRi5
    1) Git console :
        ```bash
        git clone https://github.com/glouarn/riri5
        ```
    2) installation in the conda environment (in folder `riri5`)
        ```bash
        python setup.py develop
        ```

5) Install LightVegeManager
    1) Git console :
        ```bash
        git clone https://github.com/mwoussen/lightvegemanager
        ```
    2) installation in the conda environment (in folder `lightvegemanager`)
        ```bash
        python setup.py develop
        ```

## Contributing

The documentation is made with `sphinx` and unit testing with `pytest`.

## Licence

This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
