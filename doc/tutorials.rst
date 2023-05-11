.. _tutorials:

Tutorials
==========

.. _installation:

Installation
############

Setting a Python environment
----------------------------

We propose a way to install an openalea environment with conda. 

#. Create a conda environment with miniconda3

    .. code-block:: bash
        
        conda create -n myenvname python=3.7 `
                                    openalea.mtg=2.0.5 `
                                    openalea.plantgl=3.14.1 `
                                    openalea.lpy=3.9.2 `
                                    alinea.caribu=8.0.7 `
                                    alinea.astk=2.1.0 `
                                    xlrd=2.0.1 `
                                    coverage=6.3 `
                                    nose=1.3.7 `
                                    sphinx=4.4.0 `
                                    statsmodels=0.13.1 `
                                    numpy=1.20.3 `
                                    scipy=1.7.3 `
                                    pandas=1.3.4 `
                                    progressbar2=3.37.1 `
                                    -c conda-forge -c fredboudon

#. Place yourself in the created environment: ``conda activate myenvname``

#. Temporary step: add modifications to CARIBU in order to use virtual sensors and soilmesh
    
    #. download the file https://raw.githubusercontent.com/mwoussen/caribu/master/src/alinea/caribu/CaribuScene.py and https://raw.githubusercontent.com/mwoussen/caribu/master/src/alinea/caribu/caribu_shell.py
    
    #. replace the two python files in your conda environment folder:
        
        * Windows: ``User/Anaconda3/envs/myenvname/Lib/site-packages/alinea.caribu-*/alinea/caribu/``
        
        * Linux: ``/opt/miniconda/lib/python3.7/site-packages/alinea.caribu-*/alinea/caribu/``

.. note:: ``soilmesh`` in Caribu 8.0.10 is working correctly

Installation of Needed Packages
-------------------------------

Adel-Wheat (required by WheatFspm)
***********************************

* Windows:

    #. download et unzip the archive : https://github.com/rbarillot/adel/archive/python3.zip
    
    #. run the following command in the folder ``adel-python3``
        
        .. code-block:: bash

            python setup.py develop

* Linux:

.. code-block:: bash

    wget https://github.com/rbarillot/adel/archive/python3.zip
    unzip python3.zip && cd adel-python3
    python setup.py develop

Wheat-fspm
***********
    
#. Git console:
    
    .. code-block:: bash

        git clone --recurse-submodules https://github.com/openalea-incubator/WheatFspm.git
        cd WheatFspm
        git submodule update --remote

#. Installation in the conda environment (in folder ``WheatFspm``)
    
    .. code-block:: bash

        python setup.py develop

l-egume
*******
    
#. Git console:
    
    .. code-block:: bash

        git clone -b Develop https://github.com/glouarn/l-egume.git

#. Installation in the conda environment (in folder ``l-egume``)
    
    .. code-block:: bash

        python setup.py develop

LightVegeManager
****************
    
#. Git console:
    
    .. code-block:: bash

        git clone https://github.com/mwoussen/lightvegemanager

#. Installation in the conda environment (in folder ``lightvegemanager``)
   
    .. code-block:: bash

        python setup.py develop

PyRATP
***************
    
#. Git console:
    
    .. code-block:: bash

        git clone https://github.com/mwoussen/PyRATP

#. Installation in the conda environment (in folder ``PyRATP``)
    
    .. code-block:: bash

        make mode=develop
        make clean


.. _quickstart:

Quickstart
############

Small test to quickly test the tool. LightVegeManager needs at least a geometric information and a light model to call
between ``"ratp"``, ``"caribu"`` or ``"riri"``.

.. code-block:: python

    from LightVegeManager import *

    # one triangle as a geometric element
    # we write our triangle in a CaribuScene format
    organ_id = 001
    triangle_vertices = [(0,0,0), (1,0,0), (1,1,1)]
    triangle = {organ_id : [triangle_vertices]}
    geometry = { "scenes" : [triangle] }

    # surfacic lighting with CARIBU
    lighting = LightVegeManager(lightmodel="caribu")

    # build the scene
    lighting.build(geometry)

    # compute lighting
    energy = 500
    hour = 15   
    day = 264 # 21st september
    lighting.run(energy, hour, day)

    # output
    print(lighting.elements_outputs)

.. seealso:: For more details on default values, see :mod:`LightVegeManager\_defaultvalues`


Tutorials with jupyter
#######################

TODO