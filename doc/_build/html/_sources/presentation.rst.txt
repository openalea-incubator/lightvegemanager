.. _presentation:

Presentation
================================

Context
--------

This package was born in the framework of the research project MobiDiv (https://anr.fr/ProjetIA-20-PCPA-0006), which aims to find solutions for pesticide-free agriculture.
One of the targets of this project is to model coupled crops between wheats, modelled by CN-Wheat (https://www.quantitative-plant.org/model/CN-Wheat), and alfalfas, thanks to l-egume model (https://github.com/glouarn/l-egume).
Thus, LightVegeManager was created to manage the lighting in a vegetal simulation, in formats understandable by the tools proposed on the OpenAlea platform (https://github.com/openalea).

Intentions
-----------

We want to simplify the light management when several models, plant and lighting, are involved.
The tool manages many input and output formats, in order to match a wide variety of situations.

How
----

.. figure:: _img/how.png
    :align: center

    General Functionning

It reads the following input formats:

    * PlantGL scene
    
    * MTG table
 
    * VGX file
 
    * Caribu scene, dict storing a triangulation

    * l-egume grid, dict storing a grid of voxels

More details about the input formats here :ref:`Scenes format <scenes>`.

And it can return the outputs in the following formats:

    * updates a MTG table
 
    * returns two tables compatible with l-egume
  
    * Pandas Dataframe, in multiple scales (triangles, voxels, elements) depending on the light model.
    
Lighting
--------

At the moment, LightVegeManager can handle three lighting models, CARIBU (https://github.com/openalea-incubator/caribu), PyRATP (https://github.com/openalea-incubator/PyRATP), and RiRi modified in l-egume (https://github.com/glouarn/l-egume).
Those models are ajusted for virtual plant canopies with considerations, like infinitisation of the scene, which are preferred for simulating large crop fields.

They offers two different approach for light modelling:

Surfacic approach (CARIBU)
**************************

.. figure:: _img/triangles.PNG
    :align: center

    Example of a triangulation

The radiations are computed on a set of triangles. 
This approach offers a high precision of lighting on the different organs, but has a high cost in time computation. 

Volumic approach (PyRATP, RiRi)
*******************************

.. figure:: _img/voxels.PNG
    :align: center

    Example of a grid of voxels

The geometric scene is represented by a set of voxels, which contains a leaf area and a leaf angle informations. 
Then, the models solves turbid medium equations from the rays inside the grid.
This approach has a higher approximation of geometry, but is well adapted for dense canopy. 
The computation cost is highly reduced regarding the surfacic approach but the geometric scene must matches the turbid medium representation and hypothesis.