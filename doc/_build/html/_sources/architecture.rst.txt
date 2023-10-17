.. _architecture:

Architecture of the tool
=========================

.. _package:

Package Content
-----------------

The package has 8 folders and 4 files in its root:

- data: data files used in some the use examples in the package
 
- doc: user and reference documentations
  
- examples: small python scripts with examples of LightVegeManager use

- notebooks: jupyter notebooks with documented tutorials exploring the tool features

- s2v and s5: sources of analysis tools
  
- src: sources of LightVegeManager
  
- tests: unit testing files written in pytest format
  
- .gitignore: files and folders to ignore by git

- .readthedocs.yaml: config file for read-the-docs

- requirements.txt: packages needed for read-the-docs in order to compile
  
- README.md: Read Me file

- pytest.ini: config file for pytest
  
- setup.py: setup script used for installation

.. _structure:

Code structure
---------------

The code structure of the tool is designed with a main class LightVegeManager, which calls modules for computing specific 
parts depending on the inputs.

.. figure:: _img/diagram_classmodules.png
    :align: center
    
    Dependances between modules


Modules are:

    * :mod:`LightVegeManager_tool` Main class of the tool. Calls all the other modules in ``src``.

    * :mod:`LightVegeManager_CARIBUinputs` Manages and prepares input information for CARIBU.

    * :mod:`LightVegeManager_RATPinputs` Manage vegetation and meteo input informations for RATP
    
    * :mod:`LightVegeManager_RiRi5inputs` Manage grid of voxels for RiRi5

    * :mod:`LightVegeManager_VTK` Writes VTK files from LightVegeManager geometry and lighting data. Used for visualisation

    * :mod:`LightVegeManager_basicgeometry` Provides basic geometric operations in 3D used by LightVegeManager. 

    * :mod:`LightVegeManager_buildRATPscene` Build a PyRATP.grid from inputs.

    * :mod:`LightVegeManager_defaultvalues` Default for simulation fixed parameters

    * :mod:`LightVegeManager_leafangles` Handles leaf angle distribution, both in its dynamic computing or reading a file

    * :mod:`LightVegeManager_outputs` Manages and reformats output results from the light models in pandas.Dataframe with similar columns names

    * :mod:`LightVegeManager_plantGL` Visualisation with plantGL

    * :mod:`LightVegeManager_sky` Build a sky from LightVegeManager input informations

    * :mod:`LightVegeManager_stems` Manages stems element

    * :mod:`LightVegeManager_sun` Build a sun respecting each light model format

    * :mod:`LightVegeManager_tesselator` Manages subdivision of a triangulation

    * :mod:`LightVegeManager_transfer` Manages transfer of LightVegeManager results to plant Models

    * :mod:`LightVegeManager_trianglesmesh` Builds and handles triangulation mesh.

    * :mod:`LightVegeManager_voxelsmesh` Builds and handles axis oriented voxels mesh


.. _frontend:

Front End: the main commands 
------------------------------

th tool is used through the class ``LightVegeManager`` and the following methods:

    * constructor ``__init__`` builds the sky, which stays the same throughout all the simulation. It sets also default values if not precised in the inputs.

    * :meth:`build` creates a common geometric scene from inputs and set parameters for the light model.

    * :meth:`run` computes the lighting. 

The outputs from radiations are automatically gathered in a pandas Dataframe and accessible from the getters :meth:`elements_outputs`, :meth:`triangles_outputs` and :meth:`voxels_outputs`.

As part of our initial objective, we added two methods in order to convert the results in formats understandable by CN-Wheat and l-egume:

    * :meth:`to_MTG`, which updates a MTG table read by CN-Wheat
    
    * :meth:`to_l_egume`, which updates two tables read by l-egume

.. note:: l-egume needs a local information of transmitted lighting among a voxel grid. Then, you need to provide the grid dimensions to LightVegeManager.

The other getters available are:

    * :meth:`sensors_outputs`, outputs of virtual sensors, only with CARIBU
    * :meth:`soilenergy`, radiation received by the soil, only with CARIBU
    * :meth:`sun`, object containing sun position xyz
    * :meth:`maxtrianglearea`, if you entered a triangle mesh, return the largest triangle
    * :meth:`tesselationtime`, if you activated tesselation of a triangle mesh (redraw of triangles), return computation time
    * :meth:`modelruntime`, return the computation time of the light model

Finally, you also have additional tools available for analysing the inputs and visualising the outputs (:ref:`additional tools <other>`).

.. _backend:

Back End: More details about how the tool works
-------------------------------------------------

First of all, the geometric merging is set in a Caribu scene format. It is a dict where keys are indices and values are list of triangles, one triangle is list of 3 3-tuple representing the vertices. 
Here is an example:

.. code-block:: python

    # organ 1
    triangle1 = [(0,0,0), (0,1,0), (1,1,0)]
    triangle2 = [(0,0,0), (0,1,1), (1,1,1)]

    # organ 2
    triangle3 = [(0,0,2), (0,1,2), (1,1,2)]

    caribuscene = { 1 : [triangle1, triangle2] ,
                    2 : [triangle3] }

We choose this format for its low processing cost, because it uses basic python objects.

We also save and recreate a dict to organize the indices inside each input scene called ``matching_ids`` (:ref:`index managment <indexes>`). 

The input parameters defines a common way for setting each light model parameters.

Geometric merging
*****************

.. figure:: _img/merging.png
    :align: center

    Tool workflow for geometry

We built a module to tesselate one triangle into four smaller triangles. The tesselation is applied either uniformly among all the triangle mesh or on a sides of a voxels grid.
The tesselation following a grid is made in order to have a better matching of triangles in a voxels grid and attenuate the error of converting the mesh type.

We also built the possibility to compute dynamicaly the leaf angle distribution, either globally among all triangles, or locally inside each voxel.

.. _indexes:

Managing indexes
*****************

LightVegeManager expects a list of geometric scenes in its inputs. Each scene represents a plant specy. 
Each scene can also by sorted by elements, which can represent plant organs or just sets of triangles.
The tool will then reorder all the indexes in order to avoid any confusion. 
The correspondance between input indices and new indices is stored in the dict attribute ``matching_ids``.

Example of reordering:

.. figure:: _img/indices.png
    :align: center

    Reordering indices in LightVegeManager


Triangle subdivision
*********************

We implemented an algorithm for subdividing a triangle in 4 triangles according to triangle position in a grid of voxels.
The goal was to have a better matching between a triangulation and a grid of voxels. 
Triangles are subdivided if they are between several voxels. 

.. figure:: _img/tesselation.PNG
    :align: center
    :scale: 50%

    Subdivision of a triangle


CARIBU, l-egume and virtual sensors
************************************

l-egume needs two different values to understand lighting results:

    * total intercepted radiation for each plant
    
    * local transmitted radiation following a voxel grid

In order to use CARIBU with l-egume, you need to retrieve transmitted radiations for each position of a voxel grid. 
LightVegeManager implements functions which can create a set of virtual sensors following a voxel grid.
Then, with the virtual sensors and the soilmesh options activated, you can calculate the local transmitted radiation.
Make sure to have the same grid dimensions as l-egume intern grid.


