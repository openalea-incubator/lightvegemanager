.. _inputs:

Inputs
======

In LightVegeManager, we choose to organize the input parameters in 3 python dict: geometric, environment and light model parameters. 

.. _geometric:

Geometric inputs
----------------

.. code-block:: python
    
    geometry = {
			   "scenes" : [scene0, scene1, scene2, ...] ,
			   "domain" : ((xmin, ymin), (xmax, ymax)),
			   "stems id" : [(id_element, id_scene), ...],
			   "transformations" : {
									"scenes unit" : kwarg ,
									"rescale" : kwarg ,
									"translate" : kwarg ,
									"xyz orientation" : kwarg
									}
                }

It contains the geometric scenes and directives to operate on them.
You can also add descriptions on the scenes, especially if you have a large numbers of scenes and plant model evolved.
Each scene must contains only one specy of plant, except for l-egume grid format, which can contains multiple species. 
As the geometric informations are likely to change during a simulation, this dict is not read by the constructor but rather by the :meth:`build` method.

.. _scenes:

Scenes
*******

``"scenes"`` key needs a list of the geometric scenes. 
Each scene represents one plant specie, except for voxel grid format. There is no restrictions on the numbers of input scenes. 
Possible format are:
    
    * PlantGL scene
    
    * MTG table
    
    * VGX file
    
    * CaribuScene dict
         
        .. code-block:: python
     
            # a vertice is a 3-tuple of float
            triangle1 = [ (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) ]
            triangle2 = [ (x4, y4, z4), (x5, y5, z5), (x6, y6, z6) ]
            mycaribuscene = { 888 : [ triangle1, triangle2 ] }

            # 888 is the id of the element composed by triangle1 and triangle2


    * Voxels grid: a dict with two keys:
        
        * ``"LA"``: a ``scipy.array`` of dimension ``[number of species, nz, ny, nx]`` and stores for each voxel of each specie a leaf area
        
        * ``"distrib"``: a list of ``scipy.array`` for each specy storing a leaf angle distribution

``"domain"`` is necessary for CARIBU and infinite scene, it defines the xy square to replicate. 
``"stems id"`` indicates the organ id of possible stems. They will be considered as strictly opaque element.

Geometric transformations
**************************

There is an optional fourth input key in the geometric dict, used to specify geometric transformations on several scenes.
For example, you can arrange your inputs scenes in a specific pattern by translating them.
The value of each entry is a dict: 

.. code-block:: python
    
    "transformation" : { index_from_input_scenes : transformation_value}

Possible keys are the following:

``"scenes unit"``
#################

Associates a measure unit among this list ``['mm','cm','dm','m','dam','hm','km']`` with a scene. 
You can also assign a unit to the inside global scene, the tool will rescale all the scenes in the same unit.

example:

.. code-block:: python

    transformations["scenes unit"] = {0 : 'm', 
                                      1 : 'cm', 
                                      2 : 'hm'}

``"rescale"``
#############

Multiply all the lengths of a scene by a ratio.

example:

.. code-block:: python

    transformations["rescale"] = {0 : 0.5, 
                                  2 : 5}

``"translate"``
###############

Apply a translation vector to a scene. A vector is 3-tuple.

example:

.. code-block:: python

    transformations["translate"] = {0 : (5, 0.2, 1), 
                                    1 : (0, 5, 0)}

``"xyz orientation"``
######################

Change the orientation of a scene among the xyz axis. 
When a model creates a scene in a different convention used by the light models, LightVegeManager allows for a common convention.
RATP and CARIBU use the convention x+ matches to North.

Possible convention:

    * ``"x+ = S"``: x+ matches to South

    * ``"x+ = W"``: x+ matches to West

    * ``"x+ = E"``: x+ matches to East

example:

.. code-block:: python

    transformations["xyz orientation"] = {0 : "x+ = S", 
                                          1 : "x+ = W"}

.. _environment:

Environment parameters
----------------------

.. code-block:: python

    environment = {
                    "coordinates" : [latitude, longitude, timezone] ,
                    
                    "sky" : "turtle46" ,
                    "sky" : ["file", filepath] ,
                    "sky" : [nb_azimut, nb_zenith, "soc" or "uoc"] ,

                    "direct" : bool, # sun radiations
                    "diffus" : bool, # sky radiations
                    "reflected" : bool, # reflected radiation in the canopy
                    "infinite" : bool, # infinitisation of the scene
                    }


Contains parameters related to the environment of modelling and static parameters during all the simulation.

The sky is defined by the ``"sky"`` key and has 3 different possibilities:
  
    -  ``"turtle46"`` : SOC sky with 46 directions in "turtle" shape
  
    -  ``["file", filepath]`` : sky defined by the file ``filepath``
  
    -  ``[nb_azimut, nb_zenith, "soc" ou "uoc"]``: computes dynamically the source directions 


Light model parameters
----------------------

Specific parameters for each light model, as they use very different syntax and approach.

.. _caribuparam:

CARIBU: surfacic approach
**************************

.. code-block:: python

    caribu_args = {
                        "sun algo" : "ratp",
                        "sun algo" : "caribu",

                        "caribu opt" : {
                                        band0 = (reflectance, transmittance),
                                        band1 = (reflectance, transmittance),
                                        ...
                                        },
                        "debug" : bool,
                        "soil mesh" : bool,
                        "sensors" : ["grid", dxyz, nxyz, orig, vtkpath, "vtk"]
                    }


The key ``"sun algo"`` specifies how to calculate the sun position. Possibilities are ``"ratp"`` or ``"caribu"``, which will call the function related to either one of the light models.

``"caribu opt"`` defines the optical parameters in CARIBU format. For example:
        
    .. code-block:: python

        caribu_parameters["caribu opt"] = {
                                            "par" = (0.10, 0.07) ,
                                            "nir" = (0.15, 0.05)
                                            } 

For virtual sensors, we propose the construction following a grid of voxels, where a sensor is the bottom side of a voxel.
        
        .. code-block:: python
            
            caribu_parameters["sensors"] = ["grid", dxyz, nxyz, orig, path, "vtk"]
        
        Where:

        *  ``dxyz = [x, y, z]`` : list of size of each side of a voxel
        
        *  ``nxyz = [nx, ny, nz]`` : list of number of voxels in each direction
        
        *  ``path`` : string of file path where to save the vtk file of the sensors
        
        *  ``"vtk"`` : flag to specify you want to write sensors in VTK format for visualisation
        
        ``path`` and ``"vtk"`` are optional elements in the list
    
    .. note:: At the moment, only grid option is available for sensors
    

.. _ratpparam:

RATP: volumic approach
***********************

.. code-block:: python

    ratp_args = {
                    # Grid specifications
                    "voxel size" : [dx, dy, dz],
                    "voxel size" : "dynamic",
                    
                    "origin" : [xorigin, yorigin, zorigin],
                    "origin" : [xorigin, yorigin],

                    "number voxels" : [nx, ny, nz],
                    "grid slicing" : "ground = 0."
                    "tesselation level" : int

                    # Leaf angle distribution
                    "angle distrib algo" : "compute global",
                    "angle distrib algo" : "compute voxel",
                    "angle distrib algo" : "file",

                    "nb angle classes" : int,
                    "angle distrib file" : filepath,

                    # Vegetation type
                    "soil reflectance" : [reflectance_band0, reflectance_band1, ...],
                    "reflectance coefficients" : [reflectance_band0, reflectance_band1, ...],
                    "mu" : [mu_scene0, mu_scene1, ...]
                }


Grid specifications
###################

* ``"voxel size"``: defines the size of one voxel, 2 possibilities
    
    * ``[float, float, float]``: for a length on each direction (xyz)
    
    * ``"dynamic"``: create squared voxels where the length of a side is l = 5.max(area(triangle))

* ``"number voxels"``: ``[nx, ny, nz]``, defines the number of voxels on each xyz axis. If not specified, the tool computes dynamically the grid size according to input geometries.

* ``"origin"``: possibility to define the grid origin.
    
    * ``[xorigin, yorigin, zorigin]``
    
    * ``[xorigin, yorigin]`` and ``zorigin = -zmin``

* ``"grid slicing"``: if ``"grid slicing" : "ground = 0."`` and ``"number voxels"`` is not specified and there are triangles in the input geometries, it rescales the grid to avoid the soil layer < 0

* ``"tesselation level"``: if there are triangles in the input geometries, you can tesselate them a certain number of times to match the voxel grid.

Leaf angle distribution
########################

You can specify how the leaf angle distribution is computed. It can be computed dynamically from a triangulation.

* ``"angle distrib algo"``: 3 possibilities 
    
    * ``"compute global"``: LightVegeManager computes a global leaf angle distribution
    
    * ``"compute voxel"``: LightVegeManager computes a local leaf angle distribution in each voxel of the grid
    
    * ``file``: LightVegeManager reads a global leaf angle distribution in a file

* ``"nb angle classes"``: if ``"angle distrib algo" = "compute global"`` or ``"angle distrib algo" = "compute voxel"``, number of angle classes in the distribution

* ``"angle distrib file"``: if ``"angle distrib algo" = file"`` string precising the file path to read

Vegetation type
###############

* ``"soil reflectance"``: if precised, soil reflects rays by a coefficient, one float by wavelength ``[coef_par, coef_nir]``

* ``"reflectance coefficients"``: if precised, leaf reflects rays by a coefficient, one float by wavelength ``[coef_par, coef_nir]``

* ``"mu"``: clumping effect for each input species