'''
    LightVegeManager
    ****************

    Main class of the tool. Calls all the other modules in ``src``.

    3 inputs dict for setting all parameters:

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

    Currently LightVegeManager handles the light models RATP and CARIBU:

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

    .. seealso:: For more details :ref:`Inputs description <inputs>`
'''

from alinea.pyratp.runratp import runRATP
from alinea.caribu.CaribuScene import CaribuScene
# import riri5

import time
import os
import subprocess

from lightvegemanager.outputs import *
from lightvegemanager.trianglesmesh import *
from lightvegemanager.voxelsmesh import *
from lightvegemanager.stems import *
from lightvegemanager.leafangles import *
from lightvegemanager.tesselator import *
from lightvegemanager.sun import *
from lightvegemanager.sky import *
from lightvegemanager.RATPinputs import *
from lightvegemanager.CARIBUinputs import *
from lightvegemanager.buildRATPscene import *
from lightvegemanager.transfer import *
from lightvegemanager.VTK import *
from lightvegemanager.plantGL import *
from lightvegemanager.defaultvalues import *

class LightVegeManager(object) :
    """Main class for the tool LightVegeManager
    
    Common simulation order: 
    
    ``input geometries -> build and prepare data -> call light model -> transfer results to plant models``

    It includes:

        Main methods:

        * ``__init__``: initializes and builds static object for the rest of simulation
        * :meth:`build`: builds and prepare all geometric meshes
        * :meth:`run`: calls a light model and manages its inputs and outputs

        Transfer methods:

        * :meth:`to_MTG`: transfers ligthing results to a MTG table
        * :meth:`to_l_egume`: transfers ligthing results to l-egume by creating two arrays as inputs for the plant model

        Analysis tools: analyses a set of triangles and formats them as a turbid medium inputs

        * :meth:`s5`: fortran tool
        * :meth:`s2v`: c++ tool

        Visualisation tools:

        * :meth:`VTK_nolight`: write VTK file with only geometric informations
        * :meth:`VTK_light`: write VTK file with geometric informations and associated light results
        * :meth:`VTK_sun`: : write VTK file representing the sun as a line

        Getters: to use light results with external routines

        * :meth:`legume_transmitted_light`
        * :meth:`legume_intercepted_light`
        * :meth:`elements_outputs`
        * :meth:`triangles_outputs`
        * :meth:`voxels_outputs`
        * :meth:`sensors_outputs`
        * :meth:`sun`
        * :meth:`soilenergy`
        * :meth:`maxtrianglearea`
        * :meth:`legume_empty_layers`
        * :meth:`tesselationtime`
        * :meth:`modelruntime`

    :param environment: Environment parameters, defaults to {}
    :type environment: dict, optional
    :param lightmodel: either ``"ratp"`` or ``"caribu"``, defaults to ""
    :type lightmodel: str, optional
    :param lightmodel_parameters: light model parameters, defaults to {}
    :type lightmodel_parameters: dict, optional
    :param main_unit: measure unit for the global scene where the light will be computed, defaults to "m"
    :type main_unit: str, optional
    :raises ValueError: lightmodel entry not valid, either ``'ratp'`` or ``'caribu'``
    """     
    ## MAIN METHODS ##
    def __init__(self,
                    environment={},
                    lightmodel="",
                    lightmodel_parameters={},
                    main_unit="m") :        
        # gets default variables from LightVegeManager_defaultvalues.py
        default_environnement,    \
        default_ratp_parameters,  \
        default_caribu_parameters = default_LightVegeManager_inputs()

        # check if choosen light model is known by the tool
        if lightmodel != "ratp" and lightmodel != "caribu" :
            raise ValueError("Unknown lightmodel: can be either 'ratp' or 'caribu' ")

        # save inputs in the instance
        self.__environment = environment
        self.__lightmodel = lightmodel
        self.__lightmodel_parameters = lightmodel_parameters
        self.__main_unit = main_unit

        # copy of default values in input parameters on not initialized keys
        for key, value in default_environnement.items() :
            if key not in self.__environment : 
                self.__environment[key] = value
        
        if lightmodel == "caribu" : lghtdict = default_caribu_parameters
        elif lightmodel == "ratp" : lghtdict = default_ratp_parameters
        for key, value in lghtdict.items() :
            if key not in self.__lightmodel_parameters : 
                self.__lightmodel_parameters[key] = value

        # sky building
        skytype = self.__environment["sky"]
        if lightmodel == "caribu":
            if self.__environment["diffus"] : self.__sky = CARIBUsky(skytype)
        elif lightmodel == "ratp" : self.__sky = RATPsky(skytype)     

    def build(self, geometry = {}, global_scene_tesselate_level=0) :
        """Builds a mesh of the simulation scene in the right light model format

        :param geometry: geometric parameters, contains geometric scenes, defaults to {}
        :type geometry: dict, optional
        :param global_scene_tesselate_level: option to subdivide all triangles of the mesh a certain number of times (to fine tuning the mesh), defaults to 0
        :type global_scene_tesselate_level: int, optional
        :raises ValueError: Currently, converting voxels mesh to triangles mesh is not possible
        """
        self.__geometry = geometry

        # First process of the scenes list, it gathers all triangulations
        self.__complete_trimesh, \
        self.__matching_ids,      \
        legume_grid,               \
        id_legume_scene = chain_triangulations(self.__geometry["scenes"])

        # applies geometric transformations on inputs scenes if precised
        if "transformations" in self.__geometry :
            apply_transformations(self.__complete_trimesh, 
                                    self.__matching_ids, 
                                    self.__geometry["transformations"], 
                                    self.__main_unit)

        self.__areamax = compute_area_max(self.__complete_trimesh)

        # global tesselation of triangulation
        if self.__matching_ids and global_scene_tesselate_level > 0 :
            new_trimesh = {}
            for id_ele, triangles in self.__complete_trimesh.items() :
                new_tr_scene = []
                for t in triangles :
                    level = 0
                    iterate_triangles(t, level, global_scene_tesselate_level, 
                                                                new_tr_scene)
                new_trimesh[id_ele] = new_tr_scene
            self.__complete_trimesh = new_trimesh

        self.__pmin, \
        self.__pmax = compute_minmax_coord(self.__complete_trimesh)

        self.__triangleLmax = compute_trilenght_max(self.__complete_trimesh)

        if self.__lightmodel == "caribu":
            if legume_grid :
                raise ValueError("Conversion from voxels grid to triangles \
                                is not possible yet")
        
        # Builds voxels grid from input geometry
        elif self.__lightmodel == "ratp":
            # triangles in the inputs
            if self.__matching_ids :
                # separates stem elements in a new specy
                if "stems id" in self.__geometry :
                    manage_stems_for_ratp(self.__geometry["stems id"], 
                                            self.__matching_ids, 
                                            self.__lightmodel_parameters)
                else : self.__geometry["stems id"] = None

                arg = (self.__complete_trimesh,        
                        [self.__pmin, self.__pmax],     
                        self.__triangleLmax,             
                        self.__matching_ids,              
                        self.__lightmodel_parameters,      
                        self.__environment["coordinates"], 
                        self.__environment["reflected"], 
                        self.__environment["infinite"],      
                        self.__geometry["stems id"],          
                        len(self.__geometry["scenes"]))
    
                self.__complete_voxmesh, \
                self.__matching_tri_vox, \
                self.__angle_distrib =  build_RATPscene_from_trimesh(*arg)

            # creates an empty RATP grid of voxels if geometric inputs are empty
            else :
                arg = (self.__lightmodel_parameters,      
                        [self.__pmin, self.__pmax],    
                        self.__environment["coordinates"],  
                        self.__environment["infinite"])

                self.__complete_voxmesh, \
                self.__angle_distrib =  build_RATPscene_empty(*arg)
                self.__matching_tri_vox = {}
            
            # if there is a grid of voxels in the inputs, we converts it in a RATP grid
            if legume_grid and not self.__matching_ids :
                arg = (self.__geometry["scenes"][0],
                        self.__lightmodel_parameters,       
                        self.__environment["coordinates"],  
                        self.__environment["reflected"], 
                        self.__environment["infinite"])

                self.__complete_voxmesh, \
                self.__angle_distrib, \
                self.__nb0 =  legumescene_to_RATPscene(*arg)
                    
    def run(self, 
            energy=0., 
            day=0, 
            hour=0, 
            parunit="micromol.m-2.s-1", 
            truesolartime=False,
            id_sensors=None) :
        """Calls the light model and formats lighting results

        :param energy: input radiation energy, defaults to 0
        :type energy: float, optional
        :param day: simulation day, defaults to 0
        :type day: int, optional
        :param hour: simulation hour, defaults to 0
        :type hour: int, optional
        :param parunit: input energy unit, light models manages radiations in different unit, you can precise input unit and LightVegeManager will convert it in the right unit, defaults to "micromol.m-2.s-1"
        :type parunit: str, optional
        :param truesolartime: simulation hour is a true solar time or local time (depending on simulation coordinates), defaults to False
        :type truesolartime: bool, optional
        :param id_sensors: if you use CARIBU with a grid of virtual sensors, you have to precise which input scenes the grid must match, defaults to None
        :type id_sensors: list, optional
        :raises ValueError: with CARIBU you can precise the sun algorithm to calculate sun position, can be either ``"caribu"`` or ``"ratp"``
        :raises ValueError: valid radiations are ``"direct"``, ``"diffuse"``, ``"reflected"``
        """            
        self.__energy = energy
        
        ## RATP ##
        if self.__lightmodel == "ratp" :
            vegetation = RATP_vegetation(self.__lightmodel_parameters, 
                                            self.__angle_distrib, 
                                            self.__environment["reflected"])
            meteo = RATP_meteo(energy,    
                                day, 
                                hour, 
                                self.__environment["coordinates"],
                                parunit, 
                                truesolartime,
                                self.__environment["direct"],
                                self.__environment["diffus"])

            if self.__complete_voxmesh.nveg > 0 :
                # Run of RATP
                print("debut ratp")
                start=time.time()
                res = runRATP.DoIrradiation(self.__complete_voxmesh, 
                                            vegetation, 
                                            self.__sky, 
                                            meteo)
                self.__time_runmodel = time.time() - start
                print("fin ratp")
                # output management
                self.__voxels_outputs = out_ratp_voxels(
                                                    self.__complete_voxmesh,
                                                    res,
                                                    parunit)

                # if there are triangulations in the inputs
                if self.__matching_ids :
                    arg = (self.__complete_trimesh,
                                self.__matching_ids,
                                self.__matching_tri_vox,
                                self.__voxels_outputs)
                    self.__triangles_outputs = out_ratp_triangles(*arg)

                    arg = (self.__matching_ids,
                    self.__environment["reflected"],
                    self.__lightmodel_parameters["reflectance coefficients"],
                    self.__triangles_outputs)
                    self.__elements_outputs = out_ratp_elements(*arg)

            else :
                print("--- Empty RATP grid")
                self.__voxels_outputs = out_ratp_empty_grid(day, hour)

        ## CARIBU ##
        elif self.__lightmodel == "caribu" :
            sun_up = False
            if self.__environment["direct"] :
                # computes sun position
                arg = (day, 
                        hour, 
                        self.__environment["coordinates"], 
                        truesolartime)
                if self.__lightmodel_parameters["sun algo"]=="ratp":
                    self.__sun = ratp_sun(*arg)
                elif self.__lightmodel_parameters["sun algo"]=="caribu":
                    self.__sun = caribu_sun(*arg)
                else :
                    raise ValueError("sun algo not recognize")
            
                # check if sun is up
                # criteria is set to elevation > 2°, from alinea shortwave_balance.F90 -> DirectBeam_Interception, l.68
                sun_up = math.asin(-self.__sun[2])*(180/math.pi) > 2.
            
            compute = (sun_up and self.__environment["direct"]) or \
                    self.__environment["diffus"]
            if compute:             
                # CARIBU preparations
                arg = (self.__complete_trimesh, 
                        self.__geometry, 
                        self.__matching_ids, 
                        (self.__pmin, self.__pmax), 
                        self.__lightmodel_parameters,
                        self.__environment["infinite"], 
                        id_sensors)
                opt, sensors_caribu, debug = Prepare_CARIBU(*arg)
                issensors = "sensors" in self.__lightmodel_parameters
                issoilmesh = self.__lightmodel_parameters["soil mesh"] != -1
                
                # Initialize a CaribuScene
                if self.__environment["diffus"]  and \
                    self.__environment["direct"] :
                    c_scene_sky = CaribuScene(
                        scene=self.__complete_trimesh, 
                        light=self.__sky, 
                        opt=opt, 
                        scene_unit=self.__main_unit, 
                        pattern=self.__geometry["domain"], 
                        soil_mesh = self.__lightmodel_parameters["soil mesh"], 
                        debug = debug
                        )
                    sun_sky_option = "mix"
                    light = [tuple((1., self.__sun))]

                else :
                    if self.__environment["diffus"] :
                        light = self.__sky
                        sun_sky_option = "sky"
                    
                    elif self.__environment["direct"] :
                        light = [tuple((1., self.__sun))]
                        sun_sky_option = "sun"

                    else :
                        raise ValueError("Error with radiative inputs")

                c_scene = CaribuScene(
                    scene=self.__complete_trimesh, 
                    light=light, 
                    opt=opt, 
                    scene_unit=self.__main_unit, 
                    pattern=self.__geometry["domain"], 
                    soil_mesh = self.__lightmodel_parameters["soil mesh"], 
                    debug = debug
                    )
                
                # Runs CARIBU
                arg = [c_scene, \
                        not self.__environment["reflected"], \
                        self.__environment["infinite"], \
                        sensors_caribu]
                if sun_sky_option == "mix":
                    start = time.time()
                    raw_sun, aggregated_sun = run_caribu(*arg)
                    arg[0] = c_scene_sky
                    raw_sky, aggregated_sky = run_caribu(*arg)
                    self.__time_runmodel = time.time() - start

                    #: Spitters's model estimating for the diffuse:direct ratio
                    # % de sky dans la valeur d'énergie finale
                    Rg = energy / 2.02  #: Global Radiation (W.m-2)
                    #: Diffuse fraction of the global irradiance
                    rdrs = spitters_horaire.RdRsH(
                                Rg=Rg, 
                                DOY=day, 
                                heureTU=hour, 
                                latitude=self.__environment["coordinates"][0])

                    raw,                    \
                    aggregated,             \
                    self.__sensors_outputs, \
                    self.__soilenergy = out_caribu_mix(rdrs,
                                                        c_scene,
                                                        c_scene_sky,
                                                        raw_sun, 
                                                        aggregated_sun, 
                                                        raw_sky, 
                                                        aggregated_sky,
                                                        issensors,
                                                        issoilmesh)
                else :
                    start = time.time()
                    raw, aggregated = run_caribu(*arg)
                    self.__time_runmodel = time.time() - start

                    self.__sensors_outputs, \
                    self.__soilenergy = out_caribu_nomix(c_scene,
                                                            aggregated,
                                                            issensors,
                                                            issoilmesh)
                    
            # Outputs management
            arg = [day,
                hour,
                self.__complete_trimesh,
                self.__matching_ids,
                aggregated,
                compute]
            self.__elements_outputs = out_caribu_elements(*arg)
            arg[4] = raw
            self.__triangles_outputs = out_caribu_triangles(*arg)

            if "sensors" in self.__lightmodel_parameters and self.__lightmodel_parameters["sensors"][-1] == "vtk":
                # create list with radiative value for each sensor triangle (2 triangles per sensor)
                var=[]
                for id, triangles in sensors_caribu.items():
                    for t in triangles :
                        var.append(self.__sensors_outputs["par"][id])

                sensor_path = self.__lightmodel_parameters["sensors"][-2] + "sensors_h"+str(int(hour))+"_d"+str(int(day)) + ".vtk"
                VTKtriangles(sensors_caribu, [var], ["intercepted"], sensor_path)
            
        ## RiRi (l-egume) ##
        elif self.__lightmodel == "riri" :
            dS = self.__lightmodel_parameters["voxel size"][0] * \
                                self.__lightmodel_parameters["voxel size"][1]
            
            tag_light_inputs = [self.__geometry["scenes"][0]["LA"] / dS, 
                                self.__geometry["scenes"][0]["triplets"], 
                                self.__geometry["scenes"][0]["distrib"], 
                                energy * dS]

            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            start = time.time()
            self.__legume_transmitted_light, \
            self.__legume_intercepted_light = \
                        riri.calc_extinc_allray_multi_reduced(
                                        *tag_light_inputs, 
                                        optsky=self.__in_environment["sky"][0], 
                                        opt=self.__in_environment["sky"][1]
                                                                )
            self.__time_runmodel = time.time() - start
    
    ## TRANSFER METHODS ##
    def to_MTG(self, energy=1., mtg=None, id=None) :
        """Transfers lighting results to a MTG table.
        
        .. warning:: The :meth:`run` must have been called before to have results dataframes.
        
        Results are ``pandas.Dataframe`` stored in the ``self``

        :param energy: input energy, defaults to 1
        :type energy: float, optional
        :param mtg: MTG table with ``"PARa"`` and ``"Erel"`` entries in its properties, defaults to None
        :type mtg: MTG, optional
        :param id: you can precise to which input scenes the MTG table corresponds, defaults to None
        :type id: list or tuple, optional
        :raises AttributeError: you need to call :func:run first
        """        
        if not hasattr(self, '_LightVegeManager__elements_outputs') :
            raise AttributeError("No results yet, run a light modeling first")
        
        # crée un tableau comme dans caribu_facade de fspm-wheat
        dico_par = {}
        para_dic = {}
        erel_dic = {}
        if id is None :
            for s in self.__elements_outputs["Organ"]:
                d = self.__elements_outputs[self.__elements_outputs.Organ==s]
                para_dic[s] = d["par Eabs"].values[0] * energy
                erel_dic[s] = d["par Eabs"].values[0] / self.__energy
        
        elif type(id) == list or type(id) == tuple :
            for esp in id:
                df_outputs_esp = self.__elements_outputs[ \
                                self.__elements_outputs.VegetationType==esp]
                for s in df_outputs_esp["Organ"]:
                    d = df_outputs_esp[df_outputs_esp.Organ==s]
                    para_dic[s] = d["par Eabs"].values[0] * energy
                    erel_dic[s] = d["par Eabs"].values[0] / self.__energy
        dico_par["PARa"] = para_dic
        dico_par["Erel"] = erel_dic

        for param in dico_par:
            if param not in mtg.properties():
                mtg.add_property(param)
            # update the MTG
            mtg.property(param).update(dico_par[param])

    def to_l_egume(self, 
                    energy=1., 
                    m_lais = [], 
                    list_lstring = [], 
                    list_dicFeuilBilanR = [], 
                    list_invar = [], 
                    id=None) :
        """Transfers lighting results to l-egume

        .. warning:: The :meth:`run` must have been called before to have results dataframes.

        .. warning:: l-egume needs transmitted energy informations located in a grid of voxels. You need to have the same dimensions in the lighting results.
            
            * With RATP, RATP grid must have the same dimensions as l-egume intern grid. 
            
            * With CARIBU, you must create a grid of virtual sensors in the same dimensions as l-egume intern grid.

        Results are ``pandas.Dataframe`` stored in the ``self``

        :param energy: input energy, defaults to 1
        :type energy: float, optional
        :param m_lais: leaf area represented in a numpy.array of dimension [number of species, number of z layers, number of y layers, number of x layers], defaults to []
        :type m_lais: numpy.array, optional
        :param list_lstring: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict lstring stores the l-system of each plant, defaults to []
        :type list_lstring: list of dict, optional
        :param list_dicFeuilBilanR: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict dicFeuiBilanR stores correspondances between voxels grid and each plant, defaults to []
        :type list_dicFeuilBilanR: list of dict, optional
        :param list_invar: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict invar stores instant intern variables of l-egume., defaults to []
        :type list_invar: list of dict, optional
        :param id: list of indices from input scenes which corresponds to current l-egume instance. If you have several plantmodel among the input scenes, you need to precise which results you want to transfer to which instance of l-egume. defaults to None
        :type id: list, optional
        :raises AttributeError: you need to call :meth:`run` first
        :raises ValueError: unknown light model

        :return: 

            if light model is RATP :func:`transfer_ratp_legume` :

                * ``res_abs_i``: absorbed energy in each voxels of a grid matching the dimensions of l-egume intern grid of voxels. One value for each specy.

                * ``res_trans``: transmitted energy in each voxels of a grid matching the dimensions of l-egume intern grid of voxels

            if light model is CARIBU :func:transfer_caribu_legume :

                * update of list_invar: updates the keys ``"parap"`` and ``"parip"`` for each specy. Cumulatative energy per plant

                * ``res_trans``: transmitted energy in each voxels of a grid matching the dimensions of l-egume intern grid of voxels
        
        :rtype: numpy.array
        """                   
        if not hasattr(self, '_LightVegeManager__elements_outputs') :
            raise AttributeError("No results yet, run a light modeling first")

        epsilon = 1e-14

        if self.__lightmodel == "ratp" :
            return transfer_ratp_legume(m_lais, 
                                        energy, 
                                        self.__complete_voxmesh,
                                        self.__voxels_outputs,
                                        self.__nb0,
                                        epsilon)

        elif self.__lightmodel == "caribu" :
            skylayer = reduce_layers_from_trimesh(
                                self.__complete_trimesh, 
                                self.__pmax,
                                self.__lightmodel_parameters["sensors"][1],
                                self.__lightmodel_parameters["sensors"][2],
                                self.__matching_ids, 
                                id)
            
            return transfer_caribu_legume(energy,
                                            skylayer, 
                                            self.__elements_outputs, 
                                            self.__sensors_outputs, 
                                            self.__lightmodel_parameters["sensors"][1], 
                                            self.__lightmodel_parameters["sensors"][2], 
                                            m_lais, 
                                            list_invar, 
                                            list_lstring,
                                            list_dicFeuilBilanR, 
                                            self.__environment["infinite"],
                                            epsilon) 

        else :
            raise ValueError("Unknown light model (ratp or caribu)")

    ## EXTERN TOOLS ##
    def s5(self):
        """Creates inputs files for s5 and runs it
        s5 is an external tool made to analyse a set of triangles in order to use a grid of voxels.
        It also computes leaf angle distribution from the triangulation

        .. note:: All files are located in ``s5`` folder

        **Input files created**

            * fort.51: contains triangulation. Possibility to precise stem elements it is registered in the instance of LightVegeManager

            * s5.par: stores grid of voxels informations and number of entities

        **Output files created**

            * fort.60
            
                dimensions xy of the grid

                stats by specy

                    - total leaf area
                    - leaf area index
                    - global zenital leaf angle distribution
                    - global azimutal leaf angle distribution
                
                stats by voxels
                
                    - #specy | #ix | #iy | #iz (coordinate xyz of the voxel) | eaf area density
                    - local zenital leaf angle distribution
                    - local azimutal leaf angle distribution

            * leafarea: for each specy, for each voxel
                
                ix | iy | iz | #specy | LAD | zenital-distribution | azimutal-distribution

        **example**

        >>> myscene # a plantgl Scene
        >>> testofs5 = LightVegeManager() # create a instance
        >>> testofs5.build( geometry={ "scenes" : [myscene] } ) # build the geometry
        >>> testofs5.s5() # run of s5, creates input and output files

        """        
        currentfolder = os.path.abspath( \
                                    os.path.dirname(os.path.dirname(__file__)))
        s5folder = os.path.join(currentfolder, os.path.normpath("s5"))
        fort51 = os.path.join(s5folder, os.path.normpath("fort.51"))
        s5par = os.path.join(s5folder, os.path.normpath("s5.par"))

        # écriture du fichier fort.51 contenant la triangulation
        f=open(fort51, 'w')
        c_tr=1
        for id, triangles in self.__complete_trimesh.items() :
            for t in triangles :
                if tuple(self.__matching_ids[id]) in \
                                                self.__geometry["stems id"]:
                    stem='000'
                else:
                    stem='001'
                label = str(self.__matching_ids[id][1] + 1) +         \
                        str("%05i"%(self.__matching_ids[id][1] + 1)) + \
                        stem + '000'
                
                f.write("p\t1\t%s\t3\t"%(label))
                for i in range(3):
                    f.write("%f\t%f\t%f\t"%(t[i][0],t[i][1],t[i][2]))
                f.write("\n")
                c_tr += 1
        f.close()

        # écriture du fichier s5.par contenant les informations de la grille
        f=open(s5par, 'w')    
        f.write("%i\t%i\t%i\t%i\t0\n"%(self.__complete_voxmesh.nent, 
                                            9, 9, 
                                                self.__complete_voxmesh.njz))
        for i in range(self.__complete_voxmesh.njz):
            f.write("%f\t"%(self.__complete_voxmesh.dz[i]))
        f.write("\n")

        njx = self.__complete_voxmesh.njx
        njy = self.__complete_voxmesh.njy
        dx = self.__complete_voxmesh.dx
        dy = self.__complete_voxmesh.dy
        f.write("%f\t%i\t%f\t%f\t%i\t%f\n"%(njx*dx, njx, dx, njy*dy, njy, dy))
        f.close()

        # exécution de s5 dans un sous process
        subprocess.call(".\s5.exe", shell=True, cwd=s5folder)
        
        print("\n"+"--- Fin de s5.f")

    def s2v(self):
        """Creates inputs files for s2v and runs it
        s5 is an external tool made to analyse a set of triangles in order to use a grid of voxels.
        It also computes leaf angle distribution from the triangulation

        .. note:: All files are located in ``s2v`` folder

        **Input files created**

            * fort.51: stores triangulation. Possibility to precise stem elements it is registered in the instance of LightVegeManager

            * s2v.par: stores grid of voxels informations and number of entities

        **Output files created**

            * s2v.log

                logs about processing
                
                global statistics

                    - total leaf area per specy
                    - total leaf area
                    - leaf area index
                    - global zenital leaf angle distribution

            * s2v.can

                z layer where each triangle is located (and copy its vertices)
            
            * s2v.area

                triangle id | z layer | triangle area

            * out.dang: SAIL file

                line 1: global leaf area index for specy 1
                line 2: global zenital leaf angle distribution for specy 1

            * leafarea: SAIL file

                each line: 0 | idz | leaf area index for each slope class in current z layer | 0 | 0 | leaf area density on the layer

        **example**

        >>> myscene # a plantgl Scene
        >>> testofs2v = LightVegeManager() # create a instance
        >>> testofs2v.build( geometry={ "scenes" : [myscene] } ) # build the geometry
        >>> testofs2v.s2v() # run of s2v, creates input and output files
        
        """ 
        currentfolder = os.path.abspath( \
                                    os.path.dirname(os.path.dirname(__file__)))
        s2vfolder = os.path.join(currentfolder, os.path.normpath("s2v"))
        fort51 = os.path.join(s2vfolder, os.path.normpath("fort.51"))
        s2vpar = os.path.join(s2vfolder, os.path.normpath("s2v.par"))

        # écriture du fichier fort.51 contenant la triangulation
        f=open(fort51, 'w')
        c_tr=1
        for id, triangles in self.__complete_trimesh.items() :
            for t in triangles :
                if tuple(self.__matching_ids[id]) in \
                                                self.__geometry["stems id"]:
                    stem='000'
                else:
                    stem='001'
                label = str(self.__matching_ids[id][1] + 1) +         \
                        str("%05i"%(self.__matching_ids[id][1] + 1)) + \
                        stem + '000'
                
                f.write("p\t1\t%s\t3\t"%(label))
                for i in range(3):
                    f.write("%f\t%f\t%f\t"%(t[i][0],t[i][1],t[i][2]))
                f.write("\n")
                c_tr += 1
        f.close()

        f=open(s2vpar, 'w')
        # ligne 1 
        f.write("9 9 %i\n"%(self.__complete_voxmesh.njz))

        # ligne 2
        for i in range(self.__complete_voxmesh.njz):
            f.write("%f "%(self.__complete_voxmesh.dz[i]))
        f.write("\n")
        
        # ligne 3
        njx = self.__complete_voxmesh.njx
        njy = self.__complete_voxmesh.njy
        dx = self.__complete_voxmesh.dx
        dy = self.__complete_voxmesh.dy
        f.write("%f %i %f %f %i %f %i\n"%(njx*dx, 
                                            njx, 
                                            dx, 
                                            njy*dy, 
                                            njy, 
                                            dy, 
                                            self.__complete_voxmesh.nent))

        # ligne 4
        f.write("%f %f %f\n"%(self.__complete_voxmesh.xorig, 
                                self.__complete_voxmesh.yorig, 
                                    self.__complete_voxmesh.zorig))

        f.close()

        # exécution de s5 dans un sous process
        subprocess.call(".\s2v.exe", shell=True, cwd=s2vfolder)

        print("--- Fin de s2v.cpp")

    ## VTK FILES ##
    def VTK_nolight(self, path, i=None, printtriangles=True, printvoxels=True) :
        """Writes a VTK from mesh(es) in ``self``, with only geometric informations

        :param path: file and path name
        :type path: string
        :param i: associate the created file with an indice in its filename, defaults to None
        :type i: int, optional
        :param printtriangles: write triangulation if one has been created in :func:build, defaults to True
        :type printtriangles: bool, optional
        :param printvoxels: write grid of voxels if one has been created in :func:build, defaults to True
        :type printvoxels: bool, optional
        """             
        if self.__lightmodel == "ratp" and printvoxels: 
            if i is None : filepath = path + "_voxels_nolight.vtk"
            else : filepath = path + "_voxels_nolight" + "_" + str(i) + ".vtk"
            ratp_prepareVTK(self.__complete_voxmesh, filepath)
            
        if self.__matching_ids and printtriangles:
            if i is None : filepath = path + "_triangles_nolight.vtk"
            else : filepath = path + "_triangles_nolight" + "_" + str(i)+".vtk"
            VTKtriangles(self.__complete_trimesh, [], [], filepath)

    def VTK_light(self, path, i=None, printtriangles=True, printvoxels=True) :
        """Writes a VTK from mesh(es) in ``self``, with geometric informations and lighting results

        .. warning:: The :meth:`run` must have been called before to have results dataframes.

        :param path: file and path name
        :type path: string
        :param i: associate the created file with an indice in its filename, defaults to None
        :type i: int, optional
        :param printtriangles: write triangulation if one has been created in :meth:`build`, defaults to True
        :type printtriangles: bool, optional
        :param printvoxels: write grid of voxels if one has been created in :meth:`build`, defaults to True
        :type printvoxels: bool, optional
        :raises AttributeError: you need to call :meth:`run` first
        """
        if not hasattr(self, '_LightVegeManager__elements_outputs') :
            raise AttributeError("No results yet, run a light modeling first")

        if (self.__lightmodel == "ratp" and \
                    (not hasattr(self, '_LightVegeManager__voxels_outputs'))) :
            print("--- VTK:  No light data, run the simulation")
            self.VTK_nolight(path, i, printtriangles, printvoxels)
        elif (self.__lightmodel == "caribu" and \
                    (not hasattr(self, '_LightVegeManager__triangles_outputs'))) :
            print("--- VTK:  No light data, run the simulation")
            self.VTK_nolight(path, i, printtriangles, printvoxels)
        else :
            if self.__lightmodel == "ratp" and printvoxels: 
                if i is None : filepath = path + "_voxels_light.vtk"
                else : filepath = path + "_voxels_light" + "_" + str(i) + ".vtk"

                datanames = ['ShadedPAR', 
                                'SunlitPAR', 
                                'ShadedArea',
                                'SunlitArea',
                                'PARa', 
                                'Intercepted', 
                                'Transmitted']

                ratp_prepareVTK(self.__complete_voxmesh, 
                                    filepath, 
                                    datanames, 
                                    self.__voxels_outputs)
                
            if self.__matching_ids and printtriangles:
                if i is None : filepath = path + "_triangles_light.vtk"
                else : filepath = path + "_triangles_light" + "_" + str(i) + ".vtk"

                datanames = []
                if self.__lightmodel == "ratp" :
                    datanames = ['ShadedPAR', 
                                'SunlitPAR', 
                                'ShadedArea',
                                'SunlitArea',
                                'PARa', 
                                'Intercepted', 
                                'Transmitted']
                elif self.__lightmodel == "caribu" :
                    datanames = []
                    for band in self.__lightmodel_parameters["caribu opt"].keys() :
                        datanames.append(band + " Eabs")
                        datanames.append(band + " Ei")

                data = [list(self.__triangles_outputs[name]) for name in datanames]
                VTKtriangles(self.__complete_trimesh, data, datanames, filepath)

    def VTK_sun(self, path, scale=2, orig=(0,0,0), center=True, i=None) :
        """Write a VTK file representing the sun by a simple line

        .. warning:: The :meth:`run` must have been called before to have results dataframes.    

        if center is False, orig is the starting point the line and it ends at orig + scale*sun.position

        if center is True, orig is the middle point of the line.

        :param path: file and path name
        :type path: string
        :param scale: size of the line, defaults to 2
        :type scale: int, optional
        :param orig: starting or middle point of the line, defaults to (0,0,0)
        :type orig: tuple, optional
        :param center: if True orig is the middle of the line, otherwise it is the starting point, defaults to True
        :type center: bool, optional
        :param i: associate the created file with an indice in its filename, defaults to None
        :type i: int, optional
        :raises AttributeError: you need to call :meth:`run` first
        """        
        if not hasattr(self, '_LightVegeManager__sun') : 
            raise AttributeError("No results yet, run a light modeling first")

        if i is None : filepath = path + "_sun.vtk"
        else : filepath = path + "_sun" + "_" + str(i) + ".vtk"

        if center :
            start = tuple([a + b * scale for a,b in zip(orig, self.__sun)])
            end = tuple([a + b * -scale for a,b in zip(orig, self.__sun)])
        else :
            start = orig
            end = tuple([a + b * scale for a,b in zip(orig, self.__sun)])

        VTKline(start, end, filepath)
    
    ## GETTERS ##
    @property
    def legume_transmitted_light(self):
        """transmitted results if light model is RiRi light

        :return: transmitted energy following a grid of voxels
        :rtype: numpy.array
        """        
        return self.__legume_transmitted_light 
    
    @property
    def legume_intercepted_light(self):
        """Intercepted results if light model is RiRi light, 

        :return: intercepted energy following a grid of voxels, for each specy
        :rtype: numpy.array
        """        
        return self.__legume_intercepted_light 

    @property
    def elements_outputs(self):
        """Lighting results aggregate by element

        :return: Column names can change depending on the light model. Commonly, there is element indice, its area and intercepted energy
        :rtype: pandas.Dataframe
        """        
        return self.__elements_outputs 

    @property
    def triangles_outputs(self):
        """Lighting results aggregate by triangle, if it has a triangulation in its inputs

        :return: .. seealso:: :mod:` LightVegeManager_outputs` for column names
        :rtype: pandas.Dataframe
        """        
        return self.__triangles_outputs 
    
    @property
    def voxels_outputs(self):
        """Lighting results aggregate by voxels, only with RATP as the selected light model

        :return: .. seealso:: :mod:` LightVegeManager_outputs` for column names
        :rtype: pandas.Dataframe
        """         
        return self.__voxels_outputs 

    @property
    def sensors_outputs(self):
        """Lighting results aggregate by sensors.
        Only with CARIBU if you activated the virtual sensors option

        :return: The output format is the same as CARIBU output format. Dict with ``Eabs`` key for absorbed energy and ``Ei`` for incoming energy
        :rtype: dict
        """        
        try:
            return self.__sensors_outputs
        except AttributeError:
            return None

    @property
    def sun(self):
        """Return sun position of the last :func:run call

        :return: vector (x, y, z)
        :rtype: tuple
        """        
        return self.__sun

    @property
    def soilenergy(self):
        """Return soil energy, only with CARIBU if soilmesh option is activated

        :return: The output format is the same as CARIBU output format. Dict with ``Eabs`` key for absorbed energy and ``Ei`` for incoming energy
        :rtype: dict
        """        
        return self.__soilenergy

    @property
    def maxtrianglearea(self):
        """Returns the largest triangle of triangles mesh
        Computed in :meth:`build`

        :return: area the triangle
        :rtype: float
        """        
        try :
            return self.__areamax
        except AttributeError:
            return 0.
  
    @property
    def legume_empty_layers(self):
        """Returns number of empty layers between top of the canopy and number of z layers expected by l-egume

        :return: result of :func:`fill_ratpgrid_from_legumescene`
        :rtype: int
        """        
        return self.__nb0

    @property
    def tesselationtime(self):
        """_summary_

        :return: _description_
        :rtype: _type_
        """        
        try:
            return self.__tess_time
        except AttributeError:
            return 0.
    
    @property
    def modelruntime(self):
        """Returns running time of the light model

        :return: time in s
        :rtype: float
        """        
        try:
            return self.__time_runmodel
        except AttributeError :
            return 0.
