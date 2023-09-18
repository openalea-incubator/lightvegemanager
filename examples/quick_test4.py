import os
from lightvegemanager.tool import LightVegeManager
from openalea.plantgl.all import Scene, Viewer

fet_fgeom = os.path.join(os.path.abspath(""), "data", "Fet-LD-F2.bgeom")
luz_fgeom = os.path.join(os.path.abspath(""), "data", "LD-F1.bgeom")
bgeom_files = [fet_fgeom, luz_fgeom]

# setup environment
environment = {}
environment["coordinates"] = [48.8 ,2.3 ,1] # latitude, longitude, timezone

# we compute only sun light in a finite scene
environment["diffus"] = False
environment["direct"] = True
environment["reflected"] = False
environment["infinite"] = False

# CARIBU parameters
caribu_parameters = {}
caribu_parameters["sun algo"] = "caribu"
caribu_parameters["caribu opt"] = {} 
caribu_parameters["caribu opt"]["par"] = (0.10, 0.05)

# inputs values for lighting
energy=500
day=264
hour=15

# scene generation parameters
nplants = 10
plant_density=130 
var_plant_position=110

from lightvegemanager.trianglesmesh import create_heterogeneous_canopy

# Initializing the tool
lighting = LightVegeManager(lightmodel="caribu", 
                                environment=environment, 
                                lightmodel_parameters=caribu_parameters)    

# generate random canopy from plant examples
if not isinstance(bgeom_files, list): bgeom_files = [bgeom_files]
scenes = []
for f in bgeom_files :
    plant_scene = Scene()
    plant_scene.read(f, 'BGEOM')

    # multiply a plant with variations
    canopy, domain = create_heterogeneous_canopy(plant_scene, 
                                                 nplants=nplants, 
                                                 plant_density=plant_density, 
                                                 var_plant_position=var_plant_position)

    scenes.append(canopy)

# build the scene
geometry = {"scenes" : scenes }

lighting.build(geometry)

# compute lighting
lighting.run(energy=energy, hour=hour, day=day)

s = lighting.plantGL_light(printtriangles=True, printvoxels=False)

Viewer.display(s)