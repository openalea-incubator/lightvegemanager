"""
    Test sur une canopée dense

    CARIBU
    RATP
"""
import os
from openalea.plantgl.all import Scene

from lightvegemanager.tool import LightVegeManager

print("--- START --- \n")

output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "outputs", "test_canopy")

####### CREATION OF A CANOPY #######
#### ONLY FETUQUE
print("--- ONE FETUQUE")
fet_fgeom = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "Fet-LD-F2.bgeom")
s_fet = Scene()
s_fet.read(fet_fgeom, 'BGEOM')

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
caribu_parameters["sun algo"] = "ratp"
caribu_parameters["caribu opt"] = {} 
caribu_parameters["caribu opt"]["par"] = (0.10, 0.05)

# Initializing the tool
lighting = LightVegeManager(lightmodel="caribu", 
                                environment=environment, 
                                lightmodel_parameters=caribu_parameters)

# build the scene
geometry = {"scenes" : [s_fet]} # geometry is a plantGL Scene
lighting.build(geometry)

# compute lighting
energy = 500
hour = 15
day = 264 # 21st september
lighting.run(energy=energy, hour=hour, day=day)

# print results gathered by elements (Shapes in the plantGL Scene)
print(lighting.elements_outputs)

# write outputs in a VTK file
filepath = os.path.join(output_dir, "fet_alone.vtk")
lighting.VTK_light(filepath)
print("\n")

#### MANY FETUQUES 
print("--- SET OF FETUQUES")


#### ONLY LUZERNE
luz_fgeom = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "LD-F1.bgeom")
s_luz = Scene()
s_luz.read(luz_fgeom, 'BGEOM')

#### FETUQUE + LUZERNE