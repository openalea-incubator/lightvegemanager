from LightVegeManager import *

# one triangle as a geometric element
# we write our triangle in a CaribuScene format
organ_id = 1
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