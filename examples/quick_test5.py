from lightvegemanager.tool import LightVegeManager
from lightvegemanager.trianglesmesh import random_triangle_generator
from openalea.plantgl.all import Viewer

# grid dimensions
dxyz = [1.] * 3
nxyz = [5, 5, 7]
orig = [0.] * 3

# random triangles
nb_triangles = 50
spheresize = (1., 0.3) # vertices of triangles are the sphere surface
triangles = []
for i in range(nb_triangles) :
    triangles.append(random_triangle_generator(worldsize=(0., 5.), spheresize=spheresize))

caribu_args = { "sensors" : ["grid", dxyz, nxyz, orig] }

lighting = LightVegeManager(lightmodel="caribu", lightmodel_parameters=caribu_args, environment={"infinite":True})
lighting.build(geometry={"scenes" : [triangles]})

energy = 500.
hour = 15
day = 264
lighting.run(energy=energy, hour=hour, day=day)

s = lighting.plantGL_sensors(light=True)
Viewer.display(s)