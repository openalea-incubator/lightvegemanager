# import the class tool
from lightvegemanager.tool import LightVegeManager
from lightvegemanager.trianglesmesh import random_triangle_generator
import openalea.plantgl.all as pgl
import random

spheres = []

for i in range(255):
    x, y, z = [random.choice([-1, 1]) * random.randrange(0, int(255/2)) for i in range(3)]
    m = pgl.Material(ambient=(x+int(255/2), y+int(255/2), z+int(255/2)), shininess=0.1, diffuse=1)
    spheres.append(pgl.Shape(pgl.Translated(x, y, z, pgl.Sphere(random.randrange(1, 20), 20, 20)), m))

nb_triangles = 500
spheresize = (10., 2.)
triangles = []
for i in range(nb_triangles):
    triangles.append(random_triangle_generator(spheresize=spheresize))

# initialize the instance
lighting = LightVegeManager(lightmodel="caribu")

# build the scene
lighting.build(geometry=triangles)

# compute the lighting
energy = 500.
hour = 15
day = 264
lighting.run(energy=energy, hour=hour, day=day)

pglscene = lighting.plantGL_light()
pgl.Viewer.display(pglscene)