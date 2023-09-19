from lightvegemanager.tool import LightVegeManager
from lightvegemanager.trianglesmesh import random_triangle_generator
import numpy

# grid dimensions
dxyz = [1.] * 3
nxyz = [7, 7, 7]
orig = [-1., -1., 0.]

spheresize = (1., 0.3) # vertices of triangles are the sphere surface
worldsize = (0., 5.)

nb_triangles = 10
triangles1 = [random_triangle_generator(worldsize=worldsize, spheresize=spheresize) for i in range(nb_triangles)]

nb_triangles = 9
triangles2 = [random_triangle_generator(worldsize=worldsize, spheresize=spheresize) for i in range(nb_triangles)]

nb_triangles = 8
triangles3 = [random_triangle_generator(worldsize=worldsize, spheresize=spheresize) for i in range(nb_triangles)]

scene = {0: triangles1, 1: triangles2, 2: triangles3}


riri_parameters = {"voxel size": [1., 1., 1.]}
l_scene = {"LA": numpy.ones([2, 3, 4, 4]), "distrib": [[0.1, 0.2, 0.1, 0.1, 0.05, 0.05, 0.1, 0.2, 0.1], 
                                                        [0.1, 0.1, 0.1, 0.1, 0.05, 0.1, 0.2, 0.2, 0.05]]}
for i in range(4):
    for j in range(4):
        l_scene["LA"][0][0][i][j] = 0.
        l_scene["LA"][1][0][i][j] = 0.

lighting = LightVegeManager(lightmodel="riri5", lightmodel_parameters=riri_parameters)
lighting.build(geometry={"scenes" : [l_scene] })

energy = 500.
hour = 15
day = 264
lighting.run(energy=energy, hour=hour, day=day)

res_abs_i = lighting.riri5_intercepted_light
res_trans = lighting.riri5_transmitted_light

print(res_abs_i[0])
print(res_trans)

lighting = LightVegeManager(lightmodel="riri5", lightmodel_parameters=riri_parameters)
lighting.build(geometry={"scenes" : [scene] })
lighting.run(energy=energy, hour=hour, day=day)

res_abs_i = lighting.riri5_intercepted_light
res_trans = lighting.riri5_transmitted_light