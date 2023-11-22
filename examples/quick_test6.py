from lightvegemanager.LVM import LightVegeManager
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


ratp_parameters = { "voxel size" : dxyz,
                    "origin" : orig,
                    "number voxels" : nxyz,
                    "full grid" : True}

lighting = LightVegeManager(lightmodel="ratp", lightmodel_parameters=ratp_parameters)
lighting.build(geometry={"scenes" : [scene] })

energy = 500.
hour = 15
day = 264
lighting.run(energy=energy, hour=hour, day=day)

m_lais = numpy.zeros([1] + nxyz)
for row in lighting.voxels_outputs.itertuples():
    m_lais[int(row.VegetationType)-1][row.Nz-1][row.Nx-1][row.Ny-1] = row.Area
res_abs_i, res_trans = lighting.to_l_egume(m_lais=m_lais)

print(res_abs_i[0])
print(res_trans)