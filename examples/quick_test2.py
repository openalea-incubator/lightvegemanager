from lightvegemanager.LVM import LightVegeManager
import openalea.plantgl.all as pgl

print("--- START")

# one triangle as a geometric element
triangle_vertices_1 = [[(0.,0.,0.), (1.,0.,0.), (1.,1.,1.)]]
triangle_vertices_2 = [[(0.,2.,0.), (1.,2.,0.), (1.,3.,1.)]]
stems = [(0,1)]
geometry = {
    "scenes" : [triangle_vertices_1, triangle_vertices_2],
    "stems id" : stems
}

# surfacic lighting with CARIBU
lighting = LightVegeManager(lightmodel="caribu")

# build the scene
lighting.build(geometry=geometry)

# compute lighting
energy = 500
hour = 15
day = 264 # 21st september
lighting.run(energy, hour, day)

# output
print(lighting.elements_outputs)

# visualisation
pglscene = lighting.plantGL_light()
pgl.Viewer.display(pglscene)

print("--- END")