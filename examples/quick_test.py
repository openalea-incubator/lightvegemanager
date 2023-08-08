from lightvegemanager.tool import LightVegeManager

print("--- START")

# one triangle as a geometric element
triangle_vertices = [(0.,0.,0.), (1.,0.,0.), (1.,1.,1.)]

# surfacic lighting with CARIBU
lighting = LightVegeManager(lightmodel="caribu")

# build the scene
lighting.build(geometry=triangle_vertices)

# compute lighting
energy = 500
hour = 15
day = 264 # 21st september
lighting.run(energy, hour, day)

# output
print(lighting.elements_outputs)
print("--- END")