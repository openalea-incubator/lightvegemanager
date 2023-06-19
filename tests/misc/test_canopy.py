"""
    Test on 3D objects of Alfafa and fetuque

    - one plant of each with caribu
    - random generated canopy from a fetuque model
    - mixed canopy with fetuque and alfafa with caribu and ratp
"""
import os
import time
from openalea.plantgl.all import Scene
from openalea.mtg.mtg import MTG

from lightvegemanager.tool import LightVegeManager
from lightvegemanager.trianglesmesh import create_heterogeneous_canopy


def one_plant(bgeom_file, outputs_dir, environment={}, caribu_parameters={}, energy=500, day=264, hour=15, write_vtk=False) : 
    start = time.time()
    
    # create a plantGL scene from a bgeom file
    plant_scene = Scene()
    plant_scene.read(bgeom_file, 'BGEOM')

    # Initializing the tool
    lighting = LightVegeManager(lightmodel="caribu", 
                                    environment=environment, 
                                    lightmodel_parameters=caribu_parameters)

    # build the scene
    geometry = {"scenes" : [plant_scene]} # geometry is a plantGL Scene
    lighting.build(geometry)

    # compute lighting
    lighting.run(energy=energy, hour=hour, day=day)

    # print results gathered by elements (Shapes in the plantGL Scene)
    print(lighting.elements_outputs)

    # write outputs in a VTK file
    if write_vtk :
        filepath = os.path.join(outputs_dir, bgeom_file.split('/')[-1].split('.')[0])
        lighting.VTK_light(filepath)
    
    print("running time : ", time.time() - start, " s")

    print("\n")

def canopy_from_oneplant(bgeom_files, outputs_dir, 
                        nplants=20, plant_density=100, var_plant_position=100, id_type=None,
                        light_model="caribu", environment={}, lighting_parameters={}, 
                        energy=500, day=264, hour=15, 
                        write_vtk=False) : 
    start = time.time()

    # Initializing the tool
    lighting = LightVegeManager(lightmodel=light_model, 
                                    environment=environment, 
                                    lightmodel_parameters=lighting_parameters)    

    # generate random canopy from plant examples
    if not isinstance(bgeom_files, list): bgeom_files = [bgeom_files]
    scenes = []
    for f in bgeom_files :
        plant_scene = Scene()
        plant_scene.read(f, 'BGEOM')

        # multiply a plant with variations
        canopy, domain = create_heterogeneous_canopy(plant_scene, nplants=nplants, plant_density=plant_density, var_plant_position=var_plant_position, id_type=id_type)

        scenes.append(canopy)
    
    # build the scene
    geometry = {"scenes" : scenes }
    
    lighting.build(geometry)

    # compute lighting
    lighting.run(energy=energy, hour=hour, day=day)

    # print results gathered by elements (Shapes in the plantGL Scene)
    print(lighting.elements_outputs)

    # write outputs in a VTK file
    if write_vtk :
        filename = [x.split('/')[-1].split('.')[0] for x in bgeom_files]
        filepath = os.path.join(outputs_dir, "_".join(filename) +"_"+light_model+"_canopy")
        lighting.VTK_light(filepath)
    
    print("Light model running time : ", lighting.modelruntime, " s")
    print("Total running time : ", time.time() - start, " s")
   
    print("\n")


if __name__ == '__main__':
    print("--- START --- \n")
    
    output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "outputs", "test_canopy")
    
    # inputs values for lighting
    energy=500
    day=264
    hour=15

    #### ONE FETUQUES 
    print("--- ONE FETUQUE")
    
    # file path
    fet_fgeom = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "Fet-LD-F2.bgeom")

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

    # one_plant(fet_fgeom, output_dir, environment, caribu_parameters, energy=energy, day=day, hour=hour, write_vtk=True)

    #### ONLY LUZERNE
    print("--- ONE LUZERNE")
    # file path
    luz_fgeom = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "LD-F1.bgeom")

    # We add sky modelling
    n_azimuts = 4
    n_zenits = 5
    sky_type = "soc"
    environment["sky"] = [n_azimuts, n_zenits, sky_type]
    environment["diffus"] = True

    # we change the sun computing with default caribu algorithm
    caribu_parameters["sun algo"] = "ratp"

    # one_plant(luz_fgeom, output_dir, environment, caribu_parameters, energy=energy, day=day, hour=hour, write_vtk=True)

    #### MANY FETUQUES 
    print("--- SET OF FETUQUES")
    
    # desactivate sky modelling
    environment["diffus"] = False
    # canopy_from_oneplant(fet_fgeom, output_dir, 
    #                     nplants=20, plant_density=100, var_plant_position=100, id_type="plant",
    #                     light_model="caribu", environment=environment, lighting_parameters=caribu_parameters, 
    #                     energy=energy, day=day, hour=hour, 
    #                     write_vtk=True)

    #### MIX FETUQUE + LUZERNE

    # creates grid of voxels
    dv = 20.
    dxyz = [dv, dv, 10.] # en m
    nxyz = [10, 10, 4] # nombre de capteurs dans chaque axe
    orig = [-200.,-200.,0.]

    # ## CARIBU
    print("--- MIX FETUQUE + LUZERNE --- CARIBU")
    
    environment["diffus"] = True
    caribu_parameters["sensors"] = ["grid", dxyz, nxyz, orig, output_dir, "vtk"]

    # canopy_from_oneplant([fet_fgeom, luz_fgeom], output_dir, 
    #                     nplants=25, plant_density=120, var_plant_position=130,
    #                     light_model="caribu", environment=environment, lighting_parameters=caribu_parameters, 
    #                     energy=energy, day=day, hour=hour, 
    #                     write_vtk=True)

    ## RATP
    print("--- MIX FETUQUE + LUZERNE --- RATP")
    
    environment["diffus"] = False

    # Paramètres RATP #
    ratp_parameters = {}
    ratp_parameters["voxel size"] = dxyz
    ratp_parameters["angle distrib algo"] = "compute global"
    ratp_parameters["nb angle classes"] = 30
    canopy_from_oneplant([fet_fgeom, luz_fgeom], output_dir, 
                        nplants=25, plant_density=130, var_plant_position=110,
                        light_model="ratp", environment=environment, lighting_parameters=ratp_parameters, 
                        energy=energy, day=day, hour=hour, 
                        write_vtk=True)

    print("--- END")