from src.LightVegeManager import *

if __name__ == "__main__":
    # les scènes en entrée
    s1 = pgl.Scene(r'scenes/scene_1ble.bgeom')
    s2 = pgl.Scene(r'scenes/scene_1luzerne_feuilles.geom')
    
    in_scenes = [s1, s2]
    # une liste [rescale, translation] par scène, "none" si pas de transformation
    in_trans = [[300, "none"], 
                ["none", Vector3(-30, -14, 0)]] 
    in_names = ["fspm-wheat", "l-egume"]
    
    # paramètres de RATP
    dx, dy, dz = 4, 4, 2 # cm
    latitude, longitude, timezone = 48.6, 7.76, 0.
    rs=[0.075, 0.2]
    tesselate_level = 3
    ratp_parameters = [dx, dy, dz, latitude, longitude, timezone, rs, tesselate_level]
    
    # Fichiers météo et données physiques des plantes
    meteopath = "inputs/meteo_exemple.mto"
    vegepath = "inputs/myveg_exemple.veg" 
    elevationpath = "inputs/distrib_leafinclination_exemples.dat"

    # Objet calcul de la lumière
    lght = LightVegeManager(in_scenes, in_trans, in_names, "ratp", ratp_parameters, tesselation=False)
    print(lght)
    # imprime la géométrie des scènes couplés (triangulation et voxels)
    lght.VTKinit("outputs/static_simple/lvm_")
    
    # lance le bilan radiatif
    lght.run(meteopath, vegepath, compute_classelevations=True)

    # imprime les résultats sur la triangulation
    lght.VTKout("outputs/static_simple/lvm_")
    
    # outils calculs des distributions sur la scène
    #lght.s5()
    #
    #lght.s2v()