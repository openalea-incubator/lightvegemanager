'''
valeurs par défaut de toutes les entrées des dico d'entrées
'''

def default_LightVegeManager_inputs() :
    default_environnement = {}
    # INRAE Lusignan
    default_environnement["coordinates"] = [46.4, 0., 1.]
    default_environnement["sky"] = "turtle46"
    default_environnement["diffus"] = True
    default_environnement["direct"] = True
    default_environnement["reflected"] = False
    default_environnement["infinite"] = True

    # pas encore de donnée géométrique par défaut
    # default_geometry = {}
    
    default_ratp_parameters = {}
    default_ratp_parameters["voxel size"] = [0.1, 0.1, 0.1]
    default_ratp_parameters["soil reflectance"] = [0., 0.]
    default_ratp_parameters["mu"] = [1.]
    default_ratp_parameters["tesselation level"] = [0]
    default_ratp_parameters["angle distrib algo"] = "compute global"
    default_ratp_parameters["nb angle class"] = 9
    
    default_caribu_parameters = {}
    default_caribu_parameters["sun algo"] = "caribu"
    default_caribu_parameters["soil mesh"] = -1
    default_caribu_parameters["debug"] = False

    return default_environnement,     \
            default_ratp_parameters,  \
            default_caribu_parameters \

