import openalea.plantgl.all as pgl

from PyRATP.pyratp import grid
from PyRATP.pyratp import RATP2VTK
from PyRATP.pyratp.vegetation import Vegetation
from PyRATP.pyratp.skyvault import Skyvault
from PyRATP.pyratp.micrometeo import MicroMeteo
from PyRATP.pyratp.runratp import runRATP

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky
from alinea.caribu.sky_tools import GetLight
from alinea.caribu.sky_tools import Gensun
from alinea.caribu.sky_tools import GetLightsSun
from alinea.caribu.sky_tools import spitters_horaire
from alinea.caribu.sky_tools import turtle
from alinea.caribu.sky_tools import Sun

import os 
import subprocess
import itertools
import numpy as np
import pandas
import time
import scipy

from src.Polygons import *
from src.MyTesseletor import *

# Pour utilisation de RiRi (à partir de l-egume)
try:
    import RIRI5 as riri
except:
    print("Error : l-egume not installed")
    



'''
Fonctions de gestion
'''


def PlantGL_translation(s, t):
    """Translation d'une scène plantgl
    Args :
        s(Scene plantgl) : scène en entrée/sortie
        t : Vector3
    return :
        scene PlantGL
    """
    shapes_res = []

    for shape in s:
        shapes_res.append(pgl.Shape(pgl.Translated(t[0], t[1], t[2], shape.geometry), shape.appearance))
    return pgl.Scene(shapes_res)


def pgl_to_triangles(pgl_object, tesselator=None):
    """Transforme une scene PlantGL en triangulation
    inspiré de pgl_to_triangles dans CARIBU
    https://github.com/openalea-incubator/caribu : src/PyRATP/caribu/plantgl_adaptor.py
    
    Args :
        pgl_object : une shape plantGL
    
    return :
        liste de Triangle3
    """
    triangles = []
    if tesselator is None:
        tesselator = pgl.Tesselator()
    pgl_object.apply(tesselator)
    mesh = tesselator.triangulation
    if mesh:
        indices = mesh.indexList
        pts = list(map(tuple,mesh.pointList))
        triangles = [Triangle3(*(Vector3(*pts[itri[0]]),Vector3(*pts[itri[1]]),Vector3(*pts[itri[2]]))) for itri in indices]
    return triangles


def pgl_to_caribu(pgl_object, tesselator=None):
    """Transforme une scene PlantGL en triangulation
    inspiré de pgl_to_triangles dans CARIBU
    https://github.com/openalea-incubator/caribu : src/PyRATP/caribu/plantgl_adaptor.py
    
    Args :
        pgl_object : une shape plantGL
    
    return :
        liste de Triangle3
    """
    triangles = []
    if tesselator is None:
        tesselator = pgl.Tesselator()
    pgl_object.apply(tesselator)
    mesh = tesselator.triangulation
    if mesh:
        indices = mesh.indexList
        pts = list(map(tuple,mesh.pointList))
        triangles = [(pts[itri[0]],pts[itri[1]],pts[itri[2]]) for itri in indices]
    return triangles


def whichstems_MTG(MTG, nent):
    stems=[]
    geom = MTG.property('geometry')
    for vid in geom.keys():
        if MTG.class_name(vid) == 'StemElement':
            stems.append((vid, nent))
    
    return stems


def VTKtriangles(triangles, var, varname, filename):
    """Ecriture d'un fichier VTK à partir d'un maillage triangulation
    possibilité d'associer des grandeurs aux triangles
    
    Args :
        triangles : liste de Triangle3
        var : liste de grandeurs associées aux triangles, pour n triangles
             [
                 [var1_1, ..., var1_n], ... , [varm_1, ..., varm_n]
                 ]
        varname : liste de string, liste des noms pour chaque grandeurs
        filename : string, chemin du fichier à écrire

    """

    nbtr=0
    for tr in triangles:
        nbtr +=1
    
    f=open(filename, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS '+str(nbtr * 3)+' float\n')

    for tr in triangles:
        for i in range(3):
            f.write(str(tr[i][0])+' '+str(tr[i][1])+' '+str(tr[i][2])+'\n')
    f.write('\n')

    f.write('CELLS '+str(nbtr)+' '+str(nbtr*4)+'\n')
    for i in range(nbtr):
        f.write('3 '+str(3*i)+' '+str(1+3*i)+' '+str(2+3*i)+'\n')
    f.write('\n')

    f.write('CELL_TYPES '+str(nbtr)+'\n')
    for i in range(nbtr):
        f.write('5\n')
    f.write('\n')

    f.write('CELL_DATA '+str(nbtr)+'\n')
    f.write('FIELD FieldData '+str(len(varname)+2)+'\n')
    for i, name in enumerate(varname):
        f.write(name+' 1 '+str(nbtr)+' float\n')
        for j in range(nbtr):
            f.write(str(var[i][j])+'\n')
        f.write('\n')

    f.write('\n')

    f.write("ID 1 "+str(nbtr)+' float\n')
    for i in range(nbtr):
        f.write(str(i)+'\n')
    f.write('\n')
    f.write("ShapeLVM 1 "+str(nbtr)+' float\n')
    for tr in triangles:
        f.write(str(tr.id)+'\n')
    f.write('\n')
    

    f.close()

    var=[]
    varname=[]


def VTKline(start, end, filename):
    """Ecriture d'une ligne VTK
    
    Args :
        start : Vector3 point de départ
        end : Vector3 point d'arrivée
        filename : string, chemin du fichier à écrire
    """

    f=open(filename, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk Polygon\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS 2 float\n')
    f.write(str(start[0])+" "+str(start[1])+" "+str(start[2])+"\n")
    f.write(str(end[0])+" "+str(end[1])+" "+str(end[2])+"\n")
    f.write("LINES 1 3\n")
    f.write("2 0 1")

    f.close()


class LightVegeManager:
    '''
    Classe gestion de la lumière et couplage de plusieurs modèles de plantes
    '''
    units = {'mm': 0.001, 'cm': 0.01, 'dm': 0.1, 'm': 1, 'dam': 10, 'hm': 100,'km': 1000}

    def __init__(self,
                    environment={},
                    lightmodel="ratp",
                    lightmodel_parameters={},
                    main_unit="m"):
        '''Constructeur
        
        Arg :
            environment : dict
                "names" : liste de str, un nom par élément de scenes (non obligatoire)
                "coordinates" : liste de float, [latitude, longitude, timezone]
                "sky" : "turtle46" ou ["file", "filepath"] ou [nb_azimut, nb_zenith, "soc" ou "uoc"]
                "diffus" : boolean, si on active le rayonnement diffus
                "direct" : boolean, si on active le rayonnement direct
                "reflected" : boolean, si on active le rayonnement réfléchi
                "reflectance coefficients" : liste de listes, pour chaque scène
                        si "ratp" : [lambertien PAR, lambertien NIR]
                        si "caribu" : [reflectance PAR, transmittance PAR]
                "infinite" : boolean, si on veut un couvert un infini
            
            lightmodel : string, "ratp" ou "caribu", nom du modèle à utiliser
            lightmodel_parameters : dict
                si lightmodel == "caribu" :
                    "sun algo" : "caribu" ou "ratp", indique quel modèle fournira la position du soleil
                
                si lightmodel == "ratp" :
                    "voxel size" : liste de float, dimension d'un voxel [dx, dy, dz]
                    "soil reflectance" : liste de float, [reflectance PAR, reflectance NIR], coefficients de réflexion du sol
                    "mu" : liste de float, pour chaque scene indique un coefficient de clumping 
                    "tesselation level" : int, niveau de découpage des triangles pour correspondre à la grille
                    "angle distrib algo" : "compute global" ou "compute voxel" ou "file", comment sera calculé la distribution d'angles des feuilles
                    "nb angle classes" : float, si "compute global" ou "compute voxel", indique le nombre de classes d'angles
                    "angle distrib file" : si "file", chemin du fichier à lire

            global_scene_tesselate_level : tesselation de la scène couplée globale, quelque soit le modèle de lumière
            main_unit : string, indique l'unité d'échelle du couplage parmi la liste self.units

        '''
        self.__in_environment = environment
        self.__lightmodel = lightmodel
        self.__in_lightmodel_parameters = lightmodel_parameters
        self.__main_unit = main_unit

        # paramètres par défaut
        if "coordinates" not in self.__in_environment : self.__in_environment["coordinates"] = [46.4, 0. , 1.] # INRAE Lusignan 
        if "sky" not in self.__in_environment : self.__in_environment["sky"] = "turtle46" # ciel turtle à 46 directions
        if "diffus" not in self.__in_environment : self.__in_environment["diffus"] = True
        if "direct" not in self.__in_environment : self.__in_environment["direct"] = True
        if "reflected" not in self.__in_environment : self.__in_environment["reflected"] = False
        if "infinite" not in self.__in_environment : self.__in_environment["infinite"] = False
        if "reflectance coefficients" not in self.__in_environment : self.__in_environment["reflectance coefficients"] = []

        # création d'une scène plantGL
        # pas encore en place (manque d'abstraction pour le rescale)
        
        # création de la scène du modèle de lumière
        # RATP
        if self.__lightmodel == "ratp":
            # paramètres par défaut :
            if "voxel size" not in self.__in_lightmodel_parameters : self.__in_lightmodel_parameters["voxel size"] = [0.1, 0.1, 0.1]
            if "soil reflectance" not in self.__in_lightmodel_parameters : self.__in_lightmodel_parameters["soil reflectance"] = [0., 0.]
            if "mu" not in self.__in_lightmodel_parameters : self.__in_lightmodel_parameters["mu"] = [1.]
            if "tesselation level" not in self.__in_lightmodel_parameters : self.__in_lightmodel_parameters["tesselation level"] = 0
            if "angle distrib algo" not in self.__in_lightmodel_parameters : self.__in_lightmodel_parameters["angle distrib algo"] = "compute global"
            if "nb angle classes" not in self.__in_lightmodel_parameters : self.__in_lightmodel_parameters["nb angle classes"] = 9

            # construction du ciel:
            # si on a une liste de paramètres
            if len(self.__in_environment["sky"]) > 1 and type(self.__in_environment["sky"]) == list :
                if self.__in_environment["sky"][0] == "file" : self.__sky = Skyvault.read(self.__in_environment["sky"][1])
                
                elif len(self.__in_environment["sky"]) == 3 :
                    
                    self.__in_environment["sky"][1]
                    ele, azi, omega, pc = self._discrete_sky(self.__in_environment["sky"][0], self.__in_environment["sky"][1], self.__in_environment["sky"][2])
                    self.__sky = Skyvault.initialise(ele, azi, omega, pc)

            # initialisation par défaut
            elif self.__in_environment["sky"] == "turtle46":
                # par défaut ciel à 46 directions ("tortue")
                self.__sky = Skyvault.initialise()                
            
    
        elif self.__lightmodel == "caribu":
            # paramètres par défaut
            if "sun algo" not in self.__in_lightmodel_parameters : self.__in_lightmodel_parameters["sun algo"] = "caribu"
            
            # création du ciel
            # par défaut turtle 46 directions, soc
            if len(self.__in_environment["sky"]) > 1 and type(self.__in_environment["sky"]) == list :
                if len(self.__in_environment["sky"]) == 3 :
                    sky_string = GetLight.GetLight(GenSky.GenSky()(1., self.__in_environment["sky"][2], 
                                                                        self.__in_environment["sky"][0], 
                                                                        self.__in_environment["sky"][1]))  #: (Energy, soc/uoc, azimuts, zenits)
    
                    # Convert string to list in order to be compatible with CaribuScene input format
                    self.__sky = []
                    for string in sky_string.split('\n'):
                        if len(string) != 0:
                            string_split = string.split(' ')
                            t = tuple((float(string_split[0]), tuple((float(string_split[1]), float(string_split[2]), float(string_split[3])))))
                            self.__sky.append(t)
            
                elif self.__in_environment["sky"][0] == "file" : 
                    listGene = []
                    f = open(self.__in_environment["sky"][1])
                    ndir=int(f.readline().strip().split('\t')[0])
                    hmoy=np.zeros(ndir)
                    azmoy=np.zeros(ndir)
                    omega=np.zeros(ndir)
                    pc=np.zeros(ndir)
                    for n in range(ndir):
                        listGene.append(f.readline().strip().split('\t'))
                    tabGene=np.array(listGene)
                    tabGene = np.cast['float64'](tabGene)
                    hmoy=np.transpose(tabGene)[0]*math.pi / 180
                    azmoy=np.transpose(tabGene)[1]*math.pi / 180
                    omega=np.transpose(tabGene)[2]
                    pc=np.transpose(tabGene)[3]
                    f.close()
                    
                    self.__sky = []
                    for i, p in enumerate(pc) :
                        dir=[0,0,0]
                        dir[0] = dir[1] = math.sin(hmoy[i])
                        dir[0] *= math.cos(azmoy[i])
                        dir[1] *= math.sin(azmoy[i])
                        dir[2] = -math.cos(hmoy[i])
                        t = tuple((float(p), tuple((dir[0], dir[1], dir[2]))))
                        self.__sky.append(t)


            elif self.__in_environment["sky"] == "turtle46":
                turtle_list = turtle.turtle()
                self.__sky=[]
                for i,e in enumerate(turtle_list[0]):
                    t = tuple((e, tuple((turtle_list[2][i][0], turtle_list[2][i][1], turtle_list[2][i][2]))))
                    self.__sky.append(t)   
            else:
                raise ValueError("Unknown sky parameters : can be either 'turtle46' or ['file', 'filepath'] or [nb_azimut, nb_zenith, 'soc' or 'uoc'] ")

        else:
            raise ValueError("Unknown lightmodel : can be either 'ratp' or 'caribu' ")
            
    def init_scenes(self, geometry = {}, global_scene_tesselate_level=0):
        '''Initialisation des scènes géométriques
        
        Arg :
            geometry : dict
                "scenes" : liste de scènes PlantGL
                "transformation" : dict, indique les transformations à effectuer sur les scènes (non obligatoire)
                        "rescale" : liste de float, un facteur d'aggrandissement par scène
                        "translate" : liste de src.Polygons.vector3, un vecteur par scène
                        "scenes unit" : liste de src, une unité de la liste self.units par scène
                        "xyz orientation" : liste de src, pour chaque scène, "x+ = N" ou "x+ = S" ou "x+ = E" ou "x+ = W"
                "stems id" : liste de tuple, pour chaque tige on a tuple(numero de l'élément, indice de la scène) (non obligatoire)
                "domain" : liste de float, si on a un couvert infini, indique le domaine à reproduire
        
        Notes :
                Convention de coordonnées (type RATP scene) : x+ = N, y+ = W 
                valable dans self.__myscene et self.__caribu_scene

        '''
        self.__in_geometry = geometry

        legume_grid = False
        self.__id_legume_scene = -9999
        
        # paramètre par défaut des tiges
        if "stems id" not in self.__in_geometry : self.__in_geometry["stems id"] = []

        # couplage des scènes dans un plantGL commun et création du tableau des ids
        # [plantGL] -(+)-> [Triangle3] -> transformation sur les Triangle3 -> [Triangle3]
        # une scène plantgl par espèce
        self.__matching_ids = {}
        self.__my_scene=[]
        count=0
        
        # def sceneplantgl_to_triangles()
        for i_esp, scene in enumerate(self.__in_geometry["scenes"]) :
            # scene planteGL
            if isinstance(scene, pgl.Scene):
                for id, pgl_objects in scene.todict().items():
                    lastid = len(self.__my_scene)
                    
                    tri_list = list(itertools.chain(*[pgl_to_triangles(pgl_object) for pgl_object in pgl_objects]))
                    # on set l'id des triangles de la shape
                    for tr in tri_list:
                        tr.set_id(count)
                    self.__my_scene.extend(tri_list)
                    # on set le tableau des indices
                    self.__matching_ids[count] = [id, i_esp, list(range(lastid,lastid+len(tri_list)))]
                    count += 1
            
            # fichier VGX
            elif isinstance(scene, str) :
                if scene.split(".")[-1] == "vgx":
                    lastid = len(self.__my_scene)
                    f = open(scene, 'r')
                    lines = f.readlines()
                    for l in lines[1:]:                    
                        l_list = l.split("\t")

                        # on considère comme feuille les éléments R (RGB) != 42
                        if l_list[10] != '42' :
                            tr = Triangle3(*(Vector3(float(l_list[13]), float(l_list[14]), float(l_list[15])),
                                                Vector3(float(l_list[16]), float(l_list[17]), float(l_list[18])),
                                                Vector3(float(l_list[19]), float(l_list[20]), float(l_list[21]))))
                            
                            if tr.area > 0:
                                tr.set_id(count)
                                self.__my_scene.append(tr)
                        
                    f.close()
                    self.__matching_ids[count] = [count, i_esp, list( range( lastid, len(self.__my_scene) - lastid ) ) ]
                    count += 1

            # Grille RiRi (de l-egume)
            else: 
                # on active la prise de la grille
                legume_grid = True 

                # on retient l'id de la scène l-egume
                self.__id_legume_scene = i_esp

        # def transform_triangulation()
        # applique les transformations sur les triangles
        if "transformations" in self.__in_geometry :
            for i_esp in range(len(self.__in_geometry["scenes"])):
                if isinstance(self.__in_geometry["scenes"][i_esp], pgl.Scene):
                    for tr in self.__my_scene:
                        if self.__matching_ids[tr.id][1] == i_esp:
                            if "rescale" in  self.__in_geometry["transformations"] : tr.rescale(self.__in_geometry["transformations"]["rescale"][i_esp])
                            
                            if "translate" in  self.__in_geometry["transformations"]  : tr.translate(self.__in_geometry["transformations"]["translate"][i_esp])
                            
                            # la convention du repère xyz par rapport aux points cardinaux est précisée
                            # on ramène la scène à la convention x+ = N
                            # si non précisé, on ne change pas l'orientation de la scène
                            if "xyz orientation" in self.__in_geometry["transformations"]:
                                if self.__in_geometry["transformations"]["xyz orientation"][i_esp] == "x+ = S":
                                    tr.zrotate(180)
                                elif self.__in_geometry["transformations"]["xyz orientation"][i_esp] == "x+ = W":
                                    tr.zrotate(90)
                                elif self.__in_geometry["transformations"]["xyz orientation"][i_esp] == "x+ = E":
                                    tr.zrotate(-90)
                                # cas particulier pour l-egume
                                elif self.__in_geometry["transformations"]["xyz orientation"][i_esp] == "y+ = y-":
                                    y_id = 1
                                    tr.transform_axis(y_id, h=-1, t=self.__in_lightmodel_parameters["xy max"])
                                    # tr.zrotate(-90)
                            
                            if "scenes unit" in  self.__in_geometry["transformations"] :
                                # recupère la variable pour la lisibilité
                                scene_unit = self.__in_geometry["transformations"]["scenes unit"][i_esp]
                                
                                if (scene_unit != self.__main_unit) and (scene_unit in self.units):
                                    tr.rescale(self.units[scene_unit]/self.units[self.__main_unit])
                            
        # enregistre l'aire du plus grand triangle (indicateur par rapport au besoin de tesselation)
        self.__maxtrarea = 0.
        for tr in self.__my_scene:
            if tr.area > self.__maxtrarea: self.__maxtrarea = tr.area
        
        for tr in self.__my_scene:
            if tr.area > self.__maxtrarea: self.__maxtrarea = tr.area
        
        # active tesselation sur la scène globale (pour avoir une triangulation plus fine)
        if global_scene_tesselate_level > 0:
            new_tr_scene=[]
            for tr in self.__my_scene:
                level = 0
                isworking = iterate_triangles(tr, level, global_scene_tesselate_level, new_tr_scene)
            # copie de la nouvelle triangulation
            self.__my_scene = new_tr_scene
        
        # min-max de la scène
        xmax, xmin, ymax, ymin, zmax, zmin = -999999,999999,-999999,999999,-999999,999999    
        self.__triangleLmax  = -9999999
        if self.__my_scene :
            for tr in self.__my_scene:
                trxmax, trxmin, trymax, trymin, trzmax, trzmin = -999999,999999,-999999,999999,-999999,999999 
                for i in range(3) :
                    p = tr[i]
                    # min-max global à la scene
                    if p[0] > xmax :
                        xmax = p[0]
                    if p[0] < xmin:
                        xmin = p[0]
                    if p[1] > ymax :
                        ymax = p[1]
                    if p[1] < ymin:
                        ymin = p[1]
                    if p[2] > zmax :
                        zmax = p[2]
                    if p[2] < zmin:
                        zmin = p[2]
                    
                    # min-max local au triangle
                    if p[0] > trxmax :
                        trxmax = p[0]
                    if p[0] < trxmin:
                        trxmin = p[0]
                    if p[1] > trymax :
                        trymax = p[1]
                    if p[1] < trymin:
                        trymin = p[1]
                    if p[2] > trzmax :
                        trzmax = p[2]
                    if p[2] < trzmin:
                        trzmin = p[2]
                
                # recherche de la différence de la longueur max d'un triangle max(xmax-xmin, ymax-ymin, zmax-zmin)
                m = max(trxmax-trxmin, trymax-trymin, trzmax-trzmin)
                if m > self.__triangleLmax : self.__triangleLmax = m
        
        self.__pmax = Vector3(xmax, ymax, zmax)
        self.__pmin = Vector3(xmin, ymin, zmin)
        
        # RATP
        if self.__lightmodel == "ratp":
            if self.__in_lightmodel_parameters["voxel size"] == "dynamic" :
                dv = 5 * self.__triangleLmax
                dx, dy, dz = dv, dv, dv

            else:
                dx = self.__in_lightmodel_parameters["voxel size"][0]
                dy = self.__in_lightmodel_parameters["voxel size"][1]
                dz = self.__in_lightmodel_parameters["voxel size"][2]

            # on ajuste au besoin les min-max si la scène est plane pour avoir un espace 3D
            if xmin == xmax:
                xmax += dx
                xmin -= dx
            if ymin == ymax:
                ymax += dy
                ymin -= dy
            if zmin == zmax:
                zmax += dz
                zmin -= dz
            
            self.__pmax = Vector3(xmax, ymax, zmax)
            self.__pmin = Vector3(xmin, ymin, zmin)

            # définit une origine en Pmin
            if "origin" not in self.__in_lightmodel_parameters : xorig, yorig, zorig = self.__pmin[0], self.__pmin[1], -self.__pmin[2]
            else : 
                if len(self.__in_lightmodel_parameters["origin"]) == 2 :
                    xorig, yorig, zorig = self.__in_lightmodel_parameters["origin"][0], self.__in_lightmodel_parameters["origin"][1], -self.__pmin[2]
                elif len(self.__in_lightmodel_parameters["origin"]) == 3 :
                    xorig, yorig, zorig = self.__in_lightmodel_parameters["origin"][0], self.__in_lightmodel_parameters["origin"][1], self.__in_lightmodel_parameters["origin"][2]
            
            if not legume_grid :
                if self.__matching_ids:

                    # on sépare les tiges dans une nouvelle entité si il y a des shapes non tiges    
                    # def process_stems()
                    if "stems id" in self.__in_geometry and len(self.__in_geometry["stems id"]) < len(self.__matching_ids) :
                        mem_mu=[]
                        for stem in self.__in_geometry["stems id"]:
                            # on recherche la shape correspondante dans le dict des id
                            ids=0
                            search=True
                            while ids < len(self.__matching_ids) and search:
                                if (self.__matching_ids[ids][0], self.__matching_ids[ids][1]) == stem:
                                    search=False
                                else:
                                    ids += 1
                            # change le numéro d'entité dans le dict
                            tuple_temp = (self.__matching_ids[ids][0], 
                                            self.__matching_ids[ids][1] + len(self.__in_geometry["scenes"]),
                                            self.__matching_ids[ids][2])
                            self.__matching_ids[ids] = tuple_temp

                            # ajoute des propriétés optiques et mu sur la nouvelle entité tige
                            if stem[1] not in mem_mu:
                                self.__in_lightmodel_parameters["mu"].append(self.__in_lightmodel_parameters["mu"][stem[1]])
                                self.__in_environment["reflectance coefficients"].append(self.__in_environment["reflectance coefficients"][stem[1]])
                                mem_mu.append(stem[1])

                    ## distribution d'angles pour chaque entité
                    distrib = {}
                    distrib_glob=[]

                    # recherche du nb d'entité
                    nent=0
                    for key, val in self.__matching_ids.items():
                        if val[1]+1 > nent:
                            nent = val[1]+1

                    # calcul en dynamique
                    # ele_option = nombre de classes
                    # on fait le calcul de la distribution global avant la tesselation des triangles
                    # pour optimiser les calculs
                    if self.__in_lightmodel_parameters["angle distrib algo"] == "compute global" or \
                        self.__in_lightmodel_parameters["angle distrib algo"] == "compute voxel" :
                        
                        # compte le nombre de triangles par entité
                        t_nent_area=[]
                        for k in range(nent):
                            totA=0
                            for t in self.__my_scene:
                                if self.__matching_ids[t.id][1] == k:
                                    totA+=t.area
                            t_nent_area.append(totA)
                        
                        # on va jusqu'à 91° pour prendre en compte les plans
                        angles = list(np.linspace(90/self.__in_lightmodel_parameters["nb angle classes"], 91, self.__in_lightmodel_parameters["nb angle classes"]))

                        # pour chaque entité
                        for k in range(nent):
                            classes = [0] * self.__in_lightmodel_parameters["nb angle classes"]
                            # parcourt les triangles
                            for t in self.__my_scene:
                                # pour chaque triangle de l'entité
                                if self.__matching_ids[t.id][1] == k:
                                    # recherche de la classe
                                    i=0
                                    while i<self.__in_lightmodel_parameters["nb angle classes"]:
                                        if t.elevation < angles[i]:
                                            classes[i] += t.area
                                            # pour sortir de la boucle
                                            i=self.__in_lightmodel_parameters["nb angle classes"]+10
                                        i+=1

                            distrib_glob.append(classes)

                        # convertit en pourcentage
                        for n in range(nent):
                            for i in range(len(distrib_glob[n])):
                                distrib_glob[n][i] *= 1/t_nent_area[n]
                    
                    # lecture du fichier
                    # ele_option = chemin du fichier
                    elif self.__in_lightmodel_parameters["angle distrib algo"] == "file" : 
                        f_angle = open(self.__in_lightmodel_parameters["angle distrib file"], 'r')
                        for i in range(len(self.__in_geometry["scenes"])):
                            line = f_angle.readline()
                            distrib_glob.append([float(x) for x in line.split(',')[1:]])

                    distrib["global"] = distrib_glob

                    # nombre de voxels
                    if "number voxels" in self.__in_lightmodel_parameters :
                        nx = self.__in_lightmodel_parameters["number voxels"][0]
                        ny = self.__in_lightmodel_parameters["number voxels"][1]
                        nz = self.__in_lightmodel_parameters["number voxels"][2]
                    else :
                        nx = int((self.__pmax[0] - self.__pmin[0]) // dx)
                        ny = int((self.__pmax[1] - self.__pmin[1]) // dy)
                        if "grid slicing" in self.__in_lightmodel_parameters :
                            if self.__in_lightmodel_parameters["grid slicing"] == "ground = 0." :
                                nz = int((self.__pmax[2] - 0.) // dz)
                        else :
                            nz = int((self.__pmax[2] - self.__pmin[2]) // dz)
                        if (self.__pmax[0] - self.__pmin[0]) % dx > 0 : nx += 1
                        if (self.__pmax[1] - self.__pmin[1]) % dy > 0 : ny += 1
                        if (self.__pmax[2] - self.__pmin[2]) % dz > 0 : nz += 1
                    
                    # création de la grille
                    # si pas de rayonnement réfléchi on annules la réflexion du sol
                    if self.__in_environment["reflected"] :
                        mygrid = grid.Grid.initialise(nx, ny, nz, 
                                                        dx, dy, dz, 
                                                        xorig, yorig, zorig, 
                                                        self.__in_environment["coordinates"][0], self.__in_environment["coordinates"][1], self.__in_environment["coordinates"][2], 
                                                        nent, 
                                                        self.__in_lightmodel_parameters["soil reflectance"], 
                                                        toric=self.__in_environment["infinite"])
                    else :
                        mygrid = grid.Grid.initialise(nx, ny, nz, 
                                                        dx, dy, dz, 
                                                        xorig, yorig, zorig, 
                                                        self.__in_environment["coordinates"][0], self.__in_environment["coordinates"][1], self.__in_environment["coordinates"][2], 
                                                        nent, 
                                                        [0., 0.], 
                                                        toric=self.__in_environment["infinite"])

                    # subdivision des triangles pour matcher la grille
                    if self.__in_lightmodel_parameters["tesselation level"]>0:
                        # traite les triangles de la shape
                        new_tr_scene=[]
                        start=time.time()
                        for tr in self.__my_scene:
                            level = 0
                            isworking = iterate_trianglesingrid(tr, mygrid, level, self.__in_lightmodel_parameters["tesselation level"], new_tr_scene)
                        # print("tesselation time : ",time.time()-start)
                        # copie de la nouvelle triangulation
                        self.__my_scene = new_tr_scene
                        self.__tess_time = time.time() - start
                    
                    # préparation du fill
                    # pour chaque triangle, indice entité, x, y, z, aire, nitro
                    entity, barx, bary, barz, a, n = [],[], [], [], [], []
                    for tr in self.__my_scene:
                        bar = tr.barycenter
                        barx.append(bar[0])
                        bary.append(bar[1])
                        barz.append(bar[2])
                        
                        # id : id de la shape, val : [id shape en input, id de l'entité]
                        # c'est une tige on divise par 2 le LAD
                        if (self.__matching_ids[tr.id][0], self.__matching_ids[tr.id][1] - len(self.__in_geometry["scenes"])) in self.__in_geometry["stems id"] or \
                            (self.__matching_ids[tr.id][0], self.__matching_ids[tr.id][1]) in self.__in_geometry["stems id"] :
                            a.append(tr.area)
                        else:
                            a.append(tr.area)
                        
                        n.append(0.)
                        entity.append(self.__matching_ids[tr.id][1])

                    mygrid, matching = grid.Grid.fill_1(entity, barx, bary, barz, a, n, mygrid)
                    self.__tr_vox = matching 

                    # si la distribution d'angle est par voxel, on est
                    # obligé de prendre en compte les triangles après
                    # la tesselation
                    if self.__in_lightmodel_parameters["angle distrib algo"] == "compute voxel":
                        distrib_vox=[]
                        angles = list(np.linspace(90/self.__in_lightmodel_parameters["nb angle classes"], 91, self.__in_lightmodel_parameters["nb angle classes"]))
                        t_area=[]
                        # pour chaque voxel
                        for k in range(mygrid.nveg):
                            t_nent_area=[]
                            distrib_nent=[]
                            # pour chaque entité
                            for n in range(nent):
                                areatot=0
                                classes = [0] * self.__in_lightmodel_parameters["nb angle classes"]
                                # on parcourt le dico des correspondances triangle->voxel
                                for idt, idv in self.__tr_vox.items():
                                    # si on est dans le bon voxel et le triangle appartient a la bonne entité
                                    if idv == k and self.__matching_ids[self.__my_scene[int(idt)].id][1] == n:
                                        areatot += self.__my_scene[int(idt)].area
                                        # recherche de la classe
                                        i=0
                                        while i<self.__in_lightmodel_parameters["nb angle classes"]:
                                            if self.__my_scene[int(idt)].elevation < angles[i]:
                                                classes[i] += self.__my_scene[int(idt)].area
                                                # pour sortir de la boucle
                                                i=self.__in_lightmodel_parameters["nb angle classes"]+10
                                            i+=1
                                t_nent_area.append(areatot)
                                distrib_nent.append(classes)
                            t_area.append(t_nent_area)
                            distrib_vox.append(distrib_nent)

                        to_remove=[]
                        # convertit en pourcentage
                        for k in range(mygrid.nveg):
                            for n in range(nent):
                                # peut survenir si une entité n'est pas dans le voxel
                                if t_area[k][n] != 0 :
                                    for i in range(self.__in_lightmodel_parameters["nb angle classes"]):
                                        distrib_vox[k][n][i] *= 1/t_area[k][n]
                                
                                # enleve les entités en trop
                                else:
                                    to_remove.append((k,n))
                        for t in to_remove:
                            del distrib_vox[t[0]][t[1]]

                        distrib["voxel"] = distrib_vox
                    self.__ratp_distrib = distrib

                # si l'entrée est vide
                else :
                    # nombre de voxels
                    if "number voxels" in self.__in_lightmodel_parameters :
                        nx = self.__in_lightmodel_parameters["number voxels"][0]
                        ny = self.__in_lightmodel_parameters["number voxels"][1]
                        nz = self.__in_lightmodel_parameters["number voxels"][2]
                    else :
                        nx,ny,nz = 1,1,1
                    mygrid = grid.Grid.initialise(nx, ny, nz,
                                                    dx, dy, dz, 
                                                    xorig, yorig, zorig,
                                                    self.__in_environment["coordinates"][0], self.__in_environment["coordinates"][1], self.__in_environment["coordinates"][2], 
                                                    1, 
                                                    self.__in_lightmodel_parameters["soil reflectance"], 
                                                    toric=self.__in_environment["infinite"])
                    self.__ratp_distrib = {"global" : [[1.]]}
            
            # si une scene en entrée vient de l-egume
            # scene est un dict et contient :
            #   scene["LAD"] = m_lais / surf_refVOX
            #   scene["distrib"] = ls_dif
            else:
                # copie les paramètres d'entrée dans l'instance
                self.__ratp_distrib = {"global" : self.__in_geometry["scenes"][self.__id_legume_scene]["distrib"]}               

                # compte sur le nombre d'entité en plus de la grille l-egume
                nent_triangles = 0
                for key, val in self.__matching_ids.items():
                    if val[1]+1 > nent_triangles:
                        nent_triangles = val[1]+1
                    # corrige la numérotation des entités pour les triangles si on a une grille l-egume en entrée
                    val[1] += self.__in_geometry["scenes"][self.__id_legume_scene]["LA"].shape[0]
                
                # on ajoute au nombre d'entités dans l-egume
                nent = self.__in_geometry["scenes"][self.__id_legume_scene]["LA"].shape[0] + nent_triangles

                # dimensions de la grille l-egume
                nx = self.__in_geometry["scenes"][self.__id_legume_scene]["LA"].shape[3]
                ny = self.__in_geometry["scenes"][self.__id_legume_scene]["LA"].shape[2]
                nz = self.__in_geometry["scenes"][self.__id_legume_scene]["LA"].shape[1]

                # on regarde si la liste de triangles rentre dans la grille
                if self.__pmax[0] > nx*dx :
                    nx += int((self.__pmax[0] - nx*dx) // dx)+1
                if self.__pmax[1] > ny*dy :
                    ny += int((self.__pmax[1] - ny*dy) // dy)+1
                
                nb0 = 0  # nb de couches sans feuilles/LAI
                # la triangulation dépasse de la grille (normalement 126 * dz)
                if self.__pmax[2] > nz*dz :
                    nz += int((self.__pmax[2] - nz*dz) // dz)+1
                
                # sinon on enlève les couches vides
                else :
                    # combien/quelles lignes == 0 dans la grille l-egume
                    laicum = np.sum(self.__in_geometry["scenes"][self.__id_legume_scene]["LA"], axis=0)
                    laicumvert = np.sum(laicum, axis=(1, 2))
                    for i in range(len(laicumvert)):
                        if laicumvert[i] == 0.:
                            nb0 += 1
                        else:
                            break
                    nz = nz - nb0
                    
                    # on réajuste si la triangulation dépasse le nouveau nz
                    if self.__pmax[2] > nz*dz :
                        nz += int((self.__pmax[2] - nz*dz) // dz)+1          

                    self.__legume_nb0 = self.__in_geometry["scenes"][self.__id_legume_scene]["LA"].shape[1] - nz   

                # initialisation de la grille
                # si pas de rayonnement réfléchi on annules la réflexion du sol
                # on ajoute une couche au-dessus
                if self.__in_environment["reflected"] :
                    mygrid = grid.Grid.initialise(nx, ny, nz, 
                                                    dx, dy, dz,
                                                    xorig, yorig, zorig,
                                                    self.__in_environment["coordinates"][0], 
                                                    self.__in_environment["coordinates"][1], 
                                                    self.__in_environment["coordinates"][2], 
                                                    nent, 
                                                    self.__in_lightmodel_parameters["soil reflectance"], 
                                                    toric=self.__in_environment["infinite"])
                else :
                    mygrid = grid.Grid.initialise(nx, ny, nz, 
                                                    dx, dy, dz, 
                                                    xorig, yorig, zorig,
                                                    self.__in_environment["coordinates"][0], 
                                                    self.__in_environment["coordinates"][1], 
                                                    self.__in_environment["coordinates"][2], 
                                                    nent, 
                                                    [0., 0.], 
                                                    toric=self.__in_environment["infinite"])

                # remplissage des voxels avec la grille 
                # test : rotation de +90 du plan xy
                n_vox_per_nent = []
                for ne in range(nent-nent_triangles):
                    k=0
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                legume_iz = iz + self.__legume_nb0

                                # on force les voxels vides dans les couches non vides à être interprétés par RATP
                                S_voxel = max(1e-14, self.__in_geometry["scenes"][self.__id_legume_scene]["LA"][ne][legume_iz][iy][ix])

                                mygrid.kxyz[ny-(iy+1), ix, iz] = k + 1 #ajouter 1 pour utilisation f90
                                mygrid.numx[k] = (ny-(iy+1)) + 1 #ajouter 1 pour utilisation f90
                                mygrid.numy[k] = ix + 1 #ajouter 1 pour utilisation f90
                                mygrid.numz[k] = iz + 1 #ajouter 1 pour utilisation f90
                                mygrid.nume[ne,k] = ne + 1
                                mygrid.nje[k] = max(ne + 1, mygrid.nje[k])
                                mygrid.nemax = max(mygrid.nemax, mygrid.nje[k])

                                mygrid.leafareadensity[ne,k] += S_voxel / (dx * dy * dz)
                                mygrid.s_vt_vx[ne,k] += S_voxel
                                mygrid.s_vx[k] += S_voxel
                                mygrid.s_vt[ne] += S_voxel
                                mygrid.s_canopy += S_voxel
                                    
                                k=k+1
                    n_vox_per_nent.append(k)

                mygrid.nveg=max(n_vox_per_nent)
                mygrid.nsol=mygrid.njx*mygrid.njy   # Numbering soil surface areas
                for jx in range(mygrid.njx):
                    for jy in range(mygrid.njy):
                        mygrid.kxyz[jx,jy,mygrid.njz] = mygrid.njy * jx + jy + 1

                for k in range(mygrid.nveg):
                    for je in range(mygrid.nje[k]):
                        if je == 0:
                            # !!! on considère que l-egume définit une hauteur fixe de voxel
                            mygrid.volume_canopy[mygrid.nent] = mygrid.volume_canopy[mygrid.nent] + dx * dy * dz  # Incrementing total canopy volume
                        if  mygrid.s_vt_vx[je,k] > 0. :
                            mygrid.volume_canopy[mygrid.nume[je,k] - 1] = mygrid.volume_canopy[mygrid.nume[je,k] - 1] + dx * dy * dz
                            mygrid.voxel_canopy[mygrid.nume[je,k] - 1] = mygrid.voxel_canopy[mygrid.nume[je,k] - 1] + 1
                                
                if self.__my_scene :
                    # complète avec la liste de triangles
                    # pour chaque triangle, indice entité, x, y, z, aire, nitro
                    entity, barx, bary, barz, a, n = [],[], [], [], [], []
                    for tr in self.__my_scene:
                        bar = tr.barycenter
                        barx.append(bar[0])
                        bary.append(bar[1])
                        barz.append(bar[2])
                        
                        # id : id de la shape, val : [id shape en input, id de l'entité]
                        # c'est une tige on divise par 2 le LAD
                        if (self.__matching_ids[tr.id][0], self.__matching_ids[tr.id][1] - len(self.__in_geometry["scenes"])) in self.__in_geometry["stems id"] or \
                            (self.__matching_ids[tr.id][0], self.__matching_ids[tr.id][1]) in self.__in_geometry["stems id"] :
                            a.append(tr.area)
                        else:
                            a.append(tr.area)
                        
                        n.append(0.)
                        entity.append(self.__matching_ids[tr.id][1])

                    mygrid, matching = grid.Grid.fill_1(entity, barx, bary, barz, a, n, mygrid)
                    self.__tr_vox = matching 

            # on enregistre la grille dans les deux cas (plantGL ou l-egume)
            self.__ratp_scene = mygrid
                        
        # elif self.__lightmodel == "riri":
        #     self.__legume_triplets = riri.get_ls_triplets(self.__in_geometry["scenes"][self.__id_legume_scene]["LA"], opt=self.__in_geometry["scenes"][self.__id_legume_scene]["sky"][0])

        elif self.__lightmodel == "caribu":
            # copîe du minmax
            self.__pmax = Vector3(xmax, ymax, zmax)
            self.__pmin = Vector3(xmin, ymin, zmin)

            # création d'une triangulation caribu
            self.__caribu_scene = {}
            for id,val in self.__matching_ids.items():
                self.__caribu_scene[id] = []

            # convention de CARIBU avec x+ = N, pas de changement de base
            for tr in self.__my_scene:
                l_tr=[]
                for i in range(3):
                    l_tr.append((tr[i][0], tr[i][1], tr[i][2]))
                self.__caribu_scene[tr.id].append(l_tr)
        
        else:
            raise ValueError("Unknown lightmodel : can be either 'ratp' or 'caribu' ")

    def run(self, meteo_path="", energy=0, day=0, hour=0, parunit="micromol.m-2.s-1", truesolartime=False, printsun=False):
        '''calcul du bilan radiatif
        Args :
            meteo_path : string, chemin du fichier meteo
            energy : float, light energy en µmol.m-2.s-1 ou W.m-2
            day : float, day of the year 
            hour : float, hour TU
            parunit : "micromol.m-2.s-1" ou "W.m-2" ou  "RG"
            truesolartime : boolean, si hour est l'heure solaire réel ou l'heure locale
            printsun : option pour imprimer la position du soleil

        return :
            enregistre les sorties dans self.__outputs sous forme de dataframes pandas
            sorties en µmol.m-2.s-1
        '''
    
        # RATP
        if self.__lightmodel == "ratp":
            # création d'un dict entity
            entities_param = []
            if self.__in_lightmodel_parameters["angle distrib algo"] != "compute voxel":
                for id, mu_ent in enumerate(self.__in_lightmodel_parameters["mu"]):
                    if self.__in_environment["reflected"] :
                        entities_param.append({
                                                'mu' : mu_ent,
                                                'distinc' : self.__ratp_distrib["global"][id],
                                                'rf' : self.__in_environment["reflectance coefficients"][id]
                                                })
                    else :
                        entities_param.append({
                                                'mu' : mu_ent,
                                                'distinc' : self.__ratp_distrib["global"][id],
                                                'rf' : [0., 0.]
                                                })
                vegetation = Vegetation.initialise(entities_param)
            else :
                for id, mu_ent in enumerate(self.__in_lightmodel_parameters["mu"]):
                    if self.__in_environment["reflected"] :
                        entities_param.append({
                                                'mu' : mu_ent,
                                                'rf' : self.__in_environment["reflectance coefficients"][id]
                                                })
                    else :
                        entities_param.append({
                                            'mu' : mu_ent,
                                            'rf' : [0., 0.]
                                            })
                vegetation = Vegetation.initialise(entities_param, pervoxel=True, distribvox=self.__ratp_distrib["voxel"])

            # init météo, energy en entrée en W.m-2
            if meteo_path == "":
                if parunit == "micromol.m-2.s-1" :
                    #: Spitters's model estimating for the diffuse:direct ratio
                    # coefficient 2.02 : 4.6 (conversion en W.m-2) x 0.439 (PAR -> global)
                    RdRs = spitters_horaire.RdRsH(Rg=energy/2.02, DOY=day, heureTU=hour, latitude=self.__in_environment["coordinates"][0])
                    
                    # coeff 4.6 : https://www.researchgate.net/post/Can-I-convert-PAR-photo-active-radiation-value-of-micro-mole-M2-S-to-Solar-radiation-in-Watt-m2
                    energy = energy / 4.6  # W.m-2
                else:
                    RdRs = spitters_horaire.RdRsH(Rg=energy/0.439, DOY=day, heureTU=hour, latitude=self.__in_environment["coordinates"][0])
                
                # PAR et Dif en W.m^-2
                if self.__in_environment["diffus"] :
                    # Direct et diffus
                    if self.__in_environment["direct"]: 
                        met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=energy, Rdif=energy*RdRs, truesolartime=truesolartime)
                    # uniquement diffus
                    else:
                        met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=energy, Rdif=energy, truesolartime=truesolartime)
                # uniquement direct
                else :  met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=energy, Rdif=0, truesolartime=truesolartime)
            
            # lecture d'un fichier
            else:
                met = MicroMeteo.read(meteo_path, truesolartime)

            # impression de la position du soleil
            if printsun: 
                from PyRATP.pyratp import pyratp    
                self._print_sun(day, hour, pyratp, truesolartime)

            # grille RATP non vide
            if  self.__ratp_scene.nveg > 0 :
                start=time.time()
                # Calcul du bilan radiatif sur chaque pas de temps du fichier météo
                res = runRATP.DoIrradiation(self.__ratp_scene, vegetation, self.__sky, met)
                self.__time_runmodel = time.time() - start

                # Mise en forme des sorties
                # création de plusieurs tableaux intermédiaires qui serviront à trier les sorties
                entity = {}
                for id, match in self.__matching_ids.items():
                    entity[id] = match[1] + 1
                
                # si il y a une triangulation en entrée (défini à partir d'une scene plantGL)
                # if isinstance(self.__in_geometry["scenes"][self.__id_legume_scene], pgl.Scene):
                if self.__my_scene:
                    if self.__matching_ids:
                        index = range(len(self.__tr_vox))
                        vox_id = [self.__tr_vox[str(i)] + 1 for i in index]
                        # and one additional map that allows retrieving shape_id from python_x_index
                        sh_id=[]
                        for tr in self.__my_scene:
                            sh_id.append(tr.id)

                        s=[]
                        for tr in self.__my_scene:
                                s.append(tr.area)
                
                # récupère les sorties de RATP
                # np.array en une dimension, de taille nbvoxels x nbiteration
                VegetationType,Iteration,day,hour,VoxelId,ShadedPAR,SunlitPAR,ShadedArea,SunlitArea, xintav, Ptransmitted= res.T

                if parunit == "RG" :
                    # ('PAR' is expected in  Watt.m-2 in RATP input, whereas output is in micromol => convert back to W.m2 (cf shortwavebalance, line 306))
                    # on reste en micromol !
                    # On enregistre tout dans une dataframe pandas
                    para_list=[]
                    for i in range(len(ShadedPAR)):
                        if (ShadedArea[i] + SunlitArea[i]) > 0 :
                            para_list.append(((ShadedPAR[i]/4.6) * ShadedArea[i] + (SunlitPAR[i]/4.6) * SunlitArea[i]) / (ShadedArea[i] + SunlitArea[i]))
                        else:
                            para_list.append(0.)
                    
                    # vérifie qu'on a pas de erel négatif
                    erel_list=[]
                    for i in range(len(xintav)):
                        if xintav[i] >= 1e-6 :
                            erel_list.append(xintav[i])
                        else:
                            erel_list.append(0.)

                    # liste des indices
                    numx=[]
                    numy=[]
                    numz=[]
                    for v in VoxelId :
                        numx.append(self.__ratp_scene.numx[int(v)-1])
                        numy.append(self.__ratp_scene.numy[int(v)-1])
                        numz.append(self.__ratp_scene.numz[int(v)-1])

                    dfvox =  pandas.DataFrame({'VegetationType':VegetationType,
                                        'Iteration':Iteration,
                                        'day':day,
                                        'hour':hour,
                                        'VoxelId':VoxelId,
                                        'Nx':numx,
                                        'Ny':numy,
                                        'Nz':numz,
                                        'ShadedPAR':ShadedPAR/4.6,
                                        'SunlitPAR':SunlitPAR/4.6,
                                        'ShadedArea':ShadedArea,
                                        'SunlitArea': SunlitArea,
                                        'Area': ShadedArea + SunlitArea,
                                        'PARa': para_list,
                                        'intercepted': erel_list, 
                                        'transmitted': Ptransmitted
                                    })
                else:
                    para_list=[]
                    for i in range(len(ShadedPAR)):
                        if (ShadedArea[i] + SunlitArea[i]) > 0 :
                            para_list.append((ShadedPAR[i] * ShadedArea[i] + SunlitPAR[i] * SunlitArea[i]) / (ShadedArea[i] + SunlitArea[i]))
                        else:
                            para_list.append(0.)
                    
                    # vérifie qu'on a pas de erel négatif
                    erel_list=[]
                    for i in range(len(xintav)):
                        if xintav[i] >= 1e-6 :
                            erel_list.append(xintav[i])
                        else:
                            erel_list.append(0.)
                    dfvox =  pandas.DataFrame({'VegetationType':VegetationType,
                                        'Iteration':Iteration,
                                        'day':day,
                                        'hour':hour,
                                        'VoxelId':VoxelId,
                                        'Nx':self.__ratp_scene.numx[:self.__ratp_scene.nveg],
                                        'Ny':self.__ratp_scene.numy[:self.__ratp_scene.nveg],
                                        'Nz':self.__ratp_scene.numz[:self.__ratp_scene.nveg],
                                        'ShadedPAR':ShadedPAR,
                                        'SunlitPAR':SunlitPAR,
                                        'ShadedArea':ShadedArea,
                                        'SunlitArea': SunlitArea,
                                        'Area': ShadedArea + SunlitArea,
                                        'PARa': para_list,
                                        'xintav': erel_list, 
                                    })
                            
                # ne prend pas le sol
                dfvox = dfvox[dfvox['VegetationType'] > 0]
            else :
                dfvox =  pandas.DataFrame({'VegetationType':[0],
                                    'Iteration':[1],
                                    'day':day,
                                    'hour':hour,
                                    'VoxelId':[0],
                                    'Nx':[0],
                                    'Ny':[0],
                                    'Nz':[0],
                                    'ShadedPAR':[0], # /4.6,
                                    'SunlitPAR':[0], # /4.6,
                                    'ShadedArea':[0],
                                    'SunlitArea': [0],
                                    'Area': [0],
                                    'PARa': [0],
                                    'xintav': [0]
                                })
                
            # enregistre la dataframe des voxels
            self.__voxels_outputs = dfvox

            # si il y a une triangulation en entrée (défini à partir d'une scene plantGL)
            # if isinstance(self.__in_geometry["scenes"][id_legume], pgl.Scene):
            if self.__my_scene:
                if self.__matching_ids:
                    # nouvelle data frame avec les triangles en index
                    dfmap = pandas.DataFrame({'primitive_index': index,'shape_id': sh_id, 'VoxelId':vox_id, 'VegetationType':[entity[id] for id in sh_id], 'primitive_area':s})

                    # supposé copie dans les index avec des colonnes en commun
                    # colonnes en commun : VegetationType, VoxelId
                    output = pandas.merge(dfmap, dfvox)        
                    # tri les lignes par ordre de triangles
                    output =  output.sort_values('primitive_index')

                    # enregistre la dataframe dans l'instance
                    self.__outputs = output

                    # enregistre les valeurs par shape et plantes
                    #nshapes = sum([len(l) for l in self.__in_scenes])
                    nshapes = len(self.__matching_ids)
                    s_shapes = []
                    s_area=[]
                    s_para=[]
                    s_pari=[]
                    s_xintav=[]
                    s_ite=[]
                    s_day=[]
                    s_hour=[]
                    s_ent=[]
                    s_parsun=[]
                    s_parsha=[]
                    s_areasun=[]
                    s_areasha=[]
                    for id in range(nshapes):
                        # itérations commencent à 1
                        nent = self.__matching_ids[id][1]
                        for ite in range(int(max(output["Iteration"]))):
                            dffil = output[(output.Iteration == ite+1) & (output.shape_id == id)]
                            
                            s_hour.append(dffil["hour"].values[0])
                            s_day.append(dffil["day"].values[0])
                            s_ite.append(ite+1)
                            s_area.append(sum(dffil["primitive_area"]))

                            # le rayonnemenbt réfléchi est calculé dans RATP
                            if self.__in_environment["reflected"] :
                                s_para.append(sum(dffil['primitive_area']*dffil['PARa']) / s_area[-1])
                            
                            # sinon on le rajoute manuellement
                            else:
                                # PAR incident
                                s_pari.append(sum(dffil['primitive_area']*dffil['PARa']) / s_area[-1])
                                s_para.append(s_pari[-1] - (s_pari[-1] * self.__in_environment["reflectance coefficients"][nent][0]))
                                
                                # # applique une réflectance et une transmittance sur le PAR incident
                                # if (self.__matching_ids[id][0], nent - len(self.__in_geometry["scenes"])) in self.__in_geometry["stems id"] or \
                                #     (self.__matching_ids[id][0], nent) in self.__in_geometry["stems id"] :
                                #     
                                # else:
                                #     s_para.append(s_pari[-1] - (s_pari[-1] * self.__in_environment["reflectance coefficients"][nent][0] + \
                                #                                 s_pari[-1] * self.__in_environment["reflectance coefficients"][nent][1]))

                            s_parsun.append(sum(dffil['primitive_area']*dffil['SunlitPAR']) / s_area[-1])
                            s_parsha.append(sum(dffil['primitive_area']*dffil['ShadedPAR']) / s_area[-1])
                            s_areasun.append(sum(dffil['primitive_area']*dffil['SunlitArea']) / s_area[-1])
                            s_areasha.append(sum(dffil['primitive_area']*dffil['ShadedArea']) / s_area[-1])
                            s_xintav.append(sum(dffil['primitive_area']*dffil['xintav']) / s_area[-1])
                            s_ent.append(dffil["VegetationType"].values[0])
                            s_shapes.append(self.__matching_ids[id][0])
                    self.__shape_outputs = pandas.DataFrame({
                        "Iteration" : s_ite,
                        "Day" : s_day,
                        "Hour" : s_hour,
                        "ShapeId" : s_shapes,
                        "VegetationType" : s_ent,
                        "Area" : s_area,
                        "PARi" : s_pari,
                        "PARa" : s_para,
                        "xintav" : s_xintav,
                        "SunlitPAR" : s_parsun,
                        "SunlitArea" : s_areasun,
                        "ShadedPAR" : s_parsha,
                        "ShadedArea" : s_areasha
                    })

        elif self.__lightmodel == "riri":
            dS = self.__in_lightmodel_parameters["voxel size"][0] * self.__in_lightmodel_parameters["voxel size"][1]
            tag_light_inputs = [self.__in_geometry["scenes"][self.__id_legume_scene]["LA"] / dS, 
                                self.__in_geometry["scenes"][self.__id_legume_scene]["triplets"], 
                                self.__in_geometry["scenes"][self.__id_legume_scene]["distrib"], 
                                energy * dS]  # input tag

            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            self.__legume_transmitted_light, \
                self.__legume_intercepted_light = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, 
                                                                            optsky=self.__in_environment["sky"][0], 
                                                                            opt=self.__in_environment["sky"][1])

        # CARIBU
        elif self.__lightmodel == "caribu":
            if self.__matching_ids:
                ## Création du ciel et du soleil
                #: Direct light sources (sun positions)
                if self.__in_lightmodel_parameters["sun algo"]=="ratp":
                    from PyRATP.pyratp import pyratp
                    
                    az, ele = 5.,9. # variables fantômes (non récupérées)
                    pyratp.shortwave_balance.sundirection(ele, az, 
                                                                self.__in_environment["coordinates"][0], 
                                                                self.__in_environment["coordinates"][1], 
                                                                self.__in_environment["coordinates"][2], 
                                                                day, hour, truesolartime)            
                    degtorad = math.pi/180

                    azrad = pyratp.shortwave_balance.azdeg
                    # peut avoir un nan à 12h 
                    if math.isnan(azrad) and hour==12. : 
                        # critère empirique basé sur l'algo de CARIBU (GenSun)
                        if self.__in_environment["coordinates"][0] >= 21.11:
                            azrad = 0.
                        else:
                            azrad = 180.
                    # passage de l'azimuth en South=0° clockwise (RATP) à North=0° clockwise (CARIBU.Gensun)
                    azrad = -azrad * degtorad

                    zenirad = pyratp.shortwave_balance.hdeg * degtorad

                    # le vecteur pointe du ciel vers le sol
                    sunx = suny = math.cos(zenirad)
                    sunx *= math.cos(azrad)
                    suny *= math.sin(azrad)
                    sunz = -math.sin(zenirad)
                    
                    # list soleil
                    self.__sun = [tuple((1., tuple((sunx, suny, sunz))))]
                
                # algo de CARIBU
                else:
                    # si l'heure est locale, calcul de l'heure solaire (algorithme pris dans RATP, sundirection: shortwave_balance.f90)
                    if not truesolartime:
                        om =0.017202*(day-3.244)
                        teta=om+0.03344*math.sin(om)*(1+0.021*math.cos(om))-1.3526
                        tphi=0.91747*math.sin(teta)/math.cos(teta)
                        dphi=math.atan(tphi)-om+1.3526
                        if dphi+1. <= 0: dphi=(dphi+1+1000.*math.pi % math.pi)-1.
                        eqntime=dphi*229.2
                        hour =hour+self.__in_environment["coordinates"][2]+self.__in_environment["coordinates"][1]/15.-eqntime/60. % 24.

                    # Conversion de la latitude en radian!
                    sun = Gensun.Gensun()(1., day, hour, self.__in_environment["coordinates"][0]*math.pi/180)
                    sun = GetLightsSun.GetLightsSun(sun)
                    sun_str_split = sun.split(' ')
                    
                    # format d'entrée de CaribuScene
                    self.__sun = [tuple((float(sun_str_split[0]), tuple((float(sun_str_split[1]), float(sun_str_split[2]), float(sun_str_split[3])))))]

                # affichage des coordonnées du soleil
                if printsun: 
                    from PyRATP.pyratp import pyratp    
                    self._print_sun(day, hour, pyratp, truesolartime)
                    
                # active le calcul si le soleil est levé
                # critère élévation min fixé à 2° comme dans RATP (shortwave_balance.F90 -> DirectBeam_Interception, l.68 ou l.148)
                sun_up = math.asin(-self.__sun[0][1][2])*(180/math.pi) > 2.
                if sun_up:
                    # conversion en µmol/m-2.s-1
                    if parunit == "W.m-2": PAR=PAR*4.6

                    # opt (dict): a {band_name: {primitive_id: material}} dict of dict
                    #             or a {band_name: material} dict of tuples.
                    #             In the second form the material is used for all primitives.
                    #             A material is a 1-, 2- or 4-tuple depending on its optical behavior.
                    #             A 1-tuple encode an opaque material characterised by its reflectance
                    #             A 2-tuple encode a symmetric translucent material defined
                    #             by a reflectance and a transmittance
                    #             A 4-tuple encode an asymmetric translucent material defined
                    #             the reflectance and transmittance
                    #             of the upper and lower side respectively
                    #             Alternatively, a list of band_name.opt files (require scene
                    #             to be given as a *.can file)
                    #             If None (default), all primitive are associated to the
                    #             default material of the class.

                    #: Optical properties
                    # opt = {'par': {}}
                    # for id, val in self.__matching_ids.items():
                    #     # id : id de la shape, val : [id shape en input, id de l'entité]
                    #     # c'est une tige il n'y a pas de transmission
                    #     if (val[0], val[1]) in self.__in_geometry["stems id"]:
                    #         opt['par'][id] = (self.__in_environment["reflectance coefficients"][val[1]][0],)  #: (reflectance,) of the stems ou organe opaque
                    #     else:
                    #         opt['par'][id] = (self.__in_environment["reflectance coefficients"][val[1]][0], self.__in_environment["reflectance coefficients"][val[1]][1]) #: (reflectance, transmittance) of the adaxial side of the leaves, élément translucide symétrique

                    opt = {}
                    for band, coef in self.__in_environment["caribu opt"].items() :
                        opt[band] = {}
                    for id, val in self.__matching_ids.items():
                        for band, coef in self.__in_environment["caribu opt"].items() :
                            # id : id de la shape, val : [id shape en input, id de l'entité]
                            # c'est une tige il n'y a pas de transmission
                            if (val[0], val[1]) in self.__in_geometry["stems id"]:
                                opt[band][id] = (coef[0],)  #: (reflectance,) of the stems ou organe opaque
                            else:
                                # A 2-tuple encode a symmetric translucent material defined
                                # by a reflectance and a transmittance
                                # A 4-tuple encode an asymmetric translucent material defined
                                # the reflectance and transmittance
                                # of the upper and lower side respectively 
                                #: (reflectance, transmittance) of the adaxial side of the leaves, élément translucide symétrique
                                opt[band][id] = coef 

                    debug = False
                    if "debug" in self.__in_lightmodel_parameters and self.__in_lightmodel_parameters["debug"] : debug = True
                    # si on souhaite construire une grille de capteurs
                    if "sensors" in self.__in_lightmodel_parameters and self.__in_lightmodel_parameters["sensors"][0] == "grid" :
                        sensors_caribu, sensors_plantgl, Pmax_capt = self._create_caribu_legume_sensors()

                    # construction d'une scène ciel et soleil
                    if self.__in_environment["infinite"] : # on ajoute un domaine pour la création du pattern
                        c_scene_sky = CaribuScene(scene=self.__caribu_scene, light=self.__sky, opt=opt, scene_unit=self.__main_unit, pattern=self.__in_geometry["domain"], debug = debug)
                        c_scene_sun = CaribuScene(scene=self.__caribu_scene, light=self.__sun, opt=opt, scene_unit=self.__main_unit, pattern=self.__in_geometry["domain"], debug = debug)
                    else:
                        c_scene_sky = CaribuScene(scene=self.__caribu_scene, light=self.__sky, opt=opt, scene_unit=self.__main_unit, debug = debug)
                        c_scene_sun = CaribuScene(scene=self.__caribu_scene, light=self.__sun, opt=opt, scene_unit=self.__main_unit, debug = debug)

                    if self.__in_environment["diffus"] :
                        # Direct et diffus
                        if self.__in_environment["direct"] :
                            sun_sky_option = "mix"
                        else:
                            sun_sky_option = "sky"
                    else:
                        sun_sky_option = "sun"

                    # rayonnement réfléchis
                    if self.__in_environment["reflected"] : direct_active = False
                    else : direct_active = True
                    
                    # choix du type de rayonnement
                    # direct=False : active la rediffusion
                    # infinite=False : désactive la répétition infinie de la scène (pas de pattern défini ici)
                    if sun_sky_option == "mix":
                        if "sensors" in self.__in_lightmodel_parameters :
                            start=time.time()
                            raw_sun, aggregated_sun = c_scene_sun.run(direct=direct_active, infinite=self.__in_environment["infinite"], sensors=sensors_caribu)
                            raw_sky, aggregated_sky = c_scene_sky.run(direct=direct_active, infinite=self.__in_environment["infinite"], sensors=sensors_caribu)
                            self.__time_runmodel = time.time() - start
                        else :
                            start=time.time()
                            raw_sun, aggregated_sun = c_scene_sun.run(direct=direct_active, infinite=self.__in_environment["infinite"])
                            raw_sky, aggregated_sky = c_scene_sky.run(direct=direct_active, infinite=self.__in_environment["infinite"])
                            self.__time_runmodel = time.time() - start
                                               
                        #: Spitters's model estimating for the diffuse:direct ratio
                        # % de sky dans la valeur d'énergie finale
                        Rg = energy / 2.02  #: Global Radiation (W.m-2)
                        RdRs = spitters_horaire.RdRsH(Rg=Rg, DOY=day, heureTU=hour, latitude=self.__in_environment["coordinates"][0])  #: Diffuse fraction of the global irradiance
                        
                        PARa_output_shape = {}
                        PARi_output_shape = {}
                        for band, coef in self.__in_environment["caribu opt"].items() : 
                            Erel = {}
                            Ei_output_shape = {}
                            for element_id, Erel_value in aggregated_sky[band]['Eabs'].items():
                                Erel[element_id] = RdRs * Erel_value + (1 - RdRs) * aggregated_sun[band]['Eabs'][element_id]
                                Ei_output_shape[element_id] = RdRs * aggregated_sky[band]['Ei'][element_id] + (1 - RdRs) * aggregated_sun[band]['Ei'][element_id]

                            PARa_output_shape[band] = {k: v * energy for k, v in Erel.items()}
                            PARi_output_shape[band] = {k: v * energy for k, v in Ei_output_shape.items()}
                        
                        # dataframe des triangles
                        # même numération que caribu
                        PARa_output_tr = {}
                        PARi_output_tr = {}
                        for band, coef in self.__in_environment["caribu opt"].items() : 
                            Erel_output_tr = {}
                            Ei_output_tr = {}
                            count=0
                            for key, val in raw_sky[band]['Eabs'].items():
                                for i,par in enumerate(val):
                                    Erel_output_tr[count] = RdRs*par+(1-RdRs)*raw_sun[band]['Eabs'][key][i]
                                    Ei_output_tr[count] = RdRs*raw_sky[band]['Ei'][key][i]+(1-RdRs)*raw_sun[band]['Ei'][key][i]
                                    count+=1
                            PARa_output_tr[band] = {k: v * energy for k, v in Erel_output_tr.items()}
                            PARi_output_tr[band] = {k: v * energy for k, v in Ei_output_tr.items()}
                    
                    elif sun_sky_option == "sun":
                        if "sensors" in self.__in_lightmodel_parameters :
                            start=time.time()
                            raw_sun, aggregated_sun = c_scene_sun.run(direct=direct_active, infinite=self.__in_environment["infinite"], sensors=sensors_caribu)
                            self.__time_runmodel = time.time() - start
                        else :
                            start=time.time()
                            raw_sun, aggregated_sun = c_scene_sun.run(direct=direct_active, infinite=self.__in_environment["infinite"])
                            self.__time_runmodel = time.time() - start
                        
                        PARa_output_shape = {}
                        PARi_output_shape = {}
                        for band, coef in self.__in_environment["caribu opt"].items() :                 
                            PARa_output_shape[band] = {k: v * energy for k, v in aggregated_sun[band]['Eabs'].items()}
                            PARi_output_shape[band] = {k: v * energy for k, v in aggregated_sun[band]['Ei'].items()}
                        
                        # dataframe des triangles
                        PARa_output_tr = {}
                        PARi_output_tr = {}
                        for band, coef in self.__in_environment["caribu opt"].items() :  
                            Erel_output_tr = {}
                            Ei_output_tr = {}
                            count=0
                            for key, val in raw_sun[band]['Eabs'].items():
                                for i,par in enumerate(val):
                                    Erel_output_tr[count] = raw_sun[band]['Eabs'][key][i]
                                    Ei_output_tr[count] = raw_sun[band]['Ei'][key][i]
                                    count+=1
                            PARa_output_tr[band] = {k: v * energy for k, v in Erel_output_tr.items()}
                            PARi_output_tr[band] = {k: v * energy for k, v in Ei_output_tr.items()}

                        # retient les résultats des capteurs
                        if "sensors" in self.__in_lightmodel_parameters :
                            self.__sensors_outputs = {}
                            for band, coef in self.__in_environment["caribu opt"].items() :  
                                self.__sensors_outputs[band] = aggregated_sun[band]['sensors']['Ei']
                    
                    elif sun_sky_option == "sky":
                        if "sensors" in self.__in_lightmodel_parameters :
                            start=time.time()
                            raw_sky, aggregated_sky = c_scene_sky.run(direct=direct_active, infinite=self.__in_environment["infinite"], sensors=sensors_caribu)
                            self.__time_runmodel = time.time() - start
                        else :
                            start=time.time()
                            raw_sky, aggregated_sky = c_scene_sky.run(direct=direct_active, infinite=self.__in_environment["infinite"])
                            self.__time_runmodel = time.time() - start
                        
                        PARa_output_shape = {}
                        PARi_output_shape = {}
                        for band, coef in self.__in_environment["caribu opt"].items() :                 
                            PARa_output_shape[band] = {k: v * energy for k, v in aggregated_sky[band]['Eabs'].items()}
                            PARi_output_shape[band] = {k: v * energy for k, v in aggregated_sky[band]['Ei'].items()}
                        
                        # dataframe des triangles
                        PARa_output_tr = {}
                        PARi_output_tr = {}
                        for band, coef in self.__in_environment["caribu opt"].items() :  
                            Erel_output_tr = {}
                            Ei_output_tr = {}
                            count=0
                            for key, val in raw_sky[band]['Eabs'].items():
                                for i,par in enumerate(val):
                                    Erel_output_tr[count] = raw_sky[band]['Eabs'][key][i]
                                    Ei_output_tr[count] = raw_sky[band]['Ei'][key][i]
                                    count+=1
                            PARa_output_tr[band] = {k: v * energy for k, v in Erel_output_tr.items()}
                            PARi_output_tr[band] = {k: v * energy for k, v in Ei_output_tr.items()}

                        # retient les résultats des capteurs
                        if "sensors" in self.__in_lightmodel_parameters :
                            self.__sensors_outputs = {}
                            for band, coef in self.__in_environment["caribu opt"].items() :  
                                self.__sensors_outputs[band] = aggregated_sky[band]['sensors']['Ei']
                        
                    else:
                        raise ValueError("Unknown sun_sky_option : can be either 'mix', 'sun' or 'sky'.")

                    # affichage VTK des capteurs
                    if "sensors" in self.__in_lightmodel_parameters and self.__in_lightmodel_parameters["sensors"][-1] == "vtk":
                        triangles_sensors = []
                        for id, s in sensors_plantgl.todict().items() :
                            tri_list = list(itertools.chain(*[pgl_to_triangles(pgl_object) for pgl_object in s]))
                            for tr in tri_list:
                                tr.set_id(id)
                            triangles_sensors.extend(tri_list)
                        
                        var=[]
                        for t in triangles_sensors:
                            var.append(aggregated_sky["par"]['sensors']['Ei'][t.id])
                        VTKtriangles(triangles_sensors, [var], ["par_t"], self.__in_lightmodel_parameters["sensors"][-2] + "sensors.vtk")

                # enregistre les valeurs par shape et plantes
                s_shapes = [0]*len(self.__matching_ids)
                s_area = [0]*len(self.__matching_ids)
                s_Eabs = {}
                s_Ei = {}
                s_day = [0]*len(self.__matching_ids)
                s_hour = [0]*len(self.__matching_ids)
                s_ent = [0]*len(self.__matching_ids)

                for key,val in self.__matching_ids.items():
                    s_shapes[key] = val[0]
                    if sun_up: 
                        tr_shape = [self.__my_scene[t] for t in val[2]]
                        s_area[key] = sum([tr.area for tr in tr_shape])
                    s_day[key] = day
                    s_hour[key] = hour
                    s_ent[key] = self.__matching_ids[key][1]
                
                for band, coef in self.__in_environment["caribu opt"].items() : 
                    s_Eabs[band] = [0]*len(self.__matching_ids)
                    s_Ei[band] = [0]*len(self.__matching_ids)
                    for key,val in self.__matching_ids.items():
                        if sun_up: s_Eabs[band][key] = PARa_output_shape[band][key]
                        if sun_up: s_Ei[band][key] = PARi_output_shape[band][key]
                dico_shape = {
                    "Day" : s_day,
                    "Hour" : s_hour,
                    "ShapeId" : s_shapes,
                    "VegetationType" : s_ent,
                    "Area" : s_area
                }
                for band, coef in self.__in_environment["caribu opt"].items() : 
                    dico_shape[band+" Eabs"] = s_Eabs[band]
                    dico_shape[band+" Ei"] = s_Ei[band]
                self.__shape_outputs = pandas.DataFrame(dico_shape)

                # enregistre les valeurs par triangles
                s_shapes = [0]*len(self.__my_scene)
                s_tr = [0]*len(self.__my_scene)
                s_area=[0]*len(self.__my_scene)
                s_Eabs = {}
                s_Ei = {}
                s_day=[0]*len(self.__my_scene)
                s_hour=[0]*len(self.__my_scene)
                s_ent=[0]*len(self.__my_scene)

                for id,tr in enumerate(self.__my_scene):
                    s_shapes[id] = self.__matching_ids[tr.id][0]
                    s_tr[id] = id
                    s_day[id] = day
                    s_hour[id] = hour
                    s_ent[id] = self.__matching_ids[tr.id][1]
                    s_area[id] = tr.area    
                for band, coef in self.__in_environment["caribu opt"].items() : 
                    s_Eabs[band] = [0]*len(self.__my_scene)
                    s_Ei[band] = [0]*len(self.__my_scene)
                    for id,tr in enumerate(self.__my_scene):
                        if sun_up: s_Eabs[band][id] = PARa_output_tr[band][id]
                        if sun_up: s_Ei[band][id] = PARi_output_tr[band][id]            

                dico_tr = {
                    "Day" : s_day,
                    "Hour" : s_hour,
                    "primitive_index" : s_tr,
                    "ShapeId" : s_shapes,
                    "VegetationType" : s_ent,
                    "Area" : s_area
                }

                for band, coef in self.__in_environment["caribu opt"].items() : 
                    dico_tr[band+" Eabs"] = s_Eabs[band]
                    dico_tr[band+" Ei"] = s_Ei[band]
                self.__outputs = pandas.DataFrame(dico_tr)
            
            # si entrée géo vide
            else:
                dico_shape = {
                    "Day" : [],
                    "Hour" : [],
                    "ShapeId" : [],
                    "VegetationType" : [],
                    "Area" : []
                }
                for band, coef in self.__in_environment["caribu opt"].items() : 
                    dico_shape[band+" Eabs"] = []
                    dico_shape[band+" Ei"] = []
                self.__shape_outputs = pandas.DataFrame(dico_shape)
                
                dico_tr = {
                    "Day" : [],
                    "Hour" : [],
                    "primitive_index" : [],
                    "ShapeId" : [],
                    "VegetationType" : [],
                    "Area" : []
                }
                for band, coef in self.__in_environment["caribu opt"].items() : 
                    dico_tr[band+" Eabs"] = []
                    dico_tr[band+" Ei"] = []
                self.__outputs = pandas.DataFrame(dico_tr)

        else:
            raise ValueError("Unknown lightmodel : can be either 'ratp' or 'caribu' ")

    def _create_caribu_legume_sensors(self):
        # récupère les dimensions de la grille
        dxyz = self.__in_lightmodel_parameters["sensors"][1]
        nxyz = self.__in_lightmodel_parameters["sensors"][2] 
        orig = self.__in_lightmodel_parameters["sensors"][3] 

        # réduction si le nombre de couches remplies < nombre de couches prévues
        skylayer = (self.__pmax[2]) // dxyz[2]
        if skylayer < nxyz[2] : skylayer = int(nxyz[2] - 1 - skylayer)
        
        # autrement on garde le nombre de voxels prévus
        else : skylayer = 0 

        # scene plantGL qui acceuillera les capteurs
        s_capt = pgl.Scene()
        
        points = [(0, 0, 0), (dxyz[0], 0, 0), (dxyz[0], dxyz[1], 0), (0, dxyz[1], 0)]  # capeur bas oriente vers le haut
        normals = [(0, 0, 1) for i in range(4)]
        indices = [(0, 1, 2, 3)]

        # carré taille d'un voxel
        carre = pgl.QuadSet(points, indices, normals, indices)
        
        ID_capt = 0
        dico_translat = {}
        for ix in range(nxyz[0]):
            for iy in range(nxyz[1]):
                for iz in range(nxyz[2] - skylayer):
                    # vecteur de translation
                    tx = orig[0] + ix * dxyz[0]
                    ty = orig[1] + iy * dxyz[1]
                    tz = orig[2] + iz * dxyz[2]

                    # retient la translation
                    dico_translat[ID_capt] = [tx, ty, tz]

                    # ajoute un voxel à la scene des capteurs
                    Vox = pgl.Translated(geometry=carre, translation=(tx, ty, tz))
                    s_capt.add(pgl.Shape(geometry=Vox, id=ID_capt))
                    ID_capt += 1

        
        Dico_Sensors = {}
        liste_hmax_capt = []
        liste_x_capt = []
        liste_y_capt = []
        # mise en forme de triangulation CARIBU
        for x in s_capt:

            # pgl_to_caribu(x)

            # Preparation capteurs
            pt_lst = x.geometry.geometry.pointList
            idx_lst = x.geometry.geometry.indexList

            for i in range(0, len(idx_lst)):
                x11 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
                y11 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
                z11 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

                liste_hmax_capt.append(z11)
                liste_x_capt.append(x11)
                liste_y_capt.append(y11)

                x12 = pt_lst[idx_lst[i][1]][0] + dico_translat[x.id][0]
                y12 = pt_lst[idx_lst[i][1]][1] + dico_translat[x.id][1]
                z12 = pt_lst[idx_lst[i][1]][2] + dico_translat[x.id][2]

                liste_hmax_capt.append(z12)
                liste_x_capt.append(x12)
                liste_y_capt.append(y12)

                x13 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
                y13 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
                z13 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

                liste_hmax_capt.append(z13)
                liste_x_capt.append(x13)
                liste_y_capt.append(y13)

                tple1 = [(x11, y11, z11), (x12, y12, z12), (x13, y13, z13)]
                triangle = []
                triangle.append(tple1)

                x21 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
                y21 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
                z21 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

                liste_hmax_capt.append(z21)
                liste_x_capt.append(x21)
                liste_y_capt.append(y21)

                x22 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
                y22 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
                z22 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

                liste_hmax_capt.append(z22)
                liste_x_capt.append(x22)
                liste_y_capt.append(y22)

                x23 = pt_lst[idx_lst[i][3]][0] + dico_translat[x.id][0]
                y23 = pt_lst[idx_lst[i][3]][1] + dico_translat[x.id][1]
                z23 = pt_lst[idx_lst[i][3]][2] + dico_translat[x.id][2]

                liste_hmax_capt.append(z23)
                liste_x_capt.append(x23)
                liste_y_capt.append(y23)

                tple2 = [(x21, y21, z21), (x22, y22, z22), (x23, y23, z23)]
                triangle.append(tple2)

                Dico_Sensors[x.id] = triangle
        
        return Dico_Sensors, s_capt, [(min(liste_x_capt) + max(liste_x_capt)) / 2,
                                                            (min(liste_y_capt) + max(liste_y_capt)) / 2,
                                                            max(liste_hmax_capt)]

    
    def _print_sun(self, day, hour, pyratp, truesolartime):
        """Méthode qui imprime les coordonnées du soleil (pour le débuggage)
            séparée pour la lisibilité dans run(...)
        """
        print("---\t SUN COORDONATES\t ---")
        print("--- Convention x+ = North, vector from sky to floor")
        print("--- azimut: south clockwise E = -90° W = 90°\t zenith: zenith = 0° horizon = 90°")
        print("--- day: %i \t hour: %i \t latitude: %0.2f °"%(day, hour, self.__in_environment["coordinates"][0]))
        
        if truesolartime : 
            caribuhour=hour
            print("--- true solar time")
        
        else :
            om =0.017202*(day-3.244)
            teta=om+0.03344*math.sin(om)*(1+0.021*math.cos(om))-1.3526
            tphi=0.91747*math.sin(teta)/math.cos(teta)
            dphi=math.atan(tphi)-om+1.3526
            if dphi+1. <= 0: dphi=(dphi+1+1000.*math.pi % math.pi)-1.
            eqntime=dphi*229.2
            caribuhour =hour+self.__in_environment["coordinates"][2]+self.__in_environment["coordinates"][1]/15.-eqntime/60. % 24. 
            print("--- local time, true solar time is %0.2f"%(caribuhour))
        
        # CARIBU et RATP sont dans la même convention de coordonnées, x+ = North
        print("--- RATP ---")
        az, ele = 5.,9. # variables fantômes (non récupérées)
        pyratp.shortwave_balance.sundirection(ele, az, self.__in_environment["coordinates"][0], self.__in_environment["coordinates"][1], self.__in_environment["coordinates"][2], day, hour, truesolartime)        
        degtorad = math.pi/180
        azrad = pyratp.shortwave_balance.azdeg
        if math.isnan(azrad) and hour==12. : 
            if self.__in_environment["coordinates"][0] >= 21.11:
                azrad = 0.
            else:
                azrad = 180.
        azrad = -azrad * degtorad
        elerad = pyratp.shortwave_balance.hdeg * degtorad

        sunx = suny = math.cos(elerad)
        sunx *= math.cos(azrad)
        suny *= math.sin(azrad)
        sunz = -math.sin(elerad)

        print("\t azimut: %0.3f \t zenith: %0.3f" % (-azrad*180/math.pi, pyratp.shortwave_balance.hdeg))
        print("\t x: %0.3f \t y: %0.3f \t z: %0.3f" % (sunx, suny, sunz))
        
        print("--- CARIBU ---")
        suncaribu = Sun.Sun()
        suncaribu._set_pos_astro(day, caribuhour, self.__in_environment["coordinates"][0]*math.pi/180)
        sun_str_split = suncaribu.toLight().split(' ')

        print("\t azimut: %0.3f \t zenith: %0.3f" % (-suncaribu.azim*180/math.pi, 90-(suncaribu.elev*180/math.pi)))
        print("\t x: %0.3f \t y: %0.3f \t z: %0.3f" % (float(sun_str_split[1]),float(sun_str_split[2]),float(sun_str_split[3])))
        print("\n")

    
    @property
    def legume_transmitted_light(self):
        return self.__legume_transmitted_light 
    
    @property
    def legume_intercepted_light(self):
        return self.__legume_intercepted_light 

    @property
    def shapes_outputs(self):
        return self.__shape_outputs 

    @property
    def triangles_outputs(self):
        return self.__outputs 
    
    @property
    def voxels_outputs(self):
        return self.__voxels_outputs 

    @property
    def sensors_outputs(self):
        return self.__sensors_outputs 

    @property
    def sun(self):
        return self.__sun

    # aire du plus grand triangle à partir des scènes en entrée
    # -> avant la tesselation
    @property
    def maxtrianglearea(self):
        return self.__maxtrarea
  
    @property
    def legume_empty_layers(self):
        return self.__legume_nb0


    @property
    def tesselationtime(self):
        try:
            return self.__tess_time
        except AttributeError:
            return 0.
    
    @property
    def modelruntime(self):
        try:
            return self.__time_runmodel
        except AttributeError:
            return 0.

    @property
    def triangleLmax(self):
        return self.__triangleLmax

    def RATP_info(self, filename="", writefile=False):
        if self.__lightmodel == "ratp":
            dict_global = {
                "Origin" : [float(self.__ratp_scene.xorig), float(self.__ratp_scene.yorig), float(self.__ratp_scene.zorig) ],
                "Global" : self.__ratp_distrib["global"]
            }

            VegetationType = []
            VoxelId = []
            t_nx = []
            t_ny = []
            t_nz = []
            area = []

            for k in range(self.__ratp_scene.nveg):
                for j in range(self.__ratp_scene.nje[k]):
                    VegetationType.append(j)
                    VoxelId.append(k)
                    t_nx.append(self.__ratp_scene.numx[k])
                    t_ny.append(self.__ratp_scene.numy[k])
                    t_nz.append(self.__ratp_scene.numz[k])
                    area.append(self.__ratp_scene.s_vt_vx[j,k])


            dfvox =  pandas.DataFrame({'VegetationType':VegetationType,
                                        'VoxelId':VoxelId,
                                        'Nx': t_nx,
                                        'Ny': t_ny,
                                        'Nz': t_nz,
                                        'Area': area
                                    })
            
            if self.__in_lightmodel_parameters["angle distrib algo"] == "compute voxel" :
                nbangles = len(self.__ratp_distrib["voxel"][0][0])
                
                # ajoute nb colonnes
                VegetationType = []
                VoxelId = []
                dict_angles={}
                for k in range(self.__ratp_scene.nveg):
                    for j in range(self.__ratp_scene.nje[k]):
                        VegetationType.append(j)
                        VoxelId.append(k)
                        for i in range(nbangles):
                            if i in dict_angles:
                                dict_angles[i].append(self.__ratp_distrib["voxel"][k][j][i])
                            else:
                                dict_angles[i] = [self.__ratp_distrib["voxel"][k][j][i]]
                dict_add={}
                dict_add['VegetationType'] = VegetationType
                dict_add['VoxelId'] = VoxelId
                for key, value in dict_angles.items():
                    title = "C"+str(key)
                    dict_add[title] = value

                dfvox = pandas.merge(dfvox, pandas.DataFrame(dict_add))  

            if writefile:
                if filename=="" : name= "RATP_info.data"
                else: name = filename+"_RATP_info.data"

                f=open(name, 'w')
                f.write("Xorigin\tYorigin\tZorigin\n")
                f.write("%f\t%f\t%f\n" % (dict_global["Origin"][0], dict_global["Origin"][1], dict_global["Origin"][2]))
                f.write("Global Foliar Angles Distribution For Each Entity\n")
                for n in range(len(dict_global["Global"])):
                    f.write("Entity %i\n" % (n))
                    s=""
                    for i in range(len(dict_global["Global"][n])):
                        s=s+str(dict_global["Global"][n][i])+'\t'
                    s+="\n"
                    f.write(s)
                f.write("\n")
                f.close()

                dfvox.to_csv(name, mode='a', index=False, header=True)
                

            return dict_global, dfvox

    def PAR_update_MTG(self, mtg):
        # crée un tableau comme dans caribu_facade de fspm-wheat
        dico_par = {}
        para_dic = {}
        erel_dic = {}
        for s in self.__shape_outputs["ShapeId"]:
            para_dic[s] = self.__shape_outputs[self.__shape_outputs.ShapeId==s]["PARa"].values[0]
            erel_dic[s] = self.__shape_outputs[self.__shape_outputs["ShapeId"].values==s]["xintav"].values[0]
        
        dico_par["PARa"] = para_dic
        dico_par["Erel"] = erel_dic

        # add the missing property
        for param in dico_par:
            if param not in mtg.properties():
                mtg.add_property(param)

            # update the MTG
            mtg.property(param).update(dico_par[param])

    def to_l_egume(self, m_lais = [], energy = 1, list_lstring = [], list_dicFeuilBilanR = [], list_invar = []) :
        if self.__lightmodel == "ratp" :
            epsilon = 1e-14

            # transfert des sorties
            res_abs_i = np.zeros((m_lais.shape[0], m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))
            # si le voxel est vide, on considère le transmis comme ce qui sort d'une de ses surfaces
            dS = self.__ratp_scene.dx * self.__ratp_scene.dy
            res_trans = np.ones((m_lais.shape[1], m_lais.shape[2], m_lais.shape[3])) * (energy * dS)
                    
            for ix in range(m_lais.shape[3]):
                for iy in range(m_lais.shape[2]):
                    for iz in range(self.__ratp_scene.njz):
                        legume_iz = iz + self.__legume_nb0

                        vox_data = self.__voxels_outputs[(self.__voxels_outputs.Nx==m_lais.shape[2] - iy) & 
                                                                        (self.__voxels_outputs.Ny==ix+1) & 
                                                                        (self.__voxels_outputs.Nz==iz+1)]
                        
                        res_trans[legume_iz, iy, ix] = energy * min(sum(vox_data["transmitted"]), dS)
                        # if sum(vox_data["transmitted"]) > dS : print("Warning : transmitted energy > dx.dy")

                        s_entity = 0
                        for k in range(m_lais.shape[0]) : s_entity+=m_lais[k][legume_iz][iy][ix]
                        if s_entity > 0. :
                            for ie in range(m_lais.shape[0]) :
                                if len(vox_data) > 0 :
                                    if vox_data[vox_data.VegetationType == ie+1]["intercepted"].values[0] > epsilon :
                                        res_abs_i[ie, legume_iz, iy, ix] = energy * vox_data[vox_data.VegetationType == ie+1]["intercepted"].values[0]
                                    # si le voxel est non vide, on fixe quand même une valeur min
                                    else:
                                        res_abs_i[ie, legume_iz, iy, ix] = epsilon
                        
            return res_trans, res_abs_i

        elif self.__lightmodel == "caribu" :
            # réécriture de calc_paraF dans ShootMorpho.py pour correspondre aux triangles
            # parip : somme des organes sur tous les statuts
            # parap : idem sauf si plante senescente alors == 0

            # récupère les données de la grille de capteurs
            dxyz = self.__in_lightmodel_parameters["sensors"][1]
            nxyz = self.__in_lightmodel_parameters["sensors"][2] 
            # surface d'une face d'un voxel
            dS = dxyz[0] * dxyz[1]

            # Calcul du rayonnement absorbé par plante et par espèce
            for k in range(len(list_invar)) :
                # on initialise la somme sur chaque plante à 0
                list_invar[k]['parap'] = scipy.array([0] * len(list_invar[k]['Hplante']))
                list_invar[k]['parip'] = scipy.array([0] * len(list_invar[k]['Hplante']))

                if self.__matching_ids :
                    for i in range(len(self.__shape_outputs)) :
                        organe_id = int(self.__shape_outputs.iloc[i]["ShapeId"])

                        # PAR en W/m²
                        par_intercept = self.__shape_outputs.iloc[i]['par Ei']
                        S_leaf = self.__shape_outputs.iloc[i]['Area']
                        
                        id_plante = list_lstring[k][organe_id][0]

                        S_plante = list_dicFeuilBilanR[k]["surf"][id_plante]

                        list_invar[k]['parip'][id_plante] +=  par_intercept * (S_leaf/S_plante)
                            
                        # on enlève les feuilles senescentes 
                        if list_lstring[k][organe_id][9] != 'sen' :
                            list_invar[k]['parap'][id_plante] +=  par_intercept * (S_leaf/S_plante)

                    # conversion
                    list_invar[k]['parap'] = list_invar[k]['parap'] * (3600*24)/1000000
                    list_invar[k]['parip'] = list_invar[k]['parip'] * (3600*24)/1000000
            
            # calcul de res_trans : rayonnement transmis à travers la grille de voxel
            res_trans = np.ones((m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))

            # si scène non vide
            if self.__matching_ids:
                # traitements des sensors différents si scène infinie ou non
                if self.__in_environment["infinite"] : 
                    print("ERROR : Doesn't work yet")
                
                else :
                    # réduction si le nombre de couches remplies < nombre de couches prévues
                    skylayer = (self.__pmax[2]) // dxyz[2]
                    if skylayer < nxyz[2] : skylayer = int(nxyz[2] - 1 - skylayer)
                    # autrement on garde le nombre de voxels prévus
                    else : skylayer = 0 
                
                    ID_capt = 0
                    for ix in range(nxyz[0]):
                        for iy in range(nxyz[1]):
                            for iz in range(nxyz[2] - skylayer):
                                res_trans[(nxyz[2]-1) - iz][iy][ix] -= self.__sensors_outputs['par'][ID_capt]
                                ID_capt += 1
            res_trans = res_trans * energy * dS

            return res_trans
            



    def VTKinit(self, path, plantnames=[], planttrianglevalues=[], printtriangles=True):
        '''construit des fichiers VTK de la triangulation et de la grille de voxels après leur construction

        Args :
            self.__my_scene : liste de Triangle3
            self.__ratp_scene : grille de voxels RATP

            path : string, chemin pour l'écriture des fichiers
            plantnames : liste de string, si l'on veut mettre un nom personnalisé pour chaque entité
            planttrianglevalues : liste de liste de float, si l'on veut mettre des grandeurs associées à chaque triangle pour chaque entité
        '''
        # obligé de forcer ces variables...
        plantnames=[]
        planttrianglevalues=[]
        
        if self.__lightmodel == "ratp" :
            nent = int(self.__ratp_scene.nent)
            # plot dans VTK
            temp1, temp2, temp3 = [], [], []
            # éviter les éléments en trop
            for i in range(self.__ratp_scene.nveg)  :
                for j in range(self.__ratp_scene.nje[i]):
                    temp1.append(self.__ratp_scene.leafareadensity[j, i])
                    # temp1.append(self.__ratp_scene.s_vt_vx[j, i])
                    temp2.append(self.__ratp_scene.nume[j,i])
                    temp3.append(int(i)+1) # kxyz sort en fortran
            lad = [np.array(temp1), np.array(temp2), np.array(temp3)]

            RATP2VTK.RATPVOXELS2VTK(self.__ratp_scene, lad, "LAD", path+"init_voxels.vtk")
        if self.__matching_ids and printtriangles:
            nent = len(self.__in_geometry["scenes"])
            if not plantnames:
                for i in range(nent):
                    plantnames.append("plant_"+str(i))
            
            # pour chaque plante on a une valeur par triangle
            if not planttrianglevalues:
                for i in range(nent):
                    planttrianglevalues.append([])
                
                for tr in self.__my_scene :
                    planttrianglevalues[self.__matching_ids[tr.id][1]].append(10)
                    for i in range(nent):
                        if i != self.__matching_ids[tr.id][1] : planttrianglevalues[i].append(0)

            VTKtriangles(self.__my_scene, planttrianglevalues, plantnames, path+"init_triangles.vtk")

    def VTKout(self, path, iteration=None, triangles=True, voxels=False, outvariables = None):
        '''construit des fichiers VTK de la triangulation avec les valeurs de PAR associées

        Args :
            self.__my_scene : liste de Triangle3
            self.__outputs : dataframe des résultats

            path : string, chemin pour l'écriture des fichiers
            iteration : None ou int, numéro de l'itération temporelle
            voxels : boolean, si ratp écrit ou non le PAR sur les voxels
        '''    
        values = []
        valuesnames = []

        if self.__matching_ids and triangles:
            # si la lumière a été calculée sur plusieurs itération
            if iteration is None:
                for ite in range(int(max(self.__outputs['Iteration']))):
                    par = []
                    idtr=0
                    for tr in self.__my_scene:
                        df = self.__outputs[(self.__outputs.Iteration == ite+1) & (self.__outputs.primitive_index == idtr)]
                        voxpar = df['PARa'].values[0]
                        
                        par.append(voxpar)
                        idtr += 1

                    VTKtriangles(self.__my_scene, [par], ['PARa'], path+"triangles_PAR_"+str(ite)+".vtk")
            
            # ou alors sur une itération en particulier
            else:
                if self.__lightmodel == "ratp" :
                    par = []
                    for i,tr in enumerate(self.__my_scene):
                        df = self.__outputs[(self.__outputs.primitive_index == i)]
                        par.append(df['PARa'].values[0])                
                    values = [par]
                    valuesnames=["PARa"]
    
                elif self.__lightmodel == "caribu" :
                    for band, coef in self.__in_environment["caribu opt"].items() : 
                        bandval = []
                        valuesnames.append(band+"_Ei")
                        for i,tr in enumerate(self.__my_scene):
                            df = self.__outputs[(self.__outputs.primitive_index == i)]
                            bandval.append(df[band+" Ei"].values[0])   
                        values.append(bandval)
                VTKtriangles(self.__my_scene, values, valuesnames, path+"triangles_PAR_"+str(iteration)+".vtk")

        # VTK des voxels
        if self.__lightmodel == "ratp" and voxels:
            if outvariables == None :
                # plot dans VTK
                temp1, temp2, temp3 = [], [], []
                # éviter les éléments en trop
                for i in range(self.__ratp_scene.nveg)  :
                    for j in range(self.__ratp_scene.nje[i]):
                        dfvox = self.__voxels_outputs[(self.__voxels_outputs.VoxelId==i+1) & (self.__voxels_outputs.VegetationType==j+1) & (self.__voxels_outputs.Iteration==1)]
                        temp1.append(dfvox["PARa"].values[0])
                        temp2.append(self.__ratp_scene.nume[j,i])
                        temp3.append(int(i)+1) # kxyz sort en fortran
                para = [np.array(temp1), np.array(temp2), np.array(temp3)]

                RATP2VTK.RATPVOXELS2VTK(self.__ratp_scene, para, "PARa", path+"PARa_voxels_"+str(iteration)+".vtk")
            else:
                for out in outvariables :
                    temp1, temp2, temp3 = [], [], []
                    # éviter les éléments en trop
                    for i in range(self.__ratp_scene.nveg)  :
                        for j in range(self.__ratp_scene.nje[i]):
                            dfvox = self.__voxels_outputs[(self.__voxels_outputs.VoxelId==i+1) & (self.__voxels_outputs.VegetationType==j+1) & (self.__voxels_outputs.Iteration==1)]
                            temp1.append(dfvox[out].values[0])
                            temp2.append(self.__ratp_scene.nume[j,i])
                            temp3.append(int(i)+1) # kxyz sort en fortran
                    para = [np.array(temp1), np.array(temp2), np.array(temp3)]

                    RATP2VTK.RATPVOXELS2VTK(self.__ratp_scene, para, out, path+out+"_voxels_"+str(iteration)+".vtk")


    @staticmethod
    def PlantGL_to_VTK(scenes, ite, path, in_unit="m", out_unit="m"):
        units = {'mm': 0.001, 'cm': 0.01, 'dm': 0.1, 'm': 1, 'dam': 10, 'hm': 100,'km': 1000}
        triangleslist=[]
        rescale=False
        if (in_unit != out_unit) : rescale=True

        if type(scenes) == list :
            for s in scenes :   
                for id, pgl_objects in s.todict().items():           
                    tri_list = list(itertools.chain(*[pgl_to_triangles(pgl_object) for pgl_object in pgl_objects]))
                    for tr in tri_list:
                        if rescale : tr.rescale(units[in_unit]/units[out_unit])
                        tr.set_id(0)
                    triangleslist.extend(tri_list)
        else :
            for id, pgl_objects in scenes.todict().items():           
                tri_list = list(itertools.chain(*[pgl_to_triangles(pgl_object) for pgl_object in pgl_objects]))
                for tr in tri_list:
                    if rescale : tr.rescale(units[in_unit]/units[out_unit])
                    tr.set_id(0)
                triangleslist.extend(tri_list)

        VTKtriangles(triangleslist, [], [], path+"triangles_plantgl_"+str(ite)+".vtk")


    def s5(self):
        '''construit les fichiers d'entrée pour s5 et l'exécute
        '''
        s5folder = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), os.path.normpath("s5"))
        fort51 = os.path.join(s5folder, os.path.normpath("fort.51"))
        s5par = os.path.join(s5folder, os.path.normpath("s5.par"))

        # écriture du fichier fort.51 contenant la triangulation
        f=open(fort51, 'w')
        c_tr=1
        for tr in self.__my_scene:
            if (self.__matching_ids[tr.id][0], self.__matching_ids[tr.id][1]) in self.__in_geometry["stems id"]:
                stem='000'
            else:
                stem='001'
            label = str(self.__matching_ids[tr.id][1]+1)+str("%05i"%(self.__matching_ids[tr.id][1]+1))+stem+'000'#+str("%03i"%(c_tr))
            f.write("p\t1\t%s\t3\t"%(label))
            for i in range(3):
                f.write("%f\t%f\t%f\t"%(tr[i][0],tr[i][1],tr[i][2]))
            f.write("\n")
            c_tr += 1
        f.close()

        # écriture du fichier s5.par contenant les informations de la grille
        f=open(s5par, 'w')    
        f.write("%i\t%i\t%i\t%i\t0\n"%(self.__ratp_scene.nent, 9, 9, self.__ratp_scene.njz))
        for i in range(self.__ratp_scene.njz):
            f.write("%f\t"%(self.__ratp_scene.dz[i]))
        f.write("\n")

        njx = self.__ratp_scene.njx
        njy = self.__ratp_scene.njy
        dx = self.__ratp_scene.dx
        dy = self.__ratp_scene.dy
        f.write("%f\t%i\t%f\t%f\t%i\t%f\n"%(njx*dx, njx, dx, njy*dy, njy, dy))
        f.close()

        # exécution de s5 dans un sous process
        subprocess.call(".\s5.exe", shell=True, cwd=s5folder)
        
        print("\n"+"--- Fin de s5.f")

    
    def s2v(self):
        '''construit les fichiers d'entrée pour s2v et l'exécute
        '''
        
        s2vfolder = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), os.path.normpath("s2v"))
        fort51 = os.path.join(s2vfolder, os.path.normpath("fort.51"))
        s2vpar = os.path.join(s2vfolder, os.path.normpath("s2v.par"))

        # ecrit dans le dossier de s2v
        f=open(fort51, 'w')
        c_tr=1
        for tr in self.__my_scene:
            if (self.__matching_ids[tr.id][0], self.__matching_ids[tr.id][1]) in self.__id_stems:
                stem='000'
            else:
                stem='001'
            label = str(self.__matching_ids[tr.id][1]+1)+str("%05i"%(self.__matching_ids[tr.id][1]+1))+stem+'000'#+str("%03i"%(c_tr))
            f.write("p\t1\t%s\t3\t"%(label))
            for i in range(3):
                f.write("%f\t%f\t%f\t"%(tr[i][0],tr[i][1],tr[i][2]))
            f.write("\n")
            c_tr += 1
        f.close()

        f=open(s2vpar, 'w')
        # ligne 1 
        f.write("9 9 %i\n"%(self.__ratp_scene.njz))

        # ligne 2
        for i in range(self.__ratp_scene.njz):
            f.write("%f "%(self.__ratp_scene.dz[i]))
        f.write("\n")
        
        # ligne 3
        njx = self.__ratp_scene.njx
        njy = self.__ratp_scene.njy
        dx = self.__ratp_scene.dx
        dy = self.__ratp_scene.dy
        f.write("%f %i %f %f %i %f %i\n"%(njx*dx, njx, dx, njy*dy, njy, dy, self.__ratp_scene.nent))

        # ligne 4
        f.write("%f %f %f\n"%(self.__ratp_scene.xorig, self.__ratp_scene.yorig, self.__ratp_scene.zorig))

        f.close()

        # exécution de s5 dans un sous process
        subprocess.call(".\s2v.exe", shell=True, cwd=s2vfolder)

        print("--- Fin de s2v.cpp")

    def __str__(self):
        '''imprime les infos de la scène couplée
        '''
        out = "\n-------------------------------------\n"
        out += "Light Vege Manager v0.0\n"
        out += str(len(self.__in_scenes))+" entités\n"
        out += "Modèle de lumière : "+self.__lightmodel
        if self.__lightmodelparam[4]>0 : out += " avec tesselation\n"
        else : out += "\n"
        out += "Scènes en entrée : \n"
        for i in range(len(self.__in_scenes)):
            out+= "\t scene "+str(i)+" : "+str(len(self.__in_scenes[i]))+" shapes, modèle : "+str(self.__in_names[i])+"\n"

        out += "\n-------------------------------------\n"
        return out


    def _discrete_sky(self, n_azimuts, n_zeniths, sky_type) :
        ele=[]
        azi=[]
        da = 2 * math.pi / n_azimuts
        dz = math.pi / 2 / n_zeniths    
        todeg = 180/math.pi          
        for j in range(n_azimuts):
            for k in range(n_zeniths):
                azi.append((j * da + da / 2)*todeg)
                ele.append((k * dz + dz / 2)*todeg)
        n = n_azimuts*n_zeniths
        
        omega=[2*math.pi/n]*n
        
        def uoc (teta, dt, dp):
            """ teta: angle zenithal; phi: angle azimutal du soleil """
            dt /= 2.
            x = math.cos(teta-dt)
            y = math.cos(teta+dt)
            E = (x*x-y*y)*dp/2./math.pi
            return E
        def soc (teta, dt, dp):
            """ teta: angle zenithal; phi: angle azimutal du soleil """
            dt /= 2.
            x = math.cos(teta-dt)
            y = math.cos(teta+dt)
            E = (3/14.*(x*x-y*y) + 6/21.*(x*x*x-y*y*y))*dp/math.pi
            return E
        pc=[]
        for j in range(n_azimuts):
            for k in range(n_zeniths):
                azim,elv = j * da + da / 2, k * dz + dz / 2
                I=0
                if(sky_type=='soc'):
                    I = soc(elv, dz, da)
                elif(sky_type=='uoc'):
                    I = uoc(elv, dz, da)

                pc.append(I) 
        return ele, azi, omega, pc