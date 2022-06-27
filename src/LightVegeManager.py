import openalea.plantgl.all as pgl

from PyRATP.pyratp import grid
from PyRATP.pyratp import RATP2VTK
from PyRATP.pyratp.vegetation import Vegetation
from PyRATP.pyratp.skyvault import Skyvault
from PyRATP.pyratp.micrometeo import MicroMeteo
from PyRATP.pyratp.runratp import runRATP

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun, spitters_horaire, turtle, Sun

import os, subprocess
import itertools
import numpy as np
import pandas
import concurrent.futures
import time

from src.Polygons import *
from src.MyTesseletor import *

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
                    environment = {},
                    lightmodel="ratp",
                    lightmodel_parameters = {},
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
                    ele=[]
                    azi=[]
                    da = 2 * math.pi / self.__in_environment["sky"][0]
                    dz = math.pi / 2 / self.__in_environment["sky"][1]      
                    todeg = 180/math.pi          
                    for j in range(self.__in_environment["sky"][0]):
                        for k in range(self.__in_environment["sky"][1]):
                            azi.append((j * da + da / 2)*todeg)
                            ele.append((k * dz + dz / 2)*todeg)
                    n = self.__in_environment["sky"][0]*self.__in_environment["sky"][1]
                    
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
                    for j in range(self.__in_environment["sky"][0]):
                        for k in range (self.__in_environment["sky"][1]):
                            azim,elv = j * da + da / 2, k * dz + dz / 2
                            I=0
                            if(self.__in_environment["sky"][2]=='soc'):
                                I = soc(elv, dz, da)
                            else:
                                I = uoc(elv, dz, da)

                            pc.append(I) 
                    
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
        
        # paramètre par défaut des tiges
        if "stems id" not in self.__in_geometry : self.__in_geometry["stems id"] = []

        # couplage des scènes dans un plantGL commun et création du tableau des ids
        # [plantGL] -(+)-> [Triangle3] -> transformation sur les Triangle3 -> [Triangle3]
        # une scène plantgl par espèce
        self.__matching_ids = {}
        self.__my_scene=[]
        count=0
        for i_esp, scene in enumerate(self.__in_geometry["scenes"]) :
            if isinstance(scene, pgl.Scene):
                for id, pgl_objects in scene.todict().items():
                    lastid = len(self.__my_scene)
                    
                    tri_list = list(itertools.chain(*[pgl_to_triangles(pgl_object) for pgl_object in pgl_objects]))
                    # on set l'id des triangles de la shape
                    for tr in tri_list:
                        tr.set_id(count)
                    self.__my_scene.extend(tri_list)
                    # on set le tableau des indices
                    self.__matching_ids[count] = (id, i_esp, list(range(lastid,lastid+len(tri_list))))
                    count += 1
        
        # applique les transformations sur les triangles
        if "transformations" in self.__in_geometry :
            for i_esp in range(len(self.__in_geometry["scenes"])):
                if isinstance(self.__in_geometry["scenes"][i_esp], pgl.Scene):
                    for tr in self.__my_scene:
                        if self.__matching_ids[tr.id][1] == i_esp:
                            if "rescale" in  self.__in_geometry["transformations"] : tr.rescale(self.__in_geometry["transformations"]["rescale"][i_esp])
                            
                            if "translate" in  self.__in_geometry["transformations"]  : tr.translate(self.__in_geometry["transformations"]["translate"][i_esp])
                            
                            if "scenes unit" in  self.__in_geometry["transformations"] :
                                # recupère la variable pour la lisibilité
                                scene_unit = self.__in_geometry["transformations"]["scenes unit"][i_esp]
                                
                                if (scene_unit != self.__main_unit) and (scene_unit in self.units):
                                    tr.rescale(self.units[scene_unit]/self.units[self.__main_unit])
                            
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
        
        # enregistre l'aire du plus grand triangle (indicateur par rapport au besoin de tesselation)
        self.__maxtrarea = 0.
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
        for tr in self.__my_scene:
            for i in range(3) :
                p = tr[i]
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

        # RATP
        if self.__lightmodel == "ratp":
            # on sépare les tiges dans une nouvelle entité si il y a des shapes non tiges
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

                    if stem[1] not in mem_mu:
                        self.__in_lightmodel_parameters["mu"].append(self.__in_lightmodel_parameters["mu"][stem[1]])
                        self.__in_environment["reflectance coefficients"].append(self.__in_environment["reflectance coefficients"][stem[1]])
                        mem_mu.append(stem[1])

            ## distribution d'angles pour chaque entité
            distrib = []

            # recherche du nb d'entité
            nent=0
            for key, val in self.__matching_ids.items():
                if val[1]+1 > nent:
                    nent = val[1]+1

            # calcul en dynamique
            # ele_option = nombre de classes
            # on fait le calcul de la distribution global avant la tesselation des triangles
            # pour optimiser les calculs
            if self.__in_lightmodel_parameters["angle distrib algo"] == "compute global" :
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

                    distrib.append(classes)

                # convertit en pourcentage
                for n in range(nent):
                    for i in range(len(distrib[n])):
                        distrib[n][i] *= 1/t_nent_area[n]
            
            # lecture du fichier
            # ele_option = chemin du fichier
            elif self.__in_lightmodel_parameters["angle distrib algo"] == "file" : 
                f_angle = open(self.__in_lightmodel_parameters["angle distrib file"], 'r')
                for i in range(len(self.__in_scenes)):
                    line = f_angle.readline()
                    distrib.append([float(x) for x in line.split(',')[1:]])

            # on ajuste au besoin les min-max si la scène est plane pour avoir un espace 3D
            if xmin == xmax:
                xmax += self.__in_lightmodel_parameters["voxel size"][0]
                xmin -= self.__in_lightmodel_parameters["voxel size"][0]
            if ymin == ymax:
                ymax += self.__in_lightmodel_parameters["voxel size"][1]
                ymin -= self.__in_lightmodel_parameters["voxel size"][1]
            if zmin == zmax:
                zmax += self.__in_lightmodel_parameters["voxel size"][2]
                zmin -= self.__in_lightmodel_parameters["voxel size"][2]
            
            self.__pmax = Vector3(xmax, ymax, zmax)
            self.__pmin = Vector3(xmin, ymin, zmin)

            # nombre de voxels
            nx = int((self.__pmax[0] - self.__pmin[0]) // self.__in_lightmodel_parameters["voxel size"][0])
            ny = int((self.__pmax[1] - self.__pmin[1]) // self.__in_lightmodel_parameters["voxel size"][1])
            nz = int((self.__pmax[2] - self.__pmin[2]) // self.__in_lightmodel_parameters["voxel size"][2])
            if (self.__pmax[0] - self.__pmin[0]) % self.__in_lightmodel_parameters["voxel size"][0] > 0 : nx += 1
            if (self.__pmax[1] - self.__pmin[1]) % self.__in_lightmodel_parameters["voxel size"][1] > 0 : ny += 1
            if (self.__pmax[2] - self.__pmin[2]) % self.__in_lightmodel_parameters["voxel size"][2] > 0 : nz += 1
            
            # définit une origine en Pmin 
            xorig, yorig, zorig = self.__pmin[0], self.__pmin[1], -self.__pmin[2]

            # création de la grille
            # si pas de rayonnement réfléchi on annules la réflexion du sol
            if self.__in_environment["reflected"] :
                mygrid = grid.Grid.initialise(nx, ny, nz, 
                                                self.__in_lightmodel_parameters["voxel size"][0], 
                                                self.__in_lightmodel_parameters["voxel size"][1], 
                                                self.__in_lightmodel_parameters["voxel size"][2], 
                                                xorig, yorig, zorig, 
                                                self.__in_environment["coordinates"][0], self.__in_environment["coordinates"][1], self.__in_environment["coordinates"][2], 
                                                nent, 
                                                self.__in_lightmodel_parameters["soil reflectance"], 
                                                toric=self.__in_environment["infinite"])
            else :
                mygrid = grid.Grid.initialise(nx, ny, nz, 
                                                self.__in_lightmodel_parameters["voxel size"][0], 
                                                self.__in_lightmodel_parameters["voxel size"][1], 
                                                self.__in_lightmodel_parameters["voxel size"][2], 
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
            self.__ratp_scene = mygrid
            self.__tr_vox = matching 

            # si la distribution d'angle est par voxel, on est
            # obligé de prendre en compte les triangles après
            # la tesselation
            if self.__in_lightmodel_parameters["angle distrib algo"] == "compute voxel":
                angles = list(np.linspace(90/self.__in_lightmodel_parameters["nb angle classes"], 91, self.__in_lightmodel_parameters["nb angle classes"]))
                t_area=[]
                # pour chaque voxel
                for k in range(self.__ratp_scene.nveg):
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
                    distrib.append(distrib_nent)

                to_remove=[]
                # convertit en pourcentage
                for k in range(self.__ratp_scene.nveg):
                    for n in range(nent):
                        # peut survenir si une entité n'est pas dans le voxel
                        if t_area[k][n] != 0 :
                            for i in range(self.__in_lightmodel_parameters["nb angle classes"]):
                                distrib[k][n][i] *= 1/t_area[k][n]
                        
                        # enleve les entités en trop
                        else:
                            to_remove.append((k,n))
                for t in to_remove:
                    del distrib[t[0]][t[1]]

            self.__ratp_distrib = distrib
            
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

    def run(self, meteo_path="", PARi=0, day=0, hour=0, parunit="micromol.m-2.s-1", truesolartime=False, printsun=False):
        '''calcul du bilan radiatif
        Args :
            meteo_path : string, chemin du fichier meteo
            PARi : float, PARi en µmol.m-2.s-1 ou W.m-2
            day : float, day of the year 
            hour : float, hour TU
            parunit : "micromol.m-2.s-1" ou "W.m-2"
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
                                                'distinc' : self.__ratp_distrib[id],
                                                'rf' : self.__in_environment["reflectance coefficients"][id]
                                                })
                    else :
                        entities_param.append({
                                                'mu' : mu_ent,
                                                'distinc' : self.__ratp_distrib[id],
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
                vegetation = Vegetation.initialise(entities_param, pervoxel=True, distribvox=self.__ratp_distrib)

            # init météo, PARi en entrée en W.m-2
            if meteo_path == "":
                if parunit == "micromol.m-2.s-1" :
                    #: Spitters's model estimating for the diffuse:direct ratio
                    # coefficient 2.02 : 4.6 (conversion en W.m-2) x 0.439 (PAR -> global)
                    RdRs = spitters_horaire.RdRsH(Rg=PARi/2.02, DOY=day, heureTU=hour, latitude=self.__in_environment["coordinates"][0])
                    
                    # coeff 4.6 : https://www.researchgate.net/post/Can-I-convert-PAR-photo-active-radiation-value-of-micro-mole-M2-S-to-Solar-radiation-in-Watt-m2
                    PARi = PARi / 4.6  # W.m-2
                else:
                    RdRs = spitters_horaire.RdRsH(Rg=PARi/0.439, DOY=day, heureTU=hour, latitude=self.__in_environment["coordinates"][0])
                
                # PAR et Dif en W.m^-2
                if self.__in_environment["diffus"] :
                    # Direct et diffus
                    if self.__in_environment["direct"]: 
                        met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=PARi, Rdif=PARi*RdRs, truesolartime=truesolartime)
                    # uniquement diffus
                    else:
                        met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=PARi, Rdif=PARi, truesolartime=truesolartime)
                # uniquement direct
                else :  met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=PARi, Rdif=0, truesolartime=truesolartime)
            
            # lecture d'un fichier
            else:
                met = MicroMeteo.read(meteo_path, truesolartime)

            # impression de la position du soleil
            if printsun: 
                from PyRATP.pyratp import pyratp    
                self._print_sun(day, hour, pyratp, truesolartime)

            start=time.time()
            # Calcul du bilan radiatif sur chaque pas de temps du fichier météo
            res = runRATP.DoIrradiation(self.__ratp_scene, vegetation, self.__sky, met)
            self.__time_runmodel = time.time() - start

            # Mise en forme des sorties
            # création de plusieurs tableaux intermédiaires qui serviront à trier les sorties
            entity = {}
            for id, match in self.__matching_ids.items():
                entity[id] = match[1] + 1
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
            VegetationType,Iteration,day,hour,VoxelId,ShadedPAR,SunlitPAR,ShadedArea,SunlitArea, xintav= res.T

            # ('PAR' is expected in  Watt.m-2 in RATP input, whereas output is in micromol => convert back to W.m2 (cf shortwavebalance, line 306))
            # on reste en micromol !
            # On enregistre tout dans une dataframe pandas
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
                                'ShadedPAR':ShadedPAR, # /4.6,
                                'SunlitPAR':SunlitPAR, # /4.6,
                                'ShadedArea':ShadedArea,
                                'SunlitArea': SunlitArea,
                                'Area': ShadedArea + SunlitArea,
                                'PARa': para_list,
                                'xintav': erel_list, 
                                })
            
            # tri de la dataframe par rapport aux shapes et triangles
            
            # ne prend pas le sol
            dfvox = dfvox[dfvox['VegetationType'] > 0]
            
            # enregistre la dataframe des voxels
            self.__voxels_outputs = dfvox

            # nouvelle data frame avec les triangles en index
            dfmap = pandas.DataFrame({'primitive_index': index,'shape_id': sh_id, 'VoxelId':vox_id, 'VegetationType':[entity[sh_id] for sh_id in sh_id], 'primitive_area':s})

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

        # CARIBU
        elif self.__lightmodel == "caribu":
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

                #: Optical properties
                opt = {'par': {}}
                for id, val in self.__matching_ids.items():
                    # id : id de la shape, val : [id shape en input, id de l'entité]
                    # c'est une tige il n'y a pas de transmission
                    if (val[0], val[1]) in self.__in_geometry["stems id"]:
                        opt['par'][id] = (self.__in_environment["reflectance coefficients"][val[1]][0],)  #: (reflectance,) of the stems ou organe opaque
                    else:
                        opt['par'][id] = (self.__in_environment["reflectance coefficients"][val[1]][0], self.__in_environment["reflectance coefficients"][val[1]][1]) #: (reflectance, transmittance) of the adaxial side of the leaves, élément translucide symétrique

                # construction d'une scène ciel et soleil
                if self.__in_environment["infinite"] : # on ajoute un domaine pour la création du pattern
                    c_scene_sky = CaribuScene(scene=self.__caribu_scene, light=self.__sky, opt=opt, scene_unit=self.__main_unit, pattern=self.__in_geometry["domain"])
                    c_scene_sun = CaribuScene(scene=self.__caribu_scene, light=self.__sun, opt=opt, scene_unit=self.__main_unit, pattern=self.__in_geometry["domain"])
                else:
                    c_scene_sky = CaribuScene(scene=self.__caribu_scene, light=self.__sky, opt=opt, scene_unit=self.__main_unit)
                    c_scene_sun = CaribuScene(scene=self.__caribu_scene, light=self.__sun, opt=opt, scene_unit=self.__main_unit)

                if self.__in_environment["diffus"] :
                    # Direct et diffus
                    if self.__in_environment["direct"]: 
                        sun_sky_option = "mix"
                    else:
                        sun_sky_option = "sky"
                else:
                    sun_sky_option = "sun"
                
                # choix du type de rayonnement
                # direct=False : active la rediffusion
                # infinite=False : désactive la répétition infinie de la scène (pas de pattern défini ici)
                if sun_sky_option == "mix":
                    start=time.time()
                    raw_sun, aggregated_sun = c_scene_sun.run(direct=True, infinite=self.__in_environment["infinite"])
                    Erel_sun = aggregated_sun['par']['Eabs']
                    Ei_sun = aggregated_sun['par']['Ei']
                    raw_sky, aggregated_sky = c_scene_sky.run(direct=True, infinite=self.__in_environment["infinite"])
                    Erel_sky = aggregated_sky['par']['Eabs']
                    Ei_sky = aggregated_sky['par']['Ei']
                    self.__time_runmodel = time.time() - start

                    #: Spitters's model estimating for the diffuse:direct ratio
                    # % de sky dans la valeur d'énergie finale
                    Rg = PARi / 2.02  #: Global Radiation (W.m-2)
                    RdRs = spitters_horaire.RdRsH(Rg=Rg, DOY=day, heureTU=hour, latitude=self.__in_environment["coordinates"][0])  #: Diffuse fraction of the global irradiance
                    Erel = {}
                    Ei_output_shape = {}
                    for element_id, Erel_value in Erel_sky.items():
                        Erel[element_id] = RdRs * Erel_value + (1 - RdRs) * Erel_sun[element_id]
                        Ei_output_shape[element_id] = RdRs * Ei_sky[element_id] + (1 - RdRs) * Ei_sun[element_id]

                    Erel_output_shape = Erel
                    PARa_output_shape = {k: v * PARi for k, v in Erel_output_shape.items()}
                    PARi_output_shape = {k: v * PARi for k, v in Ei_output_shape.items()}
                    
                    # dataframe des triangles
                    # même numération que caribu
                    Erel_output_tr = {}
                    Ei_output_tr = {}
                    count=0
                    for key, val in raw_sky['par']['Eabs'].items():
                        for i,par in enumerate(val):
                            Erel_output_tr[count] = RdRs*par+(1-RdRs)*raw_sun['par']['Eabs'][key][i]
                            Ei_output_tr[count] = RdRs*raw_sky['par']['Ei'][key][i]+(1-RdRs)*raw_sun['par']['Ei'][key][i]
                            count+=1
                    PARa_output_tr = {k: v * PARi for k, v in Erel_output_tr.items()}
                    PARi_output_tr = {k: v * PARi for k, v in Ei_output_tr.items()}
                
                elif sun_sky_option == "sun":
                    start=time.time()
                    raw_sun, aggregated_sun = c_scene_sun.run(direct=True, infinite=self.__in_environment["infinite"])
                    self.__time_runmodel = time.time() - start
                    Erel_output_shape = aggregated_sun['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                    Ei_output_shape = aggregated_sun['par']['Ei']
                    PARa_output_shape = {k: v * PARi for k, v in Erel_output_shape.items()}
                    PARi_output_shape = {k: v * PARi for k, v in Ei_output_shape.items()}
                    
                    # dataframe des triangles
                    Erel_output_tr = {}
                    Ei_output_tr = {}
                    count=0
                    for key, val in raw_sun['par']['Eabs'].items():
                        for i,par in enumerate(val):
                            Erel_output_tr[count] = raw_sun['par']['Eabs'][key][i]
                            Ei_output_tr[count] = raw_sun['par']['Ei'][key][i]
                            count+=1
                    PARa_output_tr = {k: v * PARi for k, v in Erel_output_tr.items()}
                    PARi_output_tr = {k: v * PARi for k, v in Ei_output_tr.items()}
                
                elif sun_sky_option == "sky":
                    start=time.time()
                    raw_sky, aggregated_sky = c_scene_sky.run(direct=True, infinite=self.__in_environment["infinite"])
                    self.__time_runmodel = time.time() - start
                    Erel_output_shape = aggregated_sky['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                    Ei_output_shape = aggregated_sky['par']['Ei']
                    PARa_output_shape = {k: v * PARi for k, v in Erel_output_shape.items()}
                    PARi_output_shape = {k: v * PARi for k, v in Ei_output_shape.items()}
                    
                    # dataframe des triangles
                    Erel_output_tr = {}
                    Ei_output_tr = {}
                    count=0
                    for key, val in raw_sky['par']['Eabs'].items():
                        for i,par in enumerate(val):
                            Erel_output_tr[count] = raw_sky['par']['Eabs'][key][i]
                            Ei_output_tr[count] = raw_sky['par']['Ei'][key][i]
                            count+=1
                    PARa_output_tr = {k: v * PARi for k, v in Erel_output_tr.items()}
                    PARi_output_tr = {k: v * PARi for k, v in Ei_output_tr.items()}
                    
                else:
                    raise ValueError("Unknown sun_sky_option : can be either 'mix', 'sun' or 'sky'.")
        
            # enregistre les valeurs par shape et plantes
            s_shapes = [0]*len(self.__matching_ids)
            s_area = [0]*len(self.__matching_ids)
            s_par = [0]*len(self.__matching_ids)
            s_pari = [0]*len(self.__matching_ids)
            s_xintav = [0]*len(self.__matching_ids)
            s_day = [0]*len(self.__matching_ids)
            s_hour = [0]*len(self.__matching_ids)
            s_ent = [0]*len(self.__matching_ids)

            for key,val in self.__matching_ids.items():
                s_shapes[key] = val[0]
                if sun_up: s_par[key] = PARa_output_shape[key]
                if sun_up: s_pari[key] = PARi_output_shape[key]
                if sun_up: s_xintav[key] = Erel_output_shape[key]
                if sun_up: 
                    tr_shape = [self.__my_scene[t] for t in val[2]]
                    s_area[key] = sum([tr.area for tr in tr_shape])
                s_day[key] = day
                s_hour[key] = hour
                s_ent[key] = self.__matching_ids[key][1]

            self.__shape_outputs = pandas.DataFrame({
                "Day" : s_day,
                "Hour" : s_hour,
                "ShapeId" : s_shapes,
                "VegetationType" : s_ent,
                "Area" : s_area,
                "PARi" : s_pari,
                "PARa" : s_par,
                "xintav" : s_xintav
            })

            # enregistre les valeurs par triangles
            s_shapes = [0]*len(self.__my_scene)
            s_tr = [0]*len(self.__my_scene)
            s_area=[0]*len(self.__my_scene)
            s_par=[0]*len(self.__my_scene)
            s_pari=[0]*len(self.__my_scene)
            s_xintav=[0]*len(self.__my_scene)
            s_day=[0]*len(self.__my_scene)
            s_hour=[0]*len(self.__my_scene)
            s_ent=[0]*len(self.__my_scene)

            for id,tr in enumerate(self.__my_scene):
                s_shapes[id] = self.__matching_ids[tr.id][0]
                if sun_up:s_par[id] = PARa_output_tr[id]
                if sun_up:s_pari[id] = PARi_output_tr[id]
                if sun_up:s_xintav[id] = Erel_output_tr[id]
                s_tr[id] = id
                s_day[id] = day
                s_hour[id] = hour
                s_ent[id] = self.__matching_ids[tr.id][1]
                s_area[id] = tr.area               

            self.__outputs = pandas.DataFrame({
                "Day" : s_day,
                "Hour" : s_hour,
                "primitive_index" : s_tr,
                "ShapeId" : s_shapes,
                "VegetationType" : s_ent,
                "Area" : s_area,
                "PARi" : s_pari,
                "PARa" : s_par,
                "xintav" : s_xintav
            })

        else:
            raise ValueError("Unknown lightmodel : can be either 'ratp' or 'caribu' ")

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
    def shapes_outputs(self):
        return self.__shape_outputs 

    @property
    def triangles_outputs(self):
        return self.__outputs 
    
    @property
    def voxels_outputs(self):
        return self.__voxels_outputs 

    @property
    def sun(self):
        return self.__sun

    # aire du plus grand triangle à partir des scènes en entrée
    # -> avant la tesselation
    @property
    def maxtrianglearea(self):
        return self.__maxtrarea

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

    def VTKinit(self, path, plantnames=[], planttrianglevalues=[]):
        '''construit des fichiers VTK de la triangulation et de la grille de voxels après leur construction

        Args :
            self.__my_scene : liste de Triangle3
            self.__ratp_scene : grille de voxels RATP

            path : string, chemin pour l'écriture des fichiers
            plantnames : liste de string, si l'on veut mettre un nom personnalisé pour chaque entité
            planttrianglevalues : liste de liste de float, si l'on veut mettre des grandeurs associées à chaque triangle pour chaque entité
        '''
        nent = len(self.__in_geometry["scenes"])
        if self.__lightmodel == "ratp":
            nent = int(self._LightVegeManager__ratp_scene.nent)
            # plot dans VTK
            temp1, temp2, temp3 = [], [], []
            # éviter les éléments en trop
            for i in range(self.__ratp_scene.nveg)  :
                for j in range(self.__ratp_scene.nje[i]):
                    temp1.append(self.__ratp_scene.leafareadensity[j, i])
                    temp2.append(self.__ratp_scene.nume[j,i])
                    temp3.append(int(i)+1) # kxyz sort en fortran
            lad = [np.array(temp1), np.array(temp2), np.array(temp3)]

            RATP2VTK.RATPVOXELS2VTK(self.__ratp_scene, lad, "LAD", path+"init_voxels.vtk")

        if plantnames==[]:
            for i in range(nent):
                plantnames.append("plant_"+str(i))
        
        # pour chaque plante on a une valeur par triangle
        if planttrianglevalues==[]:
            for i in range(nent):
                planttrianglevalues.append([])

            for tr in self.__my_scene :
                planttrianglevalues[self.__matching_ids[tr.id][1]].append(10)
                for i in range(nent):
                    if i != self.__matching_ids[tr.id][1] : planttrianglevalues[i].append(0)

        VTKtriangles(self.__my_scene, planttrianglevalues, plantnames, path+"init_triangles.vtk")

    def VTKout(self, path, iteration=None, voxels=False):
        '''construit des fichiers VTK de la triangulation avec les valeurs de PAR associées

        Args :
            self.__my_scene : liste de Triangle3
            self.__outputs : dataframe des résultats

            path : string, chemin pour l'écriture des fichiers
            iteration : None ou int, numéro de l'itération temporelle
            voxels : boolean, si ratp écrit ou non le PAR sur les voxels
        '''    
        par = []

        # si la lumière a été calculée sur plusieurs itération
        if iteration is None:
            for ite in range(int(max(self.__outputs['Iteration']))):
                par = []
                idtr=0
                for tr in self.__my_scene:
                    df = self.__outputs[(self.__outputs.Iteration == ite+1) & (self.__outputs.primitive_index == idtr)]
                    voxpar = df['PARa'].values[0]
                    
                    # tentative de réduction du par selon l'aire du triangle mais peut etre déjà fait dans le tri de df
                    #surfvox = mygrid.s_vt_vx[matching_id[id][1], int(d_E2[str(idtr)])]
                    #if surfvox==0 : surfvox = mygrid.s_vt_vx[(matching_id[id][1]-1)%mygrid.nent, int(d_E2[str(idtr)])]
                    #par.append(voxpar * min(area(tr)/surfvox, 1))
                    
                    par.append(voxpar)
                    idtr += 1

                VTKtriangles(self.__my_scene, [par], ['PARa'], path+"triangles_PAR_"+str(ite)+".vtk")
        
        # ou alors sur une itération en particulier
        else:
            for i,tr in enumerate(self.__my_scene):
                df = self.__outputs[(self.__outputs.primitive_index == i)]
                par.append(df['PARa'].values[0])                

            VTKtriangles(self.__my_scene, [par], ['PARa'], path+"triangles_PAR_"+str(iteration)+".vtk")

            # VTK des voxels
            if self.__lightmodel == "ratp" and voxels:
                # plot dans VTK
                temp1, temp2, temp3 = [], [], []
                # éviter les éléments en trop
                for i in range(self.__ratp_scene.nveg)  :
                    for j in range(self.__ratp_scene.nje[i]):
                        dfvox = self.__outputs[(self.__outputs.VoxelId==i+1) & (self.__outputs.VegetationType==j+1) & (self.__outputs.Iteration==1)]
                        temp1.append(dfvox["PARa"].values[0])
                        temp2.append(self.__ratp_scene.nume[j,i])
                        temp3.append(int(i)+1) # kxyz sort en fortran
                para = [np.array(temp1), np.array(temp2), np.array(temp3)]

                RATP2VTK.RATPVOXELS2VTK(self.__ratp_scene, para, "PARa", path+"PARa_voxels.vtk")


    def s5(self):
        '''construit les fichiers d'entrée pour s5 et l'exécute
        '''
        # ecrit dans le dossier de s5
        f=open("s5/fort.51", 'w')
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

        f=open("s5/s5.par", 'w')
        
        f.write("%i\t%i\t%i\t%i\t0\n"%(self.__ratp_scene.nent, 9, 9, self.__ratp_scene.njz))
        for i in range(self.__ratp_scene.njz):
            f.write("%f\t"%(self.__ratp_scene.dz[i]))
        f.write("\n")

        njx = self.__ratp_scene.njx
        njy = self.__ratp_scene.njy
        dx = self.__ratp_scene.dx
        dy = self.__ratp_scene.dy
        f.write("%f\t%i\t%f\t%f\t%i\t%f\n"%(njx*dx, njx, dx, njy*dy, njy, dy))

        # on se place dans le dossier 
        #os.path.dirname(os.path.abspath(__file__))
        path = os.path.abspath(__file__)
        lpath = path.split("\\")
        path = "\\".join(lpath[:-2])
        path = path + "\\s5\\"
        print(path)
        os.chdir(path)
        subprocess.call(path+"s5.exe", cwd=path)
        print("--- Fin de s5.f")

    
    def s2v(self):
        '''construit les fichiers d'entrée pour s2v et l'exécute
        '''
        # ecrit dans le dossier de s2v
        f=open("s2v/fort.51", 'w')
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

        f=open("s2v/s2v.par", 'w')
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

        os.system(r"s2v\s2v++.exe")
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