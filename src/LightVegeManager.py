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

    def __init__(self, in_scenes, 
                    in_transformations=[], 
                    in_names=[],
                    id_stems=[], 
                    lightmodel="ratp", 
                    lightmodelparam=[],
                    sky_parameters=[],
                    rf = [],
                    coordinates=[46.58, 0.38, 0.], # Poitiers
                    global_scene_tesselate_level=0,
                    multithread=0,
                    scene_unit="m"):
        '''Constructeur
        
        Arg :
            in_scenes : liste de scenes plantgl, contient une scène plantGL par entité
            in_transformation : liste pour chaque entité de 4 éléments :
                                * rescale factor(float)
                                * translation vector(Vector3) 
                                * échelle de la scène parmi le dico units
                                * convetion point cardinal de x+ (x+ = N par défaut) 
                            transformations à appliquer pour chaque scène de in_scenes
            in_names : liste de string, nom des modèles (ou autre informations) dont vient chaque scene
            lightmodel : string, nom du modèle de calcul de la lumière 'ratp' ou 'caribu'
            lightmodelparam : liste mixte, paramètre en entrée du modèle de lumière
                pour RATP : [dx(float), dy(float), dz(float), : dimension d'un voxel (dans l'unité de scene_unit)
                            rs(liste de 2 float), : réflectance du sol dans le PAR et le NIR
                            mu(list de float) : clumping pour chaque entité
                            tesselate_maxlevel(int), : niveau de subdivision des triangles pour rentrer dans la grille
                            distribalgo(string), "compute" : calcul la distribution d'élévation des triangles
                                                "file" : lit un fichier pour avoir la distribution d'élévation de chaque scène
                            distriboption(dépend de distribalgo) : si "compute" = nombre de classe
                                                                    si "file" = chemin du fichier
                            infinite : boolean, pour activer ou non le couvert infini 
                            ] 
                pour CARIBU : [sun_sky_option(string) : 'sun' ou 'sky' ou 'mix' pour le type de rayonnement,
                                sun_algo(string) : "caribu" ou "ratp : indique la manière de calculer la direction du soleil
                                infinite : boolean, pour activer ou non le couvert infini
                                domain_pattern : si l'option couvert infini est activé, défini le domaine (xy) du pattern à reproduire
                                ]
            id_stems : liste de tuple(nshape(int), nentité(int)), liste des shapes de l'entité d'entrée qui correspondent à une tige
            sky_parameters : liste, plusieurs possibilité selon le modele de lumière :
                                    [] : ciel par défaut est une turtle à 46 directions pour CARIBU et RATP
                                    [0] : annule de diffus pour RATP
                                    [chemin fichier ciel] : valable que pour RATP
                                    [n_azimuth, n_zenith] 
                                    [n_azimuth, n_zenith, type of sky] 
                                    [n_azimuth, n_zenith, weights]
            rf : list, [[reflectancePAR_ent1, reflectanceNIR_ent1], ..., [reflectancePAR_ent1, reflectanceNIR_ent1]] 
                 reflectances PAR;NIR ou reflectance;transmittance selon le modèle
            coordinates : list(float), latitude, longitude, timezone
            scene_unit : string, indique l'unité d'échelle du couplage
            global_scene_tesselate_level : tesselation de la scène couplée globale, quelque soit le modèle de lumière

            Notes :
                Convention de coordonnées (type RATP scene) : x+ = N, y+ = W 
                valable dans self.__myscene et self.__caribu_scene

        '''
        self.__in_scenes = in_scenes
        self.__in_transformations = in_transformations
        self.__in_names = in_names
        self.__lightmodel = lightmodel
        self.__lightmodelparam = lightmodelparam
        self.__scene_unit = scene_unit
        self.__id_stems = id_stems
        self.__rf = rf
        self.__coordinates = coordinates
        self.__diffus=True

        # couplage des scènes dans un plantGL commun et création du tableau des ids
        # [plantGL] -(+)-> [Triangle3] -> transformation sur les Triangle3 -> [Triangle3]
        # une scène plantgl par espèce
        self.__matching_ids = {}
        self.__my_scene=[]
        count=0
        for i_esp, scene in enumerate(in_scenes) :
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
        for i_esp, trans in enumerate(in_transformations):
            scale = (trans[0] != "none") and (trans[0] != "None")
            transl = (trans[1] != "none") and (trans[1] != "None")
            for tr in self.__my_scene:
                if self.__matching_ids[tr.id][1] == i_esp:
                    if scale : tr.rescale(trans[0])
                    if transl : tr.translate(trans[1])
                    if (trans[2] != scene_unit) and (trans[2] in self.units):
                        tr.rescale(self.units[trans[2]]/self.units[scene_unit])
                    # la convention du repère xyz par rapport aux points cardinaux est précisée
                    # on ramène la scène à la convention x+ = N
                    # si non précisé, on ne change pas l'orientation de la scène
                    if len(trans) > 3:
                        if trans[3] == "x+ = S":
                            tr.zrotate(180)
                        if trans[3] == "x+ = W":
                            tr.zrotate(90)
                        if trans[3] == "x+ = E":
                            tr.zrotate(-90)
        
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

        # création d'une scène plantGL
        # pas encore en place (manque d'abstraction pour le rescale)
        
        # création de la scène du modèle de lumière
        # RATP
        if lightmodel == "ratp":
            # récupère les paramètres
            dx, dy, dz, rs, mu, levelmax, ele_algo, ele_option, infinite = lightmodelparam     

            self.__ratp_mu = mu

            # on sépare les tiges dans une nouvelle entité
            mem_mu=[]
            for stem in id_stems:
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
                                self.__matching_ids[ids][1] + len(self.__in_scenes),
                                self.__matching_ids[ids][2])
                self.__matching_ids[ids] = tuple_temp

                if stem[1] not in mem_mu:
                    self.__ratp_mu.append(self.__ratp_mu[stem[1]])
                    self.__rf.append(self.__rf[stem[1]])
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
            if ele_algo == "compute global" :
                # compte le nombre de triangles par entité
                t_nent_area=[]
                for k in range(nent):
                    totA=0
                    for t in self.__my_scene:
                        if self.__matching_ids[t.id][1] == k:
                            totA+=t.area
                    t_nent_area.append(totA)
                
                # on ajoute 1 au nombre de classe avec range
                angles = list(np.linspace(90/ele_option, 90, ele_option))

                # pour chaque entité
                for k in range(nent):
                    classes = [0] * ele_option
                    # parcourt les triangles
                    for t in self.__my_scene:
                        # pour chaque triangle de l'entité
                        if self.__matching_ids[t.id][1] == k:
                            # recherche de la classe
                            i=0
                            while i<ele_option:
                                if t.elevation < angles[i]:
                                    classes[i] += t.area
                                    # pour sortir de la boucle
                                    i=ele_option+10
                                i+=1

                    distrib.append(classes)

                # convertit en pourcentage
                for k in range(nent):
                    for i in range(len(distrib[k])):
                        distrib[k][i] *= 1/t_nent_area[k]
            
            # lecture du fichier
            # ele_option = chemin du fichier
            elif ele_algo == "file" : 
                f_angle = open(ele_option, 'r')
                for i in range(len(self.__in_scenes)):
                    line = f_angle.readline()
                    distrib.append([float(x) for x in line.split(',')[1:]])

            # on ajuste au besoin les min-max si la scène est plane pour avoir un espace 3D
            if xmin==xmax:
                xmax+=dx
                xmin-=dx 
            if ymin==ymax:
                ymax+=dy
                ymin-=dy 
            if zmin==zmax:
                zmax+=dz
                zmin-=dz   
            
            self.__pmax = Vector3(xmax, ymax, zmax)
            self.__pmin = Vector3(xmin, ymin, zmin)

            # nombre de voxels
            nx = int((self.__pmax[0] - self.__pmin[0]) // dx)
            ny = int((self.__pmax[1] - self.__pmin[1]) // dy)
            nz = int((self.__pmax[2] - self.__pmin[2]) // dz)
            if (self.__pmax[0] - self.__pmin[0]) % dx > 0 : nx += 1
            if (self.__pmax[1] - self.__pmin[1]) % dy > 0 : ny += 1
            if (self.__pmax[2] - self.__pmin[2]) % dz > 0 : nz += 1
            
            # définit une origine en Pmin 
            xorig, yorig, zorig = self.__pmin[0], self.__pmin[1], -self.__pmin[2]

            # création de la grille
            mygrid = grid.Grid.initialise(nx, ny, nz, dx, dy, dz, xorig, yorig, zorig, coordinates[0], coordinates[1], coordinates[2], nent, rs, toric=infinite)

            # subdivision des triangles pour matcher la grille
            if levelmax>0:
                # traite les triangles de la shape
                new_tr_scene=[]
                start=time.time()
                # traitement multithreading
                if multithread > 0:
                    level = 0
                    tesseletorgrid = MyTessetelorGrid(mygrid, level, levelmax, new_tr_scene)
                    with concurrent.futures.ThreadPoolExecutor(multithread) as executor:
                        executor.map(tesseletorgrid._iterate_trianglesingrid, self.__my_scene)
                else :
                    for tr in self.__my_scene:
                        level = 0
                        isworking = iterate_trianglesingrid(tr, mygrid, level, levelmax, new_tr_scene)
                print("tesselation time : ",time.time()-start)
                # copie de la nouvelle triangulation
                self.__my_scene = new_tr_scene
            
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
                if (self.__matching_ids[tr.id][0], self.__matching_ids[tr.id][1] - len(self.__in_scenes)) in self.__id_stems:
                    a.append(tr.area/2)
                else:
                    a.append(tr.area)
                
                n.append(0.)
                entity.append(self.__matching_ids[tr.id][1])

            mygrid, matching = grid.Grid.fill_1(entity, barx, bary, barz, a, n, mygrid)
            self.__ratp_scene = mygrid
            self.__tr_vox = matching 

            if ele_algo == "compute voxel":
                angles = list(np.linspace(90/ele_option, 90, ele_option))
                t_area=[]
                # pour chaque voxel
                for k in range(self.__ratp_scene.nveg):
                    t_nent_area=[]
                    distrib_nent=[]
                    # pour chaque entité
                    for n in range(nent):
                        areatot=0
                        classes = [0] * ele_option
                        # on parcourt le dico des correspondances triangle->voxel
                        for idt, idv in self.__tr_vox.items():
                            if idv == k and self.__my_scene[int(idt)].id == n:
                                areatot += self.__my_scene[int(idt)].area
                                
                                # recherche de la classe
                                i=0
                                while i<ele_option:
                                    if self.__my_scene[int(idt)].elevation < angles[i]:
                                        classes[i] += self.__my_scene[int(idt)].area
                                        # pour sortir de la boucle
                                        i=ele_option+10
                                    i+=1
                        t_nent_area.append(areatot)
                        distrib_nent.append(classes)
                    t_area.append(t_nent_area)
                    distrib.append(distrib_nent)

                # convertit en pourcentage
                for k in range(self.__ratp_scene.nveg):
                    for n in range(nent):
                        for i in range(len(distrib[k])):
                            distrib[k][n][i] *= 1/t_area[k][n]

            self.__ratp_distrib = distrib

            # construction du ciel:
            if sky_parameters==[]:
                # par défaut ciel à 46 directions ("tortue")
                self.__sky = Skyvault.initialise()

            # annule le diffus
            elif sky_parameters[0]==0:
                self.__sky = Skyvault.initialise()
                self.__diffus = False
            
            # ciel dans un fichier à lire
            elif len(sky_parameters)==1:
                self.__sky = Skyvault.read(sky_parameters[0])
            
            # on utilise un outil de caribu pour subdiviser le ciel
            elif len(sky_parameters)==3:
                ele=[]
                azi=[]
                da = 2 * math.pi / sky_parameters[0]
                dz = math.pi / 2 / sky_parameters[1]      
                todeg = 180/math.pi          
                for j in range(sky_parameters[0]):
                    for k in range(sky_parameters[1]):
                        azi.append((j * da + da / 2)*todeg)
                        ele.append((k * dz + dz / 2)*todeg)
                n = sky_parameters[0]*sky_parameters[1]
                
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
                for j in range(sky_parameters[0]):
                    for k in range (sky_parameters[1]):
                        azim,elv = j * da + da / 2, k * dz + dz / 2
                        I=0
                        if(sky_parameters[2]=='soc'):
                            I = soc(elv, dz, da)
                        else:
                            I = uoc(elv, dz, da)

                        pc.append(I) 
                
                
                self.__sky = Skyvault.initialise(ele, azi, omega, pc)
            
            
    
        elif lightmodel == "caribu":
            self.__pmax = Vector3(xmax, ymax, zmax)
            self.__pmin = Vector3(xmin, ymin, zmin)

            # création d'une triangulation caribu
            self.__caribu_scene = {}
            for id,val in self.__matching_ids.items():
                self.__caribu_scene[id] = []

            # la convention de CARIBU avec x+ = N, pas de changement de base
            for tr in self.__my_scene:
                l_tr=[]
                for i in range(3):
                    l_tr.append((tr[i][0], tr[i][1], tr[i][2]))
                self.__caribu_scene[tr.id].append(l_tr)
            
            # création du ciel
            # par défaut turtle 46 directions, soc
            if sky_parameters==[]:
                turtle_list = turtle.turtle()
                self.__sky=[]
                for i,e in enumerate(turtle_list[0]):
                    t = tuple((e, tuple((turtle_list[2][i][0], turtle_list[2][i][1], turtle_list[2][i][2]))))
                    self.__sky.append(t)

            else:
                sky_string = GetLight.GetLight(GenSky.GenSky()(1., sky_parameters[2], sky_parameters[0], sky_parameters[1]))  #: (Energy, soc/uoc, azimuts, zenits)
            
                # Convert string to list in order to be compatible with CaribuScene input format
                self.__sky = []
                for string in sky_string.split('\n'):
                    if len(string) != 0:
                        string_split = string.split(' ')
                        t = tuple((float(string_split[0]), tuple((float(string_split[1]), float(string_split[2]), float(string_split[3])))))
                        self.__sky.append(t)

        else:
            raise ValueError("Unknown lightmodel : can be either 'ratp' or 'caribu' ")
            
    def run(self, meteo_path="", PARi=0, day=0, hour=0, parunit="micromol.m-2.s-1", domain= [], truesolartime=False, printsun=False):
        '''calcul du bilan radiatif
        Args :
            meteo_path : string, chemin du fichier meteo
            vege_path : string, chemin du fichier grandeurs associées aux plantes (pour RATP)
            elevation_path : string, chemin du fichier contenant une distribution d'angle par entité (pour RATP)
            PARi : float, PARi en µmol.m-2.s-1 ou W.m-2
            day : float, day of the year 
            hour : float, hour TU
            parunit : "micromol.m-2.s-1" ou "W.m-2"
            domain : si CARIBU est en infinitisation, définit le pattern de la CaribuScene
            truesolartime : boolean, si hour est l'heure solaire réel ou l'heure locale

        return :
            enregistre les sorties dans self.__outputs sous forme de dataframes pandas
            sorties en µmol.m-2.s-1
        '''
    
        # RATP
        if self.__lightmodel == "ratp":
            # création d'un dict entity
            entities_param = []
            if self.__lightmodelparam[6] != "compute voxel":
                for id, dist_ent in enumerate(self.__ratp_distrib):
                    entities_param.append({
                                            'mu' : self.__ratp_mu[id],
                                            'distinc' : dist_ent,
                                            'rf' : self.__rf[id]
                                            })

            vegetation = Vegetation.initialise(entities_param)

            # init météo, PARi en entrée en W.m-2
            if meteo_path == "":
                if parunit == "micromol.m-2.s-1" :
                    #: Spitters's model estimating for the diffuse:direct ratio
                    # coefficient 2.02 : 4.6 (conversion en W.m-2) x 0.439 (PAR -> global)
                    RdRs = spitters_horaire.RdRsH(Rg=PARi/2.02, DOY=day, heureTU=hour, latitude=self.__coordinates[0])
                    
                    # coeff 4.6 : https://www.researchgate.net/post/Can-I-convert-PAR-photo-active-radiation-value-of-micro-mole-M2-S-to-Solar-radiation-in-Watt-m2
                    PARi = PARi / 4.6  # W.m-2
                else:
                    RdRs = spitters_horaire.RdRsH(Rg=PARi/0.439, DOY=day, heureTU=hour, latitude=self.__coordinates[0])
                
                # PAR et Dif en W.m^-2
                if self.__diffus : 
                    met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=PARi, Rdif=PARi*RdRs, truesolartime=truesolartime)
                else :  met = MicroMeteo.initialise(doy=day, hour=hour, Rglob=PARi, Rdif=0, truesolartime=truesolartime)
            
            # lecture d'un fichier
            else:
                met = MicroMeteo.read(meteo_path, truesolartime)

            # Calcul du bilan radiatif sur chaque pas de temps du fichier météo
            res = runRATP.DoIrradiation(self.__ratp_scene, vegetation, self.__sky, met)

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
                                'PARa': (ShadedPAR * ShadedArea + SunlitPAR * SunlitArea) / (ShadedArea + SunlitArea),
                                'xintav': xintav, 
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
            s_par=[]
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
                for ite in range(int(max(output["Iteration"]))):
                    dffil = output[(output.Iteration == ite+1) & (output.shape_id == id)]
                    
                    s_hour.append(dffil["hour"].values[0])
                    s_day.append(dffil["day"].values[0])
                    s_ite.append(ite+1)
                    s_area.append(sum(dffil["primitive_area"]))
                    s_par.append(sum(dffil['primitive_area']*dffil['PARa']) / s_area[-1])
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
                "PARa" : s_par,
                "xintav" : s_xintav,
                "SunlitPAR" : s_parsun,
                "SunlitArea" : s_areasun,
                "ShadedPAR" : s_parsha,
                "ShadedArea" : s_areasha
            })

        # CARIBU
        elif self.__lightmodel == "caribu":
            if len(self.__lightmodelparam)==3 : sun_sky_option, sun_algo, infinite = self.__lightmodelparam
            elif len(self.__lightmodelparam)==4 : sun_sky_option, sun_algo, infinite, domain = self.__lightmodelparam

            ## Création du ciel et du soleil
            #: Direct light sources (sun positions)
            if sun_algo=="ratp":
                from PyRATP.pyratp import pyratp
                
                az, ele = 5.,9. # variables fantômes (non récupérées)
                pyratp.shortwave_balance.sundirection(ele, az, self.__coordinates[0], self.__coordinates[1], self.__coordinates[2], day, hour, truesolartime)            
                degtorad = math.pi/180

                azrad = pyratp.shortwave_balance.azdeg
                # peut avoir un nan à 12h 
                if math.isnan(azrad) and hour==12. : 
                    # critère empirique basé sur l'algo de CARIBU (GenSun)
                    if self.__coordinates[0] >= 21.11:
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
                    hour =hour+self.__coordinates[2]+self.__coordinates[1]/15.-eqntime/60. % 24.

                # Conversion de la latitude en radian!
                sun = Gensun.Gensun()(1., day, hour, self.__coordinates[0]*math.pi/180)
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
                    if (val[0], val[1]) in self.__id_stems:
                        opt['par'][id] = (self.__rf[val[1]][0],)  #: (reflectance,) of the stems ou organe opaque
                    else:
                        opt['par'][id] = (self.__rf[val[1]][0], self.__rf[val[1]][1]) #: (reflectance, transmittance) of the adaxial side of the leaves, élément translucide symétrique

                # construction d'une scène ciel et soleil
                if infinite : # on ajoute un domaine pour la création du pattern
                    c_scene_sky = CaribuScene(scene=self.__caribu_scene, light=self.__sky, opt=opt, scene_unit=self.__scene_unit, pattern=domain)
                    c_scene_sun = CaribuScene(scene=self.__caribu_scene, light=self.__sun, opt=opt, scene_unit=self.__scene_unit, pattern=domain)
                else:
                    c_scene_sky = CaribuScene(scene=self.__caribu_scene, light=self.__sky, opt=opt, scene_unit=self.__scene_unit)
                    c_scene_sun = CaribuScene(scene=self.__caribu_scene, light=self.__sun, opt=opt, scene_unit=self.__scene_unit)

                # choix du type de rayonnement
                # direct=False : active la rediffusion
                # infinite=False : désactive la répétition infinie de la scène (pas de pattern défini ici)
                if sun_sky_option == "mix":
                    raw_sun, aggregated_sun = c_scene_sun.run(direct=True, infinite=infinite)
                    Erel_sun = aggregated_sun['par']['Eabs']
                    Ei_sun = aggregated_sun['par']['Ei']
                    raw_sky, aggregated_sky = c_scene_sky.run(direct=True, infinite=infinite)
                    Erel_sky = aggregated_sky['par']['Eabs']
                    Ei_sky = aggregated_sky['par']['Ei']

                    #: Spitters's model estimating for the diffuse:direct ratio
                    # % de sky dans la valeur d'énergie finale
                    Rg = PARi / 2.02  #: Global Radiation (W.m-2)
                    RdRs = spitters_horaire.RdRsH(Rg=Rg, DOY=day, heureTU=hour, latitude=self.__coordinates[0])  #: Diffuse fraction of the global irradiance
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
                    raw_sun, aggregated_sun = c_scene_sun.run(direct=True, infinite=infinite)
                    Erel_output_shape = aggregated_sun['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                    Ei_output_shape = aggregated_sun['par']['Ei']
                    PARa_output_shape = {k: v * PARi for k, v in Erel_output_shape.items()}
                    PARi_output_shape = {k: v * PARi for k, v in Ei_output_shape.items()}
                    
                    # dataframe des triangles
                    Erel_output_tr = {}
                    Eri_output_tr = {}
                    count=0
                    for key, val in raw_sun['par']['Eabs'].items():
                        for i,par in enumerate(val):
                            Erel_output_tr[count] = raw_sun['par']['Eabs'][key][i]
                            Eri_output_tr[count] = raw_sun['par']['Ei'][key][i]
                            count+=1
                    PARa_output_tr = {k: v * PARi for k, v in Erel_output_tr.items()}
                    PARi_output_tr = {k: v * PARi for k, v in Eri_output_tr.items()}
                
                elif sun_sky_option == "sky":
                    raw_sky, aggregated_sky = c_scene_sky.run(direct=True, infinite=infinite)
                    Erel_output_shape = aggregated_sky['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                    Ei_output_shape = aggregated_sky['par']['Ei']
                    PARa_output_shape = {k: v * PARi for k, v in Erel_output_shape.items()}
                    PARi_output_shape = {k: v * PARi for k, v in Ei_output_shape.items()}
                    
                    # dataframe des triangles
                    Erel_output_tr = {}
                    Eri_output_tr = {}
                    count=0
                    for key, val in raw_sky['par']['Eabs'].items():
                        for i,par in enumerate(val):
                            Erel_output_tr[count] = raw_sky['par']['Eabs'][key][i]
                            Eri_output_tr[count] = raw_sky['par']['Ei'][key][i]
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
                if sun_up: s_area[key] = aggregated_sun['par']['area'][key]
                s_day[key] = day
                s_hour[key] = hour
                s_ent[key] = self.__matching_ids[key][1]

            self.__shape_outputs = pandas.DataFrame({
                "Day" : s_day,
                "Hour" : s_hour,
                "ShapeId" : s_shapes,
                "VegetationType" : s_ent,
                "Area" : s_area,
                "PARa" : s_par,
                "PARi" : s_pari,
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
                "PARa" : s_par,
                "PARi" : s_pari,
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
        print("--- day: %i \t hour: %i \t latitude: %0.2f °"%(day, hour, self.__coordinates[0]))
        
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
            caribuhour =hour+self.__coordinates[2]+self.__coordinates[1]/15.-eqntime/60. % 24. 
            print("--- local time, true solar time is %0.2f"%(caribuhour))
        
        # CARIBU et RATP sont dans la même convention de coordonnées, x+ = North
        print("--- RATP ---")
        az, ele = 5.,9. # variables fantômes (non récupérées)
        pyratp.shortwave_balance.sundirection(ele, az, self.__coordinates[0], self.__coordinates[1], self.__coordinates[2], day, hour, truesolartime)        
        degtorad = math.pi/180
        azrad = pyratp.shortwave_balance.azdeg
        if math.isnan(azrad) and hour==12. : 
            if self.__coordinates[0] >= 21.11:
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
        suncaribu._set_pos_astro(day, caribuhour, self.__coordinates[0]*math.pi/180)
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
        
        if self.__lightmodel == "ratp":
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
            for i in range(self.__ratp_scene.nent):
                plantnames.append("plant_"+str(i))
        
        # pour chaque plante on a une valeur par triangle
        if planttrianglevalues==[]:
            for i in range(self.__ratp_scene.nent):
                planttrianglevalues.append([])

            for tr in self.__my_scene :
                planttrianglevalues[self.__matching_ids[tr.id][1]].append(10)
                for i in range(self.__ratp_scene.nent):
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