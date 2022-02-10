import openalea.plantgl.all as pgl

from PyRATP.pyratp import grid
from PyRATP.pyratp import RATP2VTK
from PyRATP.pyratp.vegetation import Vegetation
from PyRATP.pyratp.skyvault import Skyvault
from PyRATP.pyratp.micrometeo import MicroMeteo
from PyRATP.pyratp.runratp import runRATP

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun, spitters_horaire

import os, subprocess
import itertools
import numpy as np
import pandas
import warnings

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
    f.write('FIELD FieldData '+str(len(varname))+'\n')
    for i, name in enumerate(varname):
        f.write(name+' 1 '+str(nbtr)+' float\n')
        for j in range(nbtr):
            f.write(str(var[i][j])+'\n')
        f.write('\n')

    f.write('\n')

    f.close()

class LightVegeManager:
    '''
    Classe qui rassemble les informations de la gestion de la lumière sur le couplage des modèles
    '''
    units = {'mm': 0.001, 'cm': 0.01, 'dm': 0.1, 'm': 1, 'dam': 10, 'hm': 100,'km': 1000}
    # transformation : [scale, translation]
    def __init__(self, in_scenes, 
                    in_transformations=[], 
                    in_names=[],
                    id_stems=[], 
                    lightmodel="ratp", 
                    lightmodelparam=[],
                    tesselation=True, 
                    pattern="none", 
                    patternwidth=0, 
                    scene_unit="m"):
        '''Constructeur
        
        Arg :
            in_scenes : liste de scenes plantgl, contient une scène plantGL par entité
            in_transformation : liste de [rescale factor(float), translation vector(Vector3), échelle de la scène parmi les units], 
                                transformations à appliquer pour chaque scène de in_scenes
            in_names : liste de string, nom des modèles dont vient chaque scene
            lightmodel : string, nom du modèle de calcul de la lumière
            lightmodelparam : liste mixte, paramètre en entrée du modèle de lumière
                pour RATP : [dx(float), dy(float), dz(float), latitude(float), longitude(float), timezone(float), rs(liste de float), tesselate_maxlevel(int), distribalgo(string), distriboption(dépend de distribalgo)]
                pour CARIBU : [diffuse_model(string), skyazimuths(float), skyzeniths(float), latitude(float), sun_sky_option(string)]
            id_stems : list(int), liste des shapes qui correspondent à une tige (pas de transmission de la lumière)
            tesselation : boolean, si on active la subdivision des triangles pour matcher une grille de voxels
            scene_unit : string, indique l'unité d'échelle du couplage
            
            pattern : string, active un algo de pattern sur les scènes en entrée (pas encore en place)
            patternwidth : distance entre les pattern (pas encore en place)

        '''
        self.__in_scenes = in_scenes
        self.__in_transformations = in_transformations
        self.__in_names = in_names
        self.__lightmodel = lightmodel
        self.__lightmodelparam = lightmodelparam
        self.__tesselation = tesselation
        self.__pattern = pattern
        self.__patternwidth = patternwidth
        self.__scene_unit = scene_unit
        self.__id_stems = id_stems

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
        
        # min-max de la scène
        xmax, xmin, ymax, ymin, zmax, zmin = 0,0,0,0,0,0            
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
        self.__pmax = Vector3(xmax, ymax, zmax)
        self.__pmin = Vector3(xmin, ymin, zmin)

        # création d'une scène plantGL
        # pas encore en place (manque d'abstraction pour le rescale)

        # création de la scène du modèle de lumière
        # RATP
        if lightmodel == "ratp":
            # récupère les paramètres
            dx, dy, dz, latitude, longitude, timezone, rs, levelmax, ele_algo, ele_option = lightmodelparam

            # nombre de voxels
            nx = int((self.__pmax[0] - self.__pmin[0]) // dx)
            ny = int((self.__pmax[1] - self.__pmin[1]) // dy)
            nz = int((self.__pmax[2] - self.__pmin[2]) // dz)
            if (self.__pmax[0] - self.__pmin[0]) % dx > 0 : nx += 1
            if (self.__pmax[1] - self.__pmin[1]) % dy > 0 : ny += 1
            if (self.__pmax[2] - self.__pmin[2]) % dz > 0 : nz += 1

            # définit une origine en Pmin 
            xorig, yorig, zorig = math.floor(self.__pmin[0]), math.floor(self.__pmin[1]), -math.floor(self.__pmin[2])

            # création de la grille
            mygrid = grid.Grid.initialise(nx, ny, nz, dx, dy, dz, xorig, yorig, zorig, latitude, longitude, timezone, len(in_scenes), rs)

            # subdivision des triangles pour matcher la grille
            if tesselation:
                # traite les triangles de la shape
                new_tr_scene=[]
                for tr in self.__my_scene:
                    level = 0
                    isworking = iterate_trianglesingrid(tr, mygrid, level, levelmax, new_tr_scene)
                
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
                a.append(tr.area)
                n.append(0.)
                entity.append(self.__matching_ids[tr.id][1])

            mygrid, matching = grid.Grid.fill_1(entity, barx, bary, barz, a, n, mygrid)
            self.__ratp_scene = mygrid
            self.__tr_vox = matching 

            ## distribution d'angles pour chaque entité
            distrib = []

            # calcul en dynamique
            # ele_option = nombre de classes
            if ele_algo == "compute" :
                # compte le nombre de triangles par entité
                t_nent_area=[]
                for k in range(len(self.__in_scenes)):
                    totA=0
                    for t in self.__my_scene:
                        if self.__matching_ids[t.id][1] == k:
                            totA+=t.area
                    t_nent_area.append(totA)
                
                # on ajoute 1 au nombre de classe avec range
                angles = list(range(10, 100, ele_option+1))

                # pour chaque entité
                for k in range(len(self.__in_scenes)):
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
                for k in range(len(self.__in_scenes)):
                    for i in range(len(distrib[k])):
                        distrib[k][i] *= 1/t_nent_area[k]
            
            # lecture du fichier
            # ele_option = chemin du fichier
            elif ele_algo == "file" : 
                f_angle = open(ele_option, 'r')
                for i in range(len(self.__in_scenes)):
                    line = f_angle.readline()
                    distrib.append([float(x) for x in line.split(',')[1:]])
            
            self.__ratp_distrib = distrib
    
        elif lightmodel == "caribu":
            self.__caribu_scene = {}
            for id,val in self.__matching_ids.items():
                self.__caribu_scene[id] = []

            for tr in self.__my_scene:
                l_tr=[]
                for i in range(3):
                    l_tr.append((tr[i][0], tr[i][1], tr[i][2]))
                self.__caribu_scene[tr.id].append(l_tr)

        else:
            raise ValueError("Unknown lightmodel : can be either 'ratp' or 'caribu' ")
            
    def __str__(self):
        '''imprime les infos de la scène couplée
        '''
        out = "\n-------------------------------------\n"
        out += "Light Vege Manager v0.0\n"
        out += str(len(self.__in_scenes))+" entités\n"
        out += "Modèle de lumière : "+self.__lightmodel
        if self.__tesselation : out += " avec tesselation\n"
        else : out += "\n"
        out += "Scènes en entrée : \n"
        for i in range(len(self.__in_scenes)):
            out+= "\t scene "+str(i)+" : "+str(len(self.__in_scenes[i]))+" shapes, modèle : "+str(self.__in_names[i])+"\n"

        out += "\n-------------------------------------\n"
        return out

    def run(self, meteo_path="", vege_path="", PARi=0, day=0, hour=0):
        '''calcul du bilan radiatif
        Args :
            meteo_path : string, chemin du fichier meteo
            vege_path : string, chemin du fichier grandeurs associées aux plantes (pour RATP)
            elevation_path : string, chemin du fichier contenant une distribution d'angle par entité (pour RATP)
            PARi : float, PARi en µmol.m-2.s-1
            day : float, day of the year 
            hour : float, hour TU

        return :
            enregistre les sorties dans self.__outputs sous forme de dataframe pandas
        '''
    
        # RATP
        if self.__lightmodel == "ratp":
            self.__inputspath = [meteo_path, vege_path]
            # fichiers d'entrée
            # initialisé en dur dans la classe
            sky = Skyvault.initialise()
            
            # initialisation des paramètres végétation
            f_veg = open(vege_path, 'r')
            # on se place sur la ligne avec les lambertiens
            f_veg.readline()
            f_veg.readline()
            f_veg.readline()
            f_veg.readline()
            f_veg.readline()
            rfline = f_veg.readline()

            # création d'un dict entity
            entities_param = []
            for dist_ent in self.__ratp_distrib:
                entities_param.append({
                                        'distinc' : dist_ent,
                                        'rf' : [float(x) for x in rfline.split('\t')[:2]]
                                        })

            vegetation = Vegetation.initialise(entities_param)

            # init météo
            met = MicroMeteo.read(meteo_path)

            # Calcul du bilan radiatif sur chaque pas de temps du fichier météo
            res = runRATP.DoIrradiation(self.__ratp_scene, vegetation, sky, met)

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
            
            # 'PAR' is expected in  Watt.m-2 in RATP input, whereas output is in micromol => convert back to W.m2 (cf shortwavebalance, line 306)
            # On enregistre tout dans une dataframe pandas
            dfvox =  pandas.DataFrame({'VegetationType':VegetationType,
                                'Iteration':Iteration,
                                'day':day,
                                'hour':hour,
                                'VoxelId':VoxelId,
                                'ShadedPAR':ShadedPAR / 4.6,
                                'SunlitPAR':SunlitPAR / 4.6,
                                'ShadedArea':ShadedArea,
                                'SunlitArea': SunlitArea,
                                'Area': ShadedArea + SunlitArea,
                                'PAR': (ShadedPAR * ShadedArea + SunlitPAR * SunlitArea) / (ShadedArea + SunlitArea) / 4.6,
                                'xintav': xintav, 
                                })
            
            # tri de la dataframe par rapport aux shapes et triangles
            
            # ne prend pas le sol
            dfvox = dfvox[dfvox['VegetationType'] > 0]
            
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
            nshapes = sum([len(l) for l in self.__in_scenes])
            s_shapes = []
            s_area=[]
            s_par=[]
            s_xintav=[]
            s_ite=[]
            s_day=[]
            s_hour=[]
            s_ent=[]
            for id in range(nshapes):
                # itérations commencent à 1
                for ite in range(int(max(output["Iteration"]))):
                    dffil = output[(output.Iteration == ite+1) & (output.shape_id == id)]
                    
                    s_hour.append(dffil["hour"].values[0])
                    s_day.append(dffil["day"].values[0])
                    s_ite.append(ite+1)
                    s_area.append(sum(dffil["Area"]))
                    s_par.append(sum(dffil['Area']*dffil['PAR']) / s_area[-1])
                    s_xintav.append(sum(dffil['Area']*dffil['xintav']) / s_area[-1])
                    s_ent.append(dffil["VegetationType"])
                    s_shapes.append(self.__matching_ids[id][0])
            self.__shape_outputs = pandas.DataFrame({
                "Iteration" : s_ite,
                "Day" : s_day,
                "Hour" : s_hour,
                "ShapeId" : s_shapes,
                "VegetationType" : s_ent,
                "Area" : s_area,
                "PAR" : s_par,
                "xintav" : s_xintav
            })

        # CARIBU
        elif self.__lightmodel == "caribu":

            ## Création du ciel et du soleil
            #: Diffuse light sources : Get the energy and positions of the source for each sector as a string
            sky_string = GetLight.GetLight(GenSky.GenSky()(1., self.__lightmodelparam[0], self.__lightmodelparam[1], self.__lightmodelparam[2]))  #: (Energy, soc/uoc, azimuts, zenits)

            # Convert string to list in order to be compatible with CaribuScene input format
            sky = []
            for string in sky_string.split('\n'):
                if len(string) != 0:
                    string_split = string.split(' ')
                    t = tuple((float(string_split[0]), tuple((float(string_split[1]), float(string_split[2]), float(string_split[3])))))
                    sky.append(t)

            #: Direct light sources (sun positions)
            sun = Gensun.Gensun()(1., day, hour, self.__lightmodelparam[3])
            sun = GetLightsSun.GetLightsSun(sun)
            sun_str_split = sun.split(' ')
            sun = [tuple((float(sun_str_split[0]), tuple((float(sun_str_split[1]), float(sun_str_split[2]), float(sun_str_split[3])))))]
            
            #: Optical properties
            opt = {'par': {}}

            for id, val in self.__matching_ids.items():
                # c'est une tige il n'y a pas de transmission
                if (val[0], val[1]) in self.__id_stems:
                    opt['par'][id] = (0.10,)  #: (reflectance,) of the stems
                else:
                    opt['par'][id] = (0.10, 0.05) #: (reflectance, transmittance) of the adaxial side of the leaves

            # construction d'une scène ciel et soleil
            c_scene_sky = CaribuScene(scene=self.__caribu_scene, light=sky, opt=opt, scene_unit=self.__scene_unit)
            c_scene_sun = CaribuScene(scene=self.__caribu_scene, light=sun, opt=opt, scene_unit=self.__scene_unit)

            # choix du type de rayonnement
            # direct=False : active la rediffusion
            # infinite=False : désactive la répétition infinie de la scène (pas de pattern défini ici)
            if self.__lightmodelparam[4] == "mix":
                raw_sun, aggregated_sun = c_scene_sun.run(direct=True)
                Erel_sun = aggregated_sun['par']['Eabs']
                raw_sky, aggregated_sky = c_scene_sky.run(direct=True)
                Erel_sky = aggregated_sky['par']['Eabs']

                #: Spitters's model estimating for the diffuse:direct ratio
                Rg = PARi / 2.02  #: Global Radiation (W.m-2)
                RdRs = spitters_horaire.RdRsH(Rg=Rg, DOY=day, heureTU=hour, latitude=self.__lightmodelparam[3])  #: Diffuse fraction of the global irradiance
                Erel = {}
                for element_id, Erel_value in Erel_sky.items():
                    Erel[element_id] = RdRs * Erel_value + (1 - RdRs) * Erel_sun[element_id]

                Erel_output_shape = Erel
                PARa_output_shape = {k: v * PARi for k, v in Erel_output_shape.items()}
                
                # même numération que caribu
                Erel_output_tr = {}
                count=0
                for key, val in raw_sky['par']['Eabs'].items():
                    for i,par in enumerate(val):
                        Erel_output_tr[count] = RdRs*par+(1-RdRs)*raw_sun['par']['Eabs'][key][i]
                        count+=1
                PARa_output_tr = {k: v * PARi for k, v in Erel_output_tr.items()}
            
            elif self.__lightmodelparam[4] == "sun":
                raw_sun, aggregated_sun = c_scene_sun.run()
                Erel_sun = aggregated_sun['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                PARa_sun = {k: v * PARi for k, v in Erel_sun.items()}
                Erel_output = Erel_sun
                PARa_output = PARa_sun
            
            elif self.__lightmodelparam[4] == "sky":
                raw_sky, aggregated_sky = c_scene_sky.run()
                Erel_sky = aggregated_sky['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                PARa_sky = {k: v * PARi for k, v in Erel_sky.items()}
                Erel_output = Erel_sky
                PARa_output = PARa_sky

            else:
                raise ValueError("Unknown sun_sky_option : can be either 'mix', 'sun' or 'sky'.")
        
            # enregistre les valeurs par shape et plantes
            s_shapes = []
            s_area=[]
            s_par=[]
            s_xintav=[]
            s_day=[]
            s_hour=[]
            s_ent=[]
            for key,val in Erel_output_shape.items():
                s_shapes.append(self.__matching_ids[key][0])
                s_par.append(PARa_output_shape[key])
                s_xintav.append(val)
                s_area.append(aggregated_sun['par']['area'][key])
                s_day.append(day)
                s_hour.append(hour)
                s_ent.append(self.__matching_ids[key][1])

            self.__shape_outputs = pandas.DataFrame({
                "Day" : s_day,
                "Hour" : s_hour,
                "ShapeId" : s_shapes,
                "VegetationType" : s_ent,
                "Area" : s_area,
                "PAR" : s_par,
                "xintav" : s_xintav
            })

            # enregistre les valeurs par triangles
            s_shapes = []
            s_tr = []
            s_area=[]
            s_par=[]
            s_xintav=[]
            s_day=[]
            s_hour=[]
            s_ent=[]
            for key,val in Erel_output_tr.items():
                s_shapes.append(self.__matching_ids[self.__my_scene[key].id][0])
                s_par.append(PARa_output_tr[key])
                s_xintav.append(val)
                s_tr.append(key)
                s_day.append(day)
                s_hour.append(hour)
                s_ent.append(self.__matching_ids[self.__my_scene[key].id][1])
                # search id in shape
                id_in_sh= self.__matching_ids[self.__my_scene[key].id][2].index(key) 
                s_area.append(raw_sun['par']['area'][self.__my_scene[key].id][id_in_sh])

            self.__outputs = pandas.DataFrame({
                "Day" : s_day,
                "Hour" : s_hour,
                "primitive_index" : s_tr,
                "ShapeId" : s_shapes,
                "VegetationType" : s_ent,
                "Area" : s_area,
                "PAR" : s_par,
                "xintav" : s_xintav
            })

        else:
            raise ValueError("Unknown lightmodel : can be either 'ratp' or 'caribu' ")

    @property
    def shapes_outputs(self):
        return self.__shape_outputs 

    @property
    def triangles_outputs(self):
        return self.__outputs 

    def PAR_update_MTG(self, mtg):
        # crée un tableau comme dans caribu_facade de fspm-wheat
        dico_par = {}
        para_dic = {}
        erel_dic = {}
        for s in self.__shape_outputs["ShapeId"]:
            para_dic[s] = self.__shape_outputs[self.__shape_outputs.ShapeId==s]["PAR"].values[0]
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

    def VTKout(self, path, iteration=None):
        '''construit des fichiers VTK de la triangulation avec les valeurs de PAR associées

        Args :
            self.__my_scene : liste de Triangle3
            self.__outputs : dataframe des résultats

            path : string, chemin pour l'écriture des fichiers
        '''    
        par = []

        # si la lumière a été calculée sur plusieurs itération
        if iteration is None:
            for ite in range(int(max(self.__outputs['Iteration']))):
                par = []
                idtr=0
                for tr in self.__my_scene:
                    df = self.__outputs[(self.__outputs.Iteration == ite+1) & (self.__outputs.primitive_index == idtr)]
                    voxpar = df['PAR'].values[0]
                    
                    # tentative de réduction du par selon l'aire du triangle mais peut etre déjà fait dans le tri de df
                    #surfvox = mygrid.s_vt_vx[matching_id[id][1], int(d_E2[str(idtr)])]
                    #if surfvox==0 : surfvox = mygrid.s_vt_vx[(matching_id[id][1]-1)%mygrid.nent, int(d_E2[str(idtr)])]
                    #par.append(voxpar * min(area(tr)/surfvox, 1))
                    
                    par.append(voxpar)
                    idtr += 1

                VTKtriangles(self.__my_scene, [par], ['PAR'], path+"triangles_PAR_"+str(ite)+".vtk")
        
        # ou alors sur une itération en particulier
        else:
            for i,tr in enumerate(self.__my_scene):
                df = self.__outputs[(self.__outputs.primitive_index == i)]
                par.append(df['PAR'].values[0])                

            VTKtriangles(self.__my_scene, [par], ['PAR'], path+"triangles_PAR_"+str(iteration)+".vtk")


    def s5(self):
        '''construit les fichiers d'entrée pour s5 et l'exécute
        '''
        # ecrit dans le dossier de s5
        f=open("s5/fort.51", 'w')
        c_tr=1
        for tr in self.__my_scene:
            label = str(self.__matching_ids[tr.id][1]+1)+str("%05i"%(self.__matching_ids[tr.id][1]+1))+'001001'#+str("%03i"%(c_tr))
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
        os.chdir(path)
        os
        subprocess.call(path+"s5.exe", cwd=path)
        print("--- Fin de s5.f")

    
    def s2v(self):
        '''construit les fichiers d'entrée pour s2v et l'exécute
        '''
        # ecrit dans le dossier de s2v
        f=open("s2v/fort.51", 'w')
        c_tr=1
        for tr in self.__my_scene:
            label = str(self.__matching_ids[tr.id][1]+1)+str("%05i"%(self.__matching_ids[tr.id][1]+1))+'001001'#+str("%03i"%(c_tr))
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
