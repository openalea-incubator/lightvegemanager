'''
classe globale de l'outil, directions principales des itérations
les autres ne sont que des outils supplémentaires

'''

from PyRATP.pyratp.runratp import runRATP
from alinea.caribu.CaribuScene import CaribuScene
import legume.RIRI5 as riri

import time
import os
import subprocess

from src.LightVegeManager_trianglesmesh import *
from src.LightVegeManager_stems import *
from src.LightVegeManager_leafangles import *
from src.LightVegeManager_sun import *
from src.LightVegeManager_sky import *
from src.LightVegeManager_RATPinputs import *
from src.LightVegeManager_CARIBUinputs import *
from src.LightVegeManager_outputs import *
from src.LightVegeManager_buildRATPscene import *
from src.LightVegeManager_transfer import *
from src.LightVegeManager_VTK import *
from src.LightVegeManager_defaultvalues import *

class LightVegeManager :
    
    ## MAIN METHODS ##
    def __init__(self,
                    environment={},
                    lightmodel="",
                    lightmodel_parameters={},
                    main_unit="m") :
        # on récupère les informations des autres modules de la classe
        default_environnement,    \
        default_ratp_parameters,  \
        default_caribu_parameters = default_LightVegeManager_inputs()

        # vérifie que le modèle de lumière choisie est valide
        if lightmodel != "ratp" and lightmodel != "caribu" :
            raise ValueError("Unknown lightmodel: can be either 'ratp' \
                                or 'caribu' ")

        # enregistre les informations fixes en entrée
        self.__environment = environment
        self.__lightmodel = lightmodel
        self.__lightmodel_parameters = lightmodel_parameters
        self.__main_unit = main_unit

        # copie des valeurs par défaut
        for key, value in default_environnement.items() :
            if key not in self.__environment : 
                self.__environment[key] = value
        
        if lightmodel == "caribu" : lghtdict = default_caribu_parameters
        elif lightmodel == "ratp" : lghtdict = default_ratp_parameters
        for key, value in lghtdict.items() :
            if key not in self.__lightmodel_parameters : 
                self.__lightmodel_parameters[key] = value

        # construction du ciel si besoin
        skytype = self.__environment["sky"]
        if lightmodel == "caribu":
            if self.__environment["diffus"] : self.__sky = CARIBUsky(skytype)
        elif lightmodel == "ratp" : self.__sky = RATPsky(skytype)     

    def build(self, geometry = {}, global_scene_tesselate_level=0) :
        self.__geometry = geometry

        # première lecture de la liste de scène
        self.__complete_trimesh, \
        self.__matching_ids,      \
        legume_grid,               \
        id_legume_scene = chain_triangulations(self.__geometry["scenes"])

        # transformations géométriques des scènes
        if "transformations" in self.__geometry :
            apply_transformations(self.__complete_trimesh, 
                                    self.__matching_ids, 
                                    self.__geometry["transformations"], 
                                    self.__main_unit)

        self.__areamax = compute_area_max(self.__complete_trimesh)

        # tesselation globale de la triangulation
        if self.__matching_ids and global_scene_tesselate_level > 0 :
            new_trimesh = {}
            for id_ele, triangles in self.__complete_trimesh.items() :
                new_tr_scene = []
                for t in triangles :
                    level = 0
                    iterate_triangles(t, level, global_scene_tesselate_level, 
                                                                new_tr_scene)
                new_trimesh[id_ele] = new_tr_scene
            self.__complete_trimesh = new_trimesh

        self.__pmin, \
        self.__pmax = compute_minmax_coord(self.__complete_trimesh)

        self.__triangleLmax = compute_trilenght_max(self.__complete_trimesh)

        # triangulation finale pour CARIBU
        if self.__lightmodel == "caribu":
            if legume_grid :
                raise ValueError("Conversion from voxels grid to triangles \
                                is not possible yet")
        
        # création d'une grille de voxels + distribution d'angles
        elif self.__lightmodel == "ratp":
            # triangulation non vide
            if self.__matching_ids :
                # création d'une nouvelle espèce pour les tiges
                if "stems id" in self.__geometry :
                    manage_stems_for_ratp(self.__geometry["stems id"], 
                                            self.__matching_ids, 
                                            self.__lightmodel_parameters)
                else : self.__geometry["stems id"] = None

                arg = (self.__complete_trimesh,        
                        [self.__pmin, self.__pmax],     
                        self.__triangleLmax,             
                        self.__matching_ids,              
                        self.__lightmodel_parameters,      
                        self.__environment["coordinates"], 
                        self.__environment["reflected"], 
                        self.__environment["infinite"],      
                        self.__geometry["stems id"],          
                        len(self.__geometry["scenes"]))
    
                self.__complete_voxmesh, \
                self.__matching_tri_vox, \
                self.__angle_distrib =  build_RATPscene_from_trimesh(*arg)

            # création d'une entrée géo RATP vide
            else :
                arg = (self.__lightmodel_parameters,      
                        [self.__pmin, self.__pmax],    
                        self.__environment["coordinates"],  
                        self.__environment["infinite"])

                self.__complete_voxmesh, \
                self.__angle_distrib =  build_RATPscene_empty(*arg)
                self.__matching_tri_vox = {}
            
            # si il y a une grille de voxels en entrée
            # et qu'elle est unique on la convertit en grille RATP
            if legume_grid and not self.__matching_ids :
                arg = (self.__geometry["scenes"][0],
                        self.__lightmodel_parameters,       
                        self.__environment["coordinates"],  
                        self.__environment["reflected"], 
                        self.__environment["infinite"])

                self.__complete_voxmesh, \
                self.__angle_distrib, \
                self.__nb0 =  legumescene_to_RATPscene(*arg)
                    
    def run(self, 
            energy=0, 
            day=0, 
            hour=0, 
            parunit="micromol.m-2.s-1", 
            truesolartime=False,
            id_sensors=None) :
        self.__energy = energy
        
        ## RATP ##
        if self.__lightmodel == "ratp" :
            vegetation = RATP_vegetation(self.__lightmodel_parameters, 
                                            self.__angle_distrib, 
                                            self.__environment["reflected"])
            meteo = RATP_meteo(energy,    
                                day, 
                                hour, 
                                self.__environment["coordinates"],
                                parunit, 
                                truesolartime,
                                self.__environment["direct"],
                                self.__environment["diffus"])

            if self.__complete_voxmesh.nveg > 0 :
                # calcul de RATP
                start=time.time()
                res = runRATP.DoIrradiation(self.__complete_voxmesh, 
                                            vegetation, 
                                            self.__sky, 
                                            meteo)
                self.__time_runmodel = time.time() - start

                self.__voxels_outputs = out_ratp_voxels(
                                                    self.__complete_voxmesh,
                                                    res,
                                                    parunit)

                # triangulation parmi les scènes d'entrée
                if self.__matching_ids :
                    arg = (self.__complete_trimesh,
                                self.__matching_ids,
                                self.__matching_tri_vox,
                                self.__voxels_outputs)
                    self.__triangles_outputs = out_ratp_triangles(*arg)

                    arg = (self.__matching_ids,
                    self.__environment["reflected"],
                    self.__lightmodel_parameters["reflectance coefficients"],
                    self.__triangles_outputs)
                    self.__elements_outputs = out_ratp_elements(*arg)

            else :
                print("--- Empty RATP grid")
                self.__voxels_outputs = out_ratp_empty_grid(day, hour)

        ## CARIBU ##
        elif self.__lightmodel == "caribu" :
            sun_up = False
            if self.__environment["direct"] :
                # calcul du soleil
                arg = (day, 
                        hour, 
                        self.__environment["coordinates"], 
                        truesolartime)
                if self.__lightmodel_parameters["sun algo"]=="ratp":
                    self.__sun = ratp_sun(*arg)
                elif self.__lightmodel_parameters["sun algo"]=="caribu":
                    self.__sun = caribu_sun(*arg)
                else :
                    raise ValueError("sun algo not recognize")
            
                # vérifie que le soleil soit levé
                # critère élévation min fixé à 2° comme dans RATP 
                # (shortwave_balance.F90 -> DirectBeam_Interception, l.68)
                sun_up = math.asin(-self.__sun[2])*(180/math.pi) > 2.
            
            compute = (sun_up and self.__environment["direct"]) or \
                    self.__environment["diffus"]
            if compute:             
                # préparatifs de CARIBU
                arg = (self.__complete_trimesh, 
                        self.__geometry, 
                        self.__matching_ids, 
                        (self.__pmin, self.__pmax), 
                        self.__lightmodel_parameters,
                        self.__environment["infinite"], 
                        id_sensors)
                opt, sensors_caribu, debug = Prepare_CARIBU(*arg)
                issensors = "sensors" in self.__lightmodel_parameters
                issoilmesh = self.__lightmodel_parameters["soil mesh"] != -1
                
                # init de la scène CARIBU
                if self.__environment["diffus"]  and \
                    self.__environment["direct"] :
                    c_scene_sky = CaribuScene(
                        scene=self.__complete_trimesh, 
                        light=self.__sky, 
                        opt=opt, 
                        scene_unit=self.__main_unit, 
                        pattern=self.__geometry["domain"], 
                        soil_mesh = self.__lightmodel_parameters["soil mesh"], 
                        debug = debug
                        )
                    sun_sky_option = "mix"
                    light = [tuple((1., self.__sun))]

                else :
                    if self.__environment["diffus"] :
                        light = self.__sky
                        sun_sky_option = "sky"
                    
                    elif self.__environment["direct"] :
                        light = [tuple((1., self.__sun))]
                        sun_sky_option = "sun"

                    else :
                        raise ValueError("Error with radiative inputs")

                c_scene = CaribuScene(
                    scene=self.__complete_trimesh, 
                    light=light, 
                    opt=opt, 
                    scene_unit=self.__main_unit, 
                    pattern=self.__geometry["domain"], 
                    soil_mesh = self.__lightmodel_parameters["soil mesh"], 
                    debug = debug
                    )
                
                # lancement
                arg = [c_scene, \
                        not self.__environment["reflected"], \
                        self.__environment["infinite"], \
                        sensors_caribu]
                if sun_sky_option == "mix":
                    start = time.time()
                    raw_sun, aggregated_sun = run_caribu(*arg)
                    arg[0] = c_scene_sky
                    raw_sky, aggregated_sky = run_caribu(*arg)
                    self.__time_runmodel = time.time() - start

                    #: Spitters's model estimating for the diffuse:direct ratio
                    # % de sky dans la valeur d'énergie finale
                    Rg = energy / 2.02  #: Global Radiation (W.m-2)
                    #: Diffuse fraction of the global irradiance
                    rdrs = spitters_horaire.RdRsH(
                                Rg=Rg, 
                                DOY=day, 
                                heureTU=hour, 
                                latitude=self.__environment["coordinates"][0])

                    raw,                    \
                    aggregated,             \
                    self.__sensors_outputs, \
                    self.__soilenergy = out_caribu_mix(rdrs,
                                                        c_scene,
                                                        c_scene_sky,
                                                        raw_sun, 
                                                        aggregated_sun, 
                                                        raw_sky, 
                                                        aggregated_sky,
                                                        issensors,
                                                        issoilmesh)
                else :
                    start = time.time()
                    raw, aggregated = run_caribu(*arg)
                    self.__time_runmodel = time.time() - start

                    self.__sensors_outputs, \
                    self.__soilenergy = out_caribu_nomix(c_scene,
                                                            aggregated,
                                                            issensors,
                                                            issoilmesh)
                    
            # gestion des sorties
            arg = [day,
                hour,
                self.__complete_trimesh,
                self.__matching_ids,
                aggregated,
                compute]
            self.__elements_outputs = out_caribu_elements(*arg)
            arg[4] = raw
            self.__triangles_outputs = out_caribu_triangles(*arg)
            
        ## RiRi (l-egume) ##
        elif self.__lightmodel == "riri" :
            dS = self.__lightmodel_parameters["voxel size"][0] * \
                                self.__lightmodel_parameters["voxel size"][1]
            
            tag_light_inputs = [self.__geometry["scenes"][0]["LA"] / dS, 
                                self.__geometry["scenes"][0]["triplets"], 
                                self.__geometry["scenes"][0]["distrib"], 
                                energy * dS]

            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            start = time.time()
            self.__legume_transmitted_light, \
            self.__legume_intercepted_light = \
                        riri.calc_extinc_allray_multi_reduced(
                                        *tag_light_inputs, 
                                        optsky=self.__in_environment["sky"][0], 
                                        opt=self.__in_environment["sky"][1]
                                                                )
            self.__time_runmodel = time.time() - start
    
    ## TRANSFER METHODS ##
    def to_MTG(self, energy=1, mtg=None, id=None) :
        # crée un tableau comme dans caribu_facade de fspm-wheat
        dico_par = {}
        para_dic = {}
        erel_dic = {}
        if id is None :
            for s in self.__elements_outputs["Organ"]:
                d = self.__elements_outputs[self.__elements_outputs.Organ==s]
                para_dic[s] = d["par Eabs"].values[0] * energy
                erel_dic[s] = d["par Eabs"].values[0] / self.__energy
        
        elif type(id) == list or type(id) == tuple :
            for esp in id:
                df_outputs_esp = self.__elements_outputs[ \
                                self.__elements_outputs.VegetationType==esp]
                for s in df_outputs_esp["Organ"]:
                    d = df_outputs_esp[df_outputs_esp.Organ==s]
                    para_dic[s] = d["par Eabs"].values[0] * energy
                    erel_dic[s] = d["par Eabs"].values[0] / self.__energy
        dico_par["PARa"] = para_dic
        dico_par["Erel"] = erel_dic

        for param in dico_par:
            if param not in mtg.properties():
                mtg.add_property(param)
            # update the MTG
            mtg.property(param).update(dico_par[param])

    def to_l_egume(self, 
                    energy=1, 
                    m_lais = [], 
                    list_lstring = [], 
                    list_dicFeuilBilanR = [], 
                    list_invar = [], 
                    id=None) :
                    
        epsilon = 1e-14

        if self.__lightmodel == "ratp" :
            return transfer_ratp_legume(m_lais, 
                                        energy, 
                                        self.__complete_voxmesh,
                                        self.__voxels_outputs,
                                        self.__nb0,
                                        epsilon)

        elif self.__lightmodel == "caribu" :
            skylayer = reduce_layers_from_trimesh(
                                self.__complete_trimesh, 
                                self.__pmax,
                                self.__lightmodel_parameters["sensors"][1],
                                self.__lightmodel_parameters["sensors"][2],
                                self.__matching_ids, 
                                id)
            
            return transfer_caribu_legume(energy,
                                skylayer, 
                                self.__elements_outputs, 
                                self.__sensors_outputs, 
                                self.__lightmodel_parameters["sensors"][1], 
                                self.__lightmodel_parameters["sensors"][2], 
                                m_lais, 
                                list_invar, 
                                list_lstring,
                                list_dicFeuilBilanR, 
                                self.__environment["infinite"],
                                epsilon) 

        else :
            raise ValueError("Unknown light model (ratp or caribu)")

    ## EXTERN TOOLS ##
    def s5(self):
        '''construit les fichiers d'entrée pour s5 et l'exécute
        '''
        currentfolder = os.path.abspath( \
                                    os.path.dirname(os.path.dirname(__file__)))
        s5folder = os.path.join(currentfolder, os.path.normpath("s5"))
        fort51 = os.path.join(s5folder, os.path.normpath("fort.51"))
        s5par = os.path.join(s5folder, os.path.normpath("s5.par"))

        # écriture du fichier fort.51 contenant la triangulation
        f=open(fort51, 'w')
        c_tr=1
        for id, triangles in self.__complete_trimesh.items() :
            for t in triangles :
                if tuple(self.__matching_ids[id]) in \
                                                self.__geometry["stems id"]:
                    stem='000'
                else:
                    stem='001'
                label = str(self.__matching_ids[id][1] + 1) +         \
                        str("%05i"%(self.__matching_ids[id][1] + 1)) + \
                        stem + '000'
                
                f.write("p\t1\t%s\t3\t"%(label))
                for i in range(3):
                    f.write("%f\t%f\t%f\t"%(t[i][0],t[i][1],t[i][2]))
                f.write("\n")
                c_tr += 1
        f.close()

        # écriture du fichier s5.par contenant les informations de la grille
        f=open(s5par, 'w')    
        f.write("%i\t%i\t%i\t%i\t0\n"%(self.__complete_voxmesh.nent, 
                                            9, 9, 
                                                self.__complete_voxmesh.njz))
        for i in range(self.__complete_voxmesh.njz):
            f.write("%f\t"%(self.__complete_voxmesh.dz[i]))
        f.write("\n")

        njx = self.__complete_voxmesh.njx
        njy = self.__complete_voxmesh.njy
        dx = self.__complete_voxmesh.dx
        dy = self.__complete_voxmesh.dy
        f.write("%f\t%i\t%f\t%f\t%i\t%f\n"%(njx*dx, njx, dx, njy*dy, njy, dy))
        f.close()

        # exécution de s5 dans un sous process
        subprocess.call(".\s5.exe", shell=True, cwd=s5folder)
        
        print("\n"+"--- Fin de s5.f")

    def s2v(self):
        '''construit les fichiers d'entrée pour s2v et l'exécute
        '''
        currentfolder = os.path.abspath( \
                                    os.path.dirname(os.path.dirname(__file__)))
        s2vfolder = os.path.join(currentfolder, os.path.normpath("s2v"))
        fort51 = os.path.join(s2vfolder, os.path.normpath("fort.51"))
        s2vpar = os.path.join(s2vfolder, os.path.normpath("s2v.par"))

        # écriture du fichier fort.51 contenant la triangulation
        f=open(fort51, 'w')
        c_tr=1
        for id, triangles in self.__complete_trimesh.items() :
            for t in triangles :
                if tuple(self.__matching_ids[id]) in \
                                                self.__geometry["stems id"]:
                    stem='000'
                else:
                    stem='001'
                label = str(self.__matching_ids[id][1] + 1) +         \
                        str("%05i"%(self.__matching_ids[id][1] + 1)) + \
                        stem + '000'
                
                f.write("p\t1\t%s\t3\t"%(label))
                for i in range(3):
                    f.write("%f\t%f\t%f\t"%(t[i][0],t[i][1],t[i][2]))
                f.write("\n")
                c_tr += 1
        f.close()

        f=open(s2vpar, 'w')
        # ligne 1 
        f.write("9 9 %i\n"%(self.__complete_voxmesh.njz))

        # ligne 2
        for i in range(self.__complete_voxmesh.njz):
            f.write("%f "%(self.__complete_voxmesh.dz[i]))
        f.write("\n")
        
        # ligne 3
        njx = self.__complete_voxmesh.njx
        njy = self.__complete_voxmesh.njy
        dx = self.__complete_voxmesh.dx
        dy = self.__complete_voxmesh.dy
        f.write("%f %i %f %f %i %f %i\n"%(njx*dx, 
                                            njx, 
                                            dx, 
                                            njy*dy, 
                                            njy, 
                                            dy, 
                                            self.__complete_voxmesh.nent))

        # ligne 4
        f.write("%f %f %f\n"%(self.__complete_voxmesh.xorig, 
                                self.__complete_voxmesh.yorig, 
                                    self.__complete_voxmesh.zorig))

        f.close()

        # exécution de s5 dans un sous process
        subprocess.call(".\s2v.exe", shell=True, cwd=s2vfolder)

        print("--- Fin de s2v.cpp")

    ## VTK FILES ##
    def VTK_nolight(self, path, i=None, printtriangles=True, printvoxels=True) :
        '''construit des fichiers VTK de la triangulation et de la grille de
         voxels après leur construction

        '''       
        if self.__lightmodel == "ratp" and printvoxels: 
            if i is None : filepath = path + "_voxels_nolight.vtk"
            else : filepath = path + "_voxels_nolight" + "_" + str(i) + ".vtk"
            ratp_prepareVTK(self.__complete_voxmesh, filepath)
            
        if self.__matching_ids and printtriangles:
            if i is None : filepath = path + "_triangles_nolight.vtk"
            else : filepath = path + "_triangles_nolight" + "_" + str(i)+".vtk"
            VTKtriangles(self.__complete_trimesh, [], [], filepath)

    def VTK_light(self, path, i=None, printtriangles=True, printvoxels=True) :
        if (self.__lightmodel == "ratp" and \
                    (not hasattr(self, '_LightVegeManager__voxels_outputs'))) :
            print("--- VTK:  No light data, run the simulation")
            self.VTK_nolight(path, i, printtriangles, printvoxels)
        elif (self.__lightmodel == "caribu" and \
                    (not hasattr(self, '_LightVegeManager__triangles_outputs'))) :
            print("--- VTK:  No light data, run the simulation")
            self.VTK_nolight(path, i, printtriangles, printvoxels)
        else :
            if self.__lightmodel == "ratp" and printvoxels: 
                if i is None : filepath = path + "_voxels_light.vtk"
                else : filepath = path + "_voxels_light" + "_" + str(i) + ".vtk"

                datanames = ['ShadedPAR', 
                                'SunlitPAR', 
                                'ShadedArea',
                                'SunlitArea',
                                'PARa', 
                                'Intercepted', 
                                'Transmitted']

                ratp_prepareVTK(self.__complete_voxmesh, 
                                    filepath, 
                                    datanames, 
                                    self.__voxels_outputs)
                
            if self.__matching_ids and printtriangles:
                if i is None : filepath = path + "_triangles_light.vtk"
                else : filepath = path + "_triangles_light" + "_" + str(i) + ".vtk"

                datanames = []
                if self.__lightmodel == "ratp" :
                    datanames = ['ShadedPAR', 
                                'SunlitPAR', 
                                'ShadedArea',
                                'SunlitArea',
                                'PARa', 
                                'Intercepted', 
                                'Transmitted']
                elif self.__lightmodel == "caribu" :
                    datanames = []
                    for band in self.__lightmodel_parameters["caribu opt"].keys() :
                        datanames.append(band + " Eabs")
                        datanames.append(band + " Ei")

                data = [list(self.__triangles_outputs[name]) for name in datanames]
                VTKtriangles(self.__complete_trimesh, data, datanames, filepath)

    def VTK_sun(self, path, scale=2, orig=(0,0,0), center=True, i=None) :
        if i is None : filepath = path + "_sun.vtk"
        else : filepath = path + "_sun" + "_" + str(i) + ".vtk"

        if center :
            start = tuple([a + b * scale for a,b in zip(orig, self.__sun)])
            end = tuple([a + b * -scale for a,b in zip(orig, self.__sun)])
        else :
            start = orig
            end = tuple([a + b * scale for a,b in zip(orig, self.__sun)])

        VTKline(start, end, filepath)
    
    ## GETTERS ##
    @property
    def legume_transmitted_light(self):
        return self.__legume_transmitted_light 
    
    @property
    def legume_intercepted_light(self):
        return self.__legume_intercepted_light 

    @property
    def elements_outputs(self):
        return self.__elements_outputs 

    @property
    def triangles_outputs(self):
        return self.__triangles_outputs 
    
    @property
    def voxels_outputs(self):
        return self.__voxels_outputs 

    @property
    def sensors_outputs(self):
        try:
            return self.__sensors_outputs
        except AttributeError:
            return None

    @property
    def sun(self):
        return self.__sun

    @property
    def soilenergy(self):
        return self.__soilenergy

    @property
    def maxtrianglearea(self):
        try :
            return self.__areamax
        except AttributeError:
            return 0.
  
    @property
    def legume_empty_layers(self):
        return self.__nb0

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
        except AttributeError :
            return 0.
