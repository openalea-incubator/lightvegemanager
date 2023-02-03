'''
Ecriture de fichiers VTK visualisation
'''

import itertools
import numpy
import string

from alinea.caribu import plantgl_adaptor

def VTKtriangles(trimesh, var, varname, filename):
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
    # correct and remove spaces from varname
    varname = [n.replace(" ", "_") for n in varname]

    if isinstance(trimesh, dict) :
        list_triangles = list(itertools.chain(*[v for v in trimesh.values()]))
        ndefaultfiedl = 2
    elif isinstance(trimesh, list) : 
        list_triangles = trimesh
        ndefaultfiedl = 1

    nbtr = len(list_triangles)
    
    f = open(filename, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS '+str(nbtr * 3)+' float\n')

    for tr in list_triangles:
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
    f.write('FIELD FieldData '+str(len(varname) + ndefaultfiedl)+'\n')

    # input fields
    for i, name in enumerate(varname):
        f.write(name+' 1 '+str(nbtr)+' float\n')
        for j in range(nbtr):
            f.write(str(var[i][j])+'\n')
        f.write('\n')
    f.write('\n')

    # ID field : 0 to len(triangles)
    f.write("ID 1 "+str(nbtr)+' float\n')
    for i in range(nbtr):
        f.write(str(i)+'\n')
    f.write('\n')

    # Element ID
    if ndefaultfiedl == 2 :
        f.write("Organ 1 "+str(nbtr)+' float\n')
        for id, triangles in trimesh.items():
            for t in triangles :
                f.write(str(id)+'\n')
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

def ratp_prepareVTK(ratpgrid, filename, columns=[], df_values=None) :
    # par défaut imprime le LAD
    kxyz, entities, lad = [], [], []
    for i in range(ratpgrid.nveg)  :
        for j in range(ratpgrid.nje[i]):
            lad.append(ratpgrid.leafareadensity[j, i])
            entities.append(ratpgrid.nume[j,i])
            kxyz.append(int(i)+1)
    datafields = [numpy.array(kxyz), numpy.array(entities), numpy.array(lad)]

    # autres variables
    if columns  and df_values is not None:
        for name in columns :
            v = []
            for i in range(ratpgrid.nveg)  :
                for j in range(ratpgrid.nje[i]):
                    filter = (df_values.Voxel == i+1) & (df_values.VegetationType == j+1)
                    v.append(df_values[filter][name].values[0])
            datafields.append(numpy.array(v))
    
    columns.insert(0, "LAD")
    # correct and remove spaces from columns
    columns = [n.replace(" ", "_") for n in columns]

    VTKvoxels(ratpgrid, datafields, columns, filename)

def VTKvoxels(grid, datafields, varnames, nomfich):
    '''    Display Voxels colored by variable with Paraview
           RATP Grid is written in VTK Format as a structured grid
           Inputs: ... variable : a list of 3 arrays composed of the a 
           RATP variable to be plotted, corresponding entities, and Voxel ID
                   ... grid : the RATP grid
                   ... varname: name of the variable to be plotted
                   ... nomfich: the VTK filename and path
           Outputs: ... a VTK file

           datafields :  [0] : kxyz
                        [1] : entities
                        [2, ..., n] : variable_0, ... variable_n-2
    '''

    f=open(nomfich,'w')

    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET RECTILINEAR_GRID\n')
    f.write('DIMENSIONS '+str(grid.njx+1)+' '+str(grid.njy+1)+' '+str(grid.njz+2)+'\n')

    f.write('Z_COORDINATES '+str(grid.njz+2)+' float\n')
    for i in range(grid.njz):
       z = grid.dz[i:grid.njz].sum()-grid.zorig
       f.write(str(z)+' ')
    f.write(str(-grid.zorig)+' ')
    f.write(str(-10*grid.dz[0]-grid.zorig)+' ')
    f.write('\n')

    f.write('Y_COORDINATES '+str(grid.njy+1)+' float\n')
    for i in range(grid.njy, -1, -1):
       y = grid.dy*i+grid.yorig
       f.write(str(y)+' ')
    f.write('\n')

    f.write('X_COORDINATES '+str(grid.njx+1)+' float\n')
    for i in  range(grid.njx+1):
       x = grid.dx*i+grid.xorig
       f.write(str(x)+' ')
    f.write('\n')

    numVoxels = (grid.njx)*(grid.njy)*(grid.njz+1)
    f.write('CELL_DATA '+str(numVoxels)+'\n')

    # Set the number of entities to write - NbScalars
    numberofentities = len(set(datafields[1]))
    numberofdatafields = len(datafields) - 2

    for ent in range(numberofentities) :
        for nvar in range(numberofdatafields) : 
            f.write('SCALARS ' + varnames[nvar] + '_entity_' + \
                                            str(int(ent)) + ' float  1 \n')
            f.write('LOOKUP_TABLE default\n')
            for ik in range(grid.njz+1):
                for ij in range(grid.njy):
                    for ii in range(grid.njx): 
                        k =grid.kxyz[ii,ij,ik]
                        if (k>0): 
                            kidDummy = numpy.where(numpy.array(datafields[0])==k)
                            kid = kidDummy[0]
                            entity = numpy.array(datafields[1])[kid]
                            ent_in_vox =numpy.where(entity == ent)
                            # pas dans le voxel
                            if numpy.alen(ent_in_vox[0])<1:
                                f.write(str(-9999.0)+'\n')
                            else:
                                # couche du sol
                                if ik == grid.njz:
                                    f.write(str(-9999.0)+'\n')
                                else:
                                    where = (kid[numpy.where(entity==ent)])
                                    value = datafields[2+nvar][where[0]]
                                    f.write(str(value)+'\n')
                        # voxel vide
                        else:
                            f.write(str(-9999.0)+'\n')
            f.write('\n')

        # indice des voxels
        f.write('SCALARS Voxel_ID float  1 \n')
        f.write('LOOKUP_TABLE default\n')
        for ik in range(grid.njz+1):
            for ij in range(grid.njy):
                for ii in range(grid.njx):
                    k =grid.kxyz[ii,ij,ik]
                    f.write(str(k)+'\n')
        f.write('\n')
    f.close()


def PlantGL_to_VTK(scenes, path, i=0, in_unit="m", out_unit="m"):
    units = {'mm': 0.001, 
                'cm': 0.01, 
                'dm': 0.1, 
                'm': 1, 
                'dam': 10, 
                'hm': 100,
                'km': 1000}
    if in_unit not in units or out_unit not in units:
        raise ValueError("Unknown unit: select one in this \
                    list ['mm','cm','dm','m','dam','hm','km']")
    
    rescale=False
    if (in_unit != out_unit) : rescale=True
    
    trimesh = {}
    # si on a plusieurs scènes
    if type(scenes) == list :
        for s in scenes :   
            trimesh = trimesh + plantgl_adaptor.scene_to_cscene(s)
    else :
       trimesh = plantgl_adaptor.scene_to_cscene(scenes)

    if rescale :
        h = units[in_unit]/units[out_unit]
        for id in trimesh.keys() : trimesh[id] = rescale(trimesh[id], h)

    VTKtriangles(trimesh, [], [], path+"triangles_plantgl_"+str(i)+".vtk")