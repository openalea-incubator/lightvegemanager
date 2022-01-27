from src.Polygons import *

'''
Set de fonctions pour subdiviser un triangle de manière à ce qu'il match une grille de voxels
'''

def whichvoxel(p, mygrid):
    '''Renvoit le voxel (indices 3D) auquel appartient p

    p : Vector3, point à tester
    mygrid : grille de voxels façon RATP

    return : [kx, ky, kz]
    '''
    vox = [(int)(p[0]/mygrid.dx)]
    vox.append((int)(p[1]/mygrid.dy))

    # recherche de le couche en z
    jz=-10
    arret=0
    j=1
    while((j < mygrid.njz) and (arret==0)):
        if(p[2]>mygrid.dz[j]):
            arret=1
        j+=1
    if arret==1:
        jz=j-1
    else:
        jz=mygrid.njz-1
    vox.append(jz)
    return vox

def samevoxel(voxels):
    '''Vérifie si la liste voxels contient le même voxel

    voxels : [[kx1, ky1, kz1], ..., [kxn, kyn, kzn]] listes de voxels en indices 3D
    
    return : boolean 
    '''
    # test si tous les sommets sont dans le même triangle
    test=0
    for i in range(len(voxels)):
        for j in range(3):
            test+=abs((float)(voxels[i][j]-voxels[(i+1)%3][j]))

    # triangle dans voxel
    return test == 0

def tesselate(mygrid, triangle):
    '''Subdivise triangle en 4 sous-triangles SI il appartient à plusieurs voxels
        SINON renvoit [triangle]

    mygrid : grille de voxels façon RATP
    triangle : Triangle3
    
    return : [Triangle3] liste de 1 triangle ou 4 triangles
    '''
    # teste si triangle dans voxel
    # voxel dans chaque sommet
    wh = []
    for i in range(3):
        wh.append(whichvoxel(triangle[i], mygrid))
    
    if samevoxel(wh):
        return [triangle]
    
    # on subdivise en 4 triangles
    else:
        smalltriangles=[]

        # récupère les milieux
        middles = []
        for i in range(3):
            middles.append(Vector3.middle(triangle[i], triangle[(i+1)%3]))

        # triangle 1
        smalltriangles.append(Triangle3(
            triangle[0],
            middles[0],
            middles[2]
        ))
        
        # triangle 2
        smalltriangles.append(Triangle3(
            middles[0],
            triangle[1],
            middles[1]
        ))

        # triangle 3
        smalltriangles.append(Triangle3(
            middles[1],
            triangle[2],
            middles[2]
        ))

        # triangle 4
        smalltriangles.append(Triangle3(
            middles[0],
            middles[1],
            middles[2]
        ))
        
        for t in smalltriangles : t.set_id(triangle.id)

        return smalltriangles

def iterate_trianglesingrid(triangle, mygrid, level, levelmax, triangles_shape):
    '''Récursion sur triangle, tant que sa subdivision ne match pas la grille ou atteint un critère d'arrêt

    triangle : Triangle3, triangle à tester 
    mygrid : grille de voxels façon RATP
    level : nombre de subdivision actuel
    levelmax : limite de subdivision
    triangles_shape : liste de Triangle3, stocke les triangles subdivisés ou non
    
    return : 1 si fonction a bien terminée et 
            met à jour la liste triangles_shape
    '''
    level+=1
    ltriangle = tesselate(mygrid, triangle)
    if len(ltriangle)==1:
        triangles_shape.append(ltriangle[0])
    else:
        if level == levelmax:
            triangles_shape.append(triangle)
        
        else:
            for subt in ltriangle:
                iterate_trianglesingrid(subt, mygrid, level, levelmax, triangles_shape)
    
    return 1
