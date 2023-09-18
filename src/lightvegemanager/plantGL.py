'''
    plantGL
    ******************************

    visualisation plantGL
'''
import itertools
import numpy

import openalea.plantgl.all as pgl

def cscene_to_plantGLScene_stems(cscene, stems_id=None, matching_ids={}):
    pglscene = pgl.Scene()
    for t in itertools.chain(*[v for k,v in cscene.items() if stems_id is None or tuple(matching_ids[k]) not in stems_id]):
        pts = []
        ind = []
        i = 0
        for p in t :
            pts.append(p)
        ind.append((i, i+1, i+2))
        i += 3

        ts = pgl.TriangleSet(pts, ind)
        mat = pgl.Material(ambient=(50, 205, 50), shininess=0.1, diffuse=1)
        pglscene.add(pgl.Shape(ts, mat))

    if stems_id is not None :
        for t in itertools.chain(*[v for k,v in cscene.items() if tuple(matching_ids[k]) in stems_id]):
            pts = []
            ind = []
            i = 0
            for p in t :
                pts.append(p)
            ind.append((i, i+1, i+2))
            i += 3

            ts = pgl.TriangleSet(pts, ind)
            mat = pgl.Material(ambient=(139, 69, 19), shininess=0.1, diffuse=1)
            pglscene.add(pgl.Shape(ts, mat))

    return pglscene

def cscene_to_plantGLScene_light(cscene, outputs={}, column_name="par Ei"):
    plt_cmap = "seismic"
    minvalue = numpy.min(outputs[column_name].values)
    maxvalue = numpy.max(outputs[column_name].values)
    colormap = pgl.PglMaterialMap(minvalue, maxvalue, plt_cmap)

    pglscene = pgl.Scene()

    for e,t in enumerate(itertools.chain(*[v for v in cscene.values()])):
        pts = []
        ind = []
        i = 0
        for p in t :
            pts.append(p)
        ind.append((i, i+1, i+2))
        i += 3

        ts = pgl.TriangleSet(pts, ind)
        value = outputs[outputs.Triangle == e][column_name].values[0]
        mat = colormap(value)
        mat.shininess = 0.1
        mat.diffuse = 1
        pglscene.add(pgl.Shape(ts, mat))

    return pglscene

def ratpgrid_to_plantGLScene(ratpgrid, transparency=0., plt_cmap="Greens", outputs={}):
    if plt_cmap == "Greens" :
        # minvalue = numpy.min(ratpgrid.s_vx)
        minvalue = 0.
        maxvalue = numpy.max(ratpgrid.s_vx)
    elif plt_cmap == "seismic" :
        minvalue = numpy.min(outputs["PARa"])
        maxvalue = numpy.max(outputs["PARa"])
    colormap = pgl.PglMaterialMap(minvalue, maxvalue, plt_cmap)
    
    scene = pgl.Scene()
    vsize = (float(ratpgrid.dx)/2, float(ratpgrid.dy)/2, float(ratpgrid.dz[0])/2)
    for z in range(ratpgrid.njz):
        for x in range(ratpgrid.njx):
            for y in range(ratpgrid.njy):
                k = ratpgrid.kxyz[x, y, z]
                if k > 0 :
                    if plt_cmap == "Greens" :
                        value = ratpgrid.s_vx[k-1]
                    elif plt_cmap == "seismic" :
                        value = outputs[outputs.Voxel==k]["PARa"].values[0]
                    
                    mat = colormap(value)
                    mat.transparency = transparency

                    vectrans = (float(ratpgrid.xorig + (0.5 + x) * ratpgrid.dx ), 
                                float(ratpgrid.yorig + (0.5 + y) * ratpgrid.dy ), 
                                float(ratpgrid.dz[z:ratpgrid.njz].sum()-ratpgrid.zorig) - 0.5 * ratpgrid.dz[z] )
                    shape = pgl.Shape(pgl.Translated(vectrans, pgl.Box(vsize)), mat)
                    scene.add(shape)
    return scene