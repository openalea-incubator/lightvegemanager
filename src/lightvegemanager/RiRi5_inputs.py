"""
    RiRi5inputs
    **********

    Translate RATP grid in RiRi5 format
"""
import numpy

def ratpgrid_to_riri5(ratpgrid):
    la = numpy.zeros([ratpgrid.nent, ratpgrid.njz + 1, ratpgrid.njy, ratpgrid.njx])
    for ne in range(ratpgrid.nent):
        for ix in range(ratpgrid.njx):
            for iy in range(ratpgrid.njy):
                for iz in range(ratpgrid.njz):
                    k = ratpgrid.kxyz[ix, iy, iz]
                    if k > 0 :
                        la[ne][iz+1][ratpgrid.njy - (iy+1)][ix] = ratpgrid.leafareadensity[ne, k]

    for ne in range(ratpgrid.nent):
        zeros_array = numpy.zeros([ratpgrid.njy, ratpgrid.njx])
        la[ne][0][:][:] = zeros_array

    return la
