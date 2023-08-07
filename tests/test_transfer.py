from lightvegemanager.transfer import transfer_ratp_legume, transfer_caribu_legume
import numpy
import pandas
from alinea.pyratp import grid


def test_transfer_ratp_legume():
    m_lais = numpy.zeros([2, 10, 4, 4])
    m_lais[0][8][0][0] = 10.0
    m_lais[0][8][0][1] = 12.0
    m_lais[1][9][1][2] = 14.0
    m_lais[1][9][3][3] = 16.0
    energy = 500.0
    pars = {"latitude": 0, "longitude": 0, "timezone": 0, "nent": 2, "rs": (0, 0), "idecaly": 0, "orientation": 0}
    pars.update(
        {"njx": 4, "njy": 4, "njz": 2, "dx": 1.5, "dy": 1.5, "dz": [1.0] * 2, "xorig": 0, "yorig": 0, "zorig": 0}
    )
    mygrid = grid.Grid.initialise(**pars)

    # we fill empty voxels in canopy with LA = 1e-14 to compute them in RATP
    voxels_outputs = pandas.DataFrame(
        {
            "VegetationType": [1, 1, 2, 2, 2, 2, 1, 1],
            "Nx": [4, 4, 3, 1, 4, 4, 3, 1],
            "Ny": [1, 2, 3, 4, 1, 2, 3, 4],
            "Nz": [1, 1, 2, 2, 1, 1, 2, 2],
            "Intercepted": [0.5, 0.4, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0],
            "Transmitted": [0.5, 0.6, 0.7, 0.8, 1.0, 1.0, 1.0, 1.0],
        }
    )
    nb0 = m_lais.shape[1] - mygrid.njz
    res_abs_i, res_trans = transfer_ratp_legume(m_lais, energy, mygrid, voxels_outputs, nb0)

    dS = mygrid.dx * mygrid.dy
    assert res_trans[0][0][0] == dS * energy
    assert res_trans[8][0][0] == (1.0 + 0.5) * energy
    assert res_trans[8][0][1] == 1.6 * energy
    assert res_trans[9][1][2] == 1.7 * energy
    assert res_trans[9][3][3] == 1.8 * energy

    assert res_abs_i[0][0][0][0] == 0.0  # fill with 0. above the canopy
    assert res_abs_i[1][0][0][0] == 0.0
    assert res_abs_i[1][8][0][0] == 1e-8  # > 0. for empty voxels inside the canopy
    assert res_abs_i[0][8][0][0] == 0.5 * energy
    assert res_abs_i[0][8][0][1] == 0.4 * energy
    assert res_abs_i[1][9][1][2] == 0.3 * energy
    assert res_abs_i[1][9][3][3] == 0.2 * energy


def test_transfer_caribu_legume():
    energy = 500.0
    skylayer = 0
    id = None
    elements_outputs = pandas.DataFrame(
        {
            "VegetationType": [0, 0, 0, 1, 1],
            "Organ": [1, 2, 22, 3, 4],
            "Area": [10.0, 12.0, 8.0, 14.0, 16.0],
            "par Ei": [0.5, 0.4, 0.9, 0.3, 0.2],
        }
    )
    sensors_outputs = {"par": {}}
    for i in range(32):
        sensors_outputs["par"][i] = 0.4

    sensors_dxyz = (1.5, 1.5, 1.0)
    sensors_nxyz = (4, 4, 2)
    m_lais = numpy.zeros([2, 10, 4, 4])
    m_lais[0][8][0][0] = 10.0
    m_lais[0][8][0][1] = 12.0
    m_lais[1][9][1][2] = 14.0
    m_lais[1][9][3][3] = 16.0
    list_invar = [{"Hplante": [0.0] * 2}, {"Hplante": [0.0] * 1}]
    list_lstring = [
        {
            1: [0, 0, 0, 0, 0, 0, 0, 0, 0, "dev"],
            2: [0, 0, 0, 0, 0, 0, 0, 0, 0, "sen"],
            22: [1, 0, 0, 0, 0, 0, 0, 0, 0, "dev"],
        },
        {3: [0, 0, 0, 0, 0, 0, 0, 0, 0, "sen"], 4: [0, 0, 0, 0, 0, 0, 0, 0, 0, "dev"]},
    ]
    list_dicFeuilBilanR = [
        {"surf": [0.5, 1e-8]},
        {"surf": [0.5]},
    ]
    infinite = True
    epsilon = 1e-8

    res_trans = transfer_caribu_legume(
        energy,
        skylayer,
        id,
        elements_outputs,
        sensors_outputs,
        sensors_dxyz,
        sensors_nxyz,
        m_lais,
        list_invar,
        list_lstring,
        list_dicFeuilBilanR,
        infinite,
        epsilon,
    )

    dS = sensors_dxyz[0] * sensors_dxyz[1]
    assert res_trans[0][0][0] == dS * energy
    assert res_trans[9][0][0] == 0.4 * dS * energy
    numpy.testing.assert_array_equal(list_invar[0]["parap"], [216.0, 311.04])
    numpy.testing.assert_array_equal(list_invar[0]["parip"], [423.36, 311.04])
    assert list_invar[1]["parap"] == 138.24
    assert list_invar[1]["parip"] == 319.68
