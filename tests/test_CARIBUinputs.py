from lightvegemanager.CARIBUinputs import Prepare_CARIBU, CARIBU_opticals, create_caribu_legume_sensors, run_caribu

import openalea.plantgl.all as pgl


def test_CARIBU_optical():
    matching_ids = {0: (123, 0), 1: (111, 0), 2: (222, 1)}
    stems = [(111, 0)]
    parameters = {"caribu opt": {"par": (0.5, 0.2), "nir": (0.1, 0.8)}}

    opt = CARIBU_opticals(matching_ids, parameters, stems)

    assert opt["par"][0] == (0.5, 0.2)
    assert opt["par"][1] == (0.5,)
    assert opt["par"][2] == (0.5, 0.2)
    assert opt["nir"][0] == (0.1, 0.8)
    assert opt["nir"][1] == (0.1,)
    assert opt["nir"][2] == (0.1, 0.8)


def test_Prepare_CARIBU():
    matching_ids = {0: (123, 0), 1: (111, 0), 2: (222, 1)}
    geometry = {"stems id": [(111, 0)]}
    parameters = {
        "caribu opt": {"par": (0.5, 0.2), "nir": (0.1, 0.8)},
        "sensors": ["grid", (1.0, 1.0, 0.5), (3, 3, 3), (0.0, 0.0, 0.0)],
        "debug": True,
    }
    triangles = {
        0: [[(0.2, 0.1, 0.1), (0.7, 0.1, 0.1), (0.7, 0.7, 0.6)]],
        1: [[(1.3, 0.1, 0.1), (1.7, 0.1, 0.1), (1.7, 0.7, 0.6)]],
        2: [[(3.4, 0.1, 0.1), (3.7, 0.1, 0.1), (3.7, 0.7, 0.7)]],
    }
    minmax = [(0.2, 0.1, 0.1), (3.7, 0.7, 0.7)]
    infinite = True
    idsensors = [0]

    opt, sensors_caribu, debug, matching_sensors_species = Prepare_CARIBU(
        triangles, geometry, matching_ids, minmax, parameters, infinite, idsensors
    )
    assert opt["par"][0] == (0.5, 0.2)
    assert opt["par"][1] == (0.5,)
    assert opt["par"][2] == (0.5, 0.2)
    assert opt["nir"][0] == (0.1, 0.8)
    assert opt["nir"][1] == (0.1,)
    assert opt["nir"][2] == (0.1, 0.8)
    assert len(sensors_caribu.items()) == 18
    assert debug


def test_create_caribu_legume_sensors():
    dxyz = (1.0, 1.0, 0.5)
    nxyz = (3, 3, 4)
    orig = (0.0, 0.0, 0.0)

    triangles = {
        0: [[(0.2, 0.1, 0.1), (0.7, 0.1, 0.1), (0.7, 0.7, 0.6)]],
        1: [[(1.3, 0.1, 0.1), (1.7, 0.1, 0.1), (1.7, 0.7, 0.6)]],
        2: [[(3.4, 0.1, 0.1), (3.7, 0.1, 0.1), (3.7, 0.7, 0.7)]],
    }
    minmax = [(0.2, 0.1, 0.1), (3.7, 0.7, 0.7)]
    matching_ids = {0: (123, 0), 1: (111, 0), 2: (222, 1)}
    id_sensors = [0]
    infinite = True

    caribu_sensors, s_capt, sensors_maxcenter = create_caribu_legume_sensors(
        dxyz, nxyz, orig, minmax[1], triangles, matching_ids, id_sensors, infinite
    )

    # only two layers of sensors are created as we focus entity 1
    assert isinstance(s_capt, pgl.Scene)
    assert len(s_capt) == 18

    assert isinstance(caribu_sensors, dict)
    assert isinstance(caribu_sensors[0], list)
    assert len(caribu_sensors.items()) == 18
    assert caribu_sensors[0][0] == [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)]
    assert caribu_sensors[17][1] == [(2.0, 2.0, 0.5), (3.0, 3.0, 0.5), (2.0, 3.0, 0.5)]

    assert sensors_maxcenter == [1.5, 1.5, 0.5]
