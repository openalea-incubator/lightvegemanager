from lightvegemanager.RATPinputs import RATP_meteo, RATP_vegetation
import pytest

parameters_global = {"angle distrib algo": ""}
parameters_voxels = {"angle distrib algo": "compute voxel"}


@pytest.mark.parametrize(
    "test_input",
    [(parameters_global), (parameters_voxels)],
)
def test_RATP_vegetation(test_input):
    # voxel distrib needs global distrib
    angle_global = [[0.5, 0.5]]
    angle_voxels = [[[0.5, 0.5]], [[0.3, 0.7]]]  # dim : voxel, ent, classe d'angle
    angles = {"global": angle_global, "voxel": angle_voxels}
    reflected = True

    test_input.update({"reflectance coefficients": [[0.0]], "mu": [1.0]})

    ratp_vege = RATP_vegetation(test_input, angles, reflected)

    if test_input["angle distrib algo"] == "compute voxel":
        # fortran float are float32, conversion adds artifacts in ~10^-10
        assert pytest.approx(ratp_vege.distincvox[1][0][0]) == 0.3
    else:
        assert ratp_vege.distincvox is None

    assert ratp_vege.distinc[0][0] == 0.5
    assert ratp_vege.mu[0] == 1.0
    assert ratp_vege.rf[0, 0] == 0.0


test1 = ("micromol.m-2.s-1", True, True)  # unit, diffus, direct
result_1 = (108.6956, 61.9996)  # Rg, Rdif
test2 = ("micromol.m-2.s-1", True, False)
result_2 = (108.6956, 108.6956)  # Rg, Rdif
test3 = ("W.m-2", False, False)
result_3 = (500.0, 0.0)
test4 = ("W.m-2", False, True)
result_4 = (500.0, 0.0)


@pytest.mark.parametrize(
    "test_input, expected",
    [(test1, result_1), (test2, result_2), (test3, result_3), (test4, result_4)],
)
def test_RATP_meteo(test_input, expected):
    day = 300
    hour = 15
    latitude = 45
    energy = 500

    ratp_meteo = RATP_meteo(
        energy, day, hour, [latitude], test_input[0], True, diffus=test_input[1], direct=test_input[2]
    )

    assert pytest.approx(ratp_meteo.tabmeteo[0][2]) == expected[0]
    assert pytest.approx(ratp_meteo.tabmeteo[0][3]) == expected[1]
