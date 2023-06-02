from lightvegemanager.sky import *
import os
import pytest


ndir_turtle46 = 46
hmoy_turtle46 = [0.16109389, 0.16109389, 0.16109389, 0.16109389, 0.16109389, 0.16109389, 0.16109389, 0.16109389, 0.16109389, 0.16109389, 0.1886701, 0.1886701, 0.1886701, 0.1886701, 0.1886701, 0.46373397, 0.46373397, 0.46373397, 0.46373397, 0.46373397, 0.54244834, 0.54244834, 0.54244834, 0.54244834, 0.54244834, 0.54244834, 0.54244834, 0.54244834, 0.54244834, 0.54244834, 0.8274606, 0.8274606, 0.8274606, 0.8274606, 0.8274606, 0.91839224, 0.91839224, 0.91839224, 0.91839224, 0.91839224, 1.2070698, 1.2070698, 1.2070698, 1.2070698, 1.2070698, 1.5707964]
azmoy_turtle46 = [0.21345377, 1.0431833, 1.4700909, 2.2998204, 2.726728, 3.5564575, 3.983365, 4.8130946, 5.240002, 6.0697317, 0.62831855, 1.8849556, 3.1415927, 4.3982296, 5.6548667, 0.0, 1.2566371, 2.5132742, 3.7699113, 5.0265484, 0.40613812, 0.8504989, 1.6627752, 2.107136, 2.9194121, 3.363773, 4.176049, 4.62041, 5.4326863, 5.877047, 0.0, 1.2566371, 2.5132742, 3.7699113, 5.0265484, 0.62831855, 1.8849556, 3.1415927, 4.3982296, 5.6548667, 0.0, 1.2566371, 2.5132742, 3.7699113, 5.0265484, 3.1415927]
ndir_file = 5
hmoy_file = [0.46373397, 0.46373397, 0.46373397, 0.46373397, 1.5707964 ]
azmoy_file = [0., 1.5707964, 3.1415927, 4.712389 , 0. ]
ndir_list = 20
hmoy_list = [0.15707964, 0.4712389 , 0.7853982 , 1.0995574 , 1.4137167 ,
                0.15707964, 0.4712389 , 0.7853982 , 1.0995574 , 1.4137167 ,
                0.15707964, 0.4712389 , 0.7853982 , 1.0995574 , 1.4137167 ,
                0.15707964, 0.4712389 , 0.7853982 , 1.0995574 , 1.4137167 ]
azmoy_list = [0.7853982, 0.7853982, 0.7853982, 0.7853982, 0.7853982, 2.3561945,
                2.3561945, 2.3561945, 2.3561945, 2.3561945, 3.9269907, 3.9269907,
                3.9269907, 3.9269907, 3.9269907, 5.497787 , 5.497787 , 5.497787 ,
                5.497787 , 5.497787 ]

datafile = os.path.join(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data"), "sky_5.data")

@pytest.mark.parametrize("test_input, expected_1, expected_2, expected_3", [
                                            ("turtle46", ndir_turtle46, hmoy_turtle46, azmoy_turtle46),
                                            (datafile, ndir_file, hmoy_file, azmoy_file), 
                                            ([4, 5, "soc"], ndir_list, hmoy_list, azmoy_list)])
def test_RATPsky(test_input, expected_1, expected_2, expected_3):
    sky = RATPsky(test_input)
    assert sky.ndir == expected_1
    for x,y in zip(sky.hmoy, expected_2) : assert pytest.approx(x) == y
    for x,y in zip(sky.azmoy, expected_3) : assert pytest.approx(x) == y

def test_ratpformat_to_caribuformat():
    sky = ratpformat_to_caribuformat([0., 180., 0., 180.], [45., 45., 90., 90.], [0.]*4, False)
    root2_2 = 2**0.5/2
    expected = [[root2_2, 0., -root2_2], [-root2_2, 0., -root2_2], [0, 0, -1], [0, 0, -1]]
    for t, c in zip(sky, expected):
        for x, y in zip(t[1],c) : assert pytest.approx(x) == y

pc_soc = [0.2918418884019217, 0.20815811159807826, 0.2918418884019217, 0.20815811159807826]
pc_uoc = [0.25, 0.25, 0.25, 0.25]
@pytest.mark.parametrize("test_input, expected", [
                                            ("soc", pc_soc),
                                            ("uoc", pc_uoc)])
def test_discrete_sky(test_input, expected):
    ele_correct = [22.5, 67.5, 22.5, 67.5]
    az_correct = [90.0, 90.0, 270.0, 270.0]
    omega_correct = [1.5707963267948966, 1.5707963267948966, 1.5707963267948966, 1.5707963267948966]
    ele, azi, omega, pc = discrete_sky(2, 2, test_input)

    assert ele == ele_correct
    assert azi == az_correct
    assert omega == omega_correct
    for i,p in enumerate(pc): assert pytest.approx(p) == expected[i]

@pytest.mark.parametrize("test_input, expected_1, expected_2", [
                                            ("turtle46", hmoy_turtle46, azmoy_turtle46),
                                            (datafile, hmoy_file, azmoy_file), 
                                            ([4, 5, "soc"], hmoy_list, azmoy_list)])
def test_CARIBUsky(test_input, expected_1, expected_2):
    sky = CARIBUsky(test_input)

    # on reuse results from RATP
    sky_correct = ratpformat_to_caribuformat(expected_2,
                                                expected_1,
                                                [0.]*len(expected_2))

    for t, c in zip(sky, sky_correct):
        for x, y in zip(t[1],c[1]) : assert pytest.approx(x, abs=1e-6) == y
