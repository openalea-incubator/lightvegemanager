from lightvegemanager.sun import *
import pytest

# global values
day = 250
hour = 10
coordinates = [60.2, 25., 3.6] # Helsinki Finland, UTC+3h
truesolartime = True
localsolartime = False

# for reference:
# https://www.suncalc.org/#/60.1712,24.9327,11/2003.09.07/10:00/1/3
expected_localtime = [0.506493, -0.755064, -0.416345]
expected_globaltime = [0.692525, 0.496989, -0.522886]
@pytest.mark.parametrize("test_input, expected", [
                                                    (truesolartime, expected_globaltime),
                                                    (localsolartime, expected_localtime)])
def test_ratp_sun(test_input, expected):
    result = ratp_sun(day, hour, coordinates, test_input)
    for x,y in zip(result, expected) : assert pytest.approx(x) == y


expected_localtime = [0.506689, -0.754871, -0.416456]
expected_globaltime = [0.692605, 0.496786, -0.522974]
@pytest.mark.parametrize("test_input, expected", [
                                                    (truesolartime, expected_globaltime),
                                                    (localsolartime, expected_localtime)])
def test_caribu_sun(test_input, expected):
    result = caribu_sun(day, hour, coordinates, test_input)
    for x,y in zip(result, expected) : assert pytest.approx(x) == y

def test_print_sun():
    print_sun(day, hour, coordinates, localsolartime)
    assert True
