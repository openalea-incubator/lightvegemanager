from LightVegeManager_sun import print_sun

"""
    A call of print_sun will create sun positions depending on time and location. 
    Sun position is computed using CARIBU and RATP intern algorithms
    Then, the function prints results for each model.

    If this test works, it means all the dependencies are correctly installed
"""

if __name__ == "__main__":
    # definition
    # print_sun(day, hour, [longitude, latitude, timezone], truesolartime)

    print_sun(201, 12, [0., 0., 0.] , True)
    print_sun(201, 12, [40., 0., 0.] , True)
    print_sun(201, 16, [40., 12., 2.] , True)
    print_sun(201, 16, [40., 12., 2.] , False)
    print_sun(347, 16, [46., 0., 0.] , True)
