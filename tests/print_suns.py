from src.LightVegeManager_sun import print_sun

if __name__ == "__main__":
    print_sun(201, 12, [0., 0., 0.] , True)
    print_sun(201, 12, [40., 0., 0.] , True)
    print_sun(201, 16, [40., 12., 2.] , True)
    print_sun(201, 16, [40., 12., 2.] , False)
    print_sun(347, 16, [46., 0., 0.] , True)
