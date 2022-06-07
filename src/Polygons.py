import math

'''
Classes Vector3 et Triangle3 pour la gestion d'une triangulation personnalisée
'''

class Vector3:
    '''Point 3D
    '''

    def __init__(self, x, y, z):
        '''
        x, y, z : float
        '''
        self.__coord = (x,y,z)
        self.__norm = math.sqrt(x**2+y**2+z**2)

    @property
    def norm(self):
        return self.__norm
    
    def __getitem__(self, id):
        return self.__coord[id]
    
    def __add__(self, v):
        '''addition de vecteurs
        v : Vector3
        return : Vector3
        '''
        return Vector3(self.__coord[0] + v.__coord[0], self.__coord[1] + v.__coord[1], self.__coord[2] + v.__coord[2])

    def __sub__(self, v):
        '''soustraction de vecteurs
        v : Vector3
        return : Vector3
        '''
        return Vector3(self.__coord[0] - v.__coord[0], self.__coord[1] - v.__coord[1], self.__coord[2] - v.__coord[2])

    # note : on peut surcharger une méthode (surcharge de __mul__ avec float et Vector3) avec singledispatchmethod de functools présent dans >= Python 3.8
    def scale(self, a):
        '''multiplication par un scalaire
        a : float
        return : Vector3
        '''
        return Vector3(self.__coord[0] * a, self.__coord[1] * a, self.__coord[2] * a)
    
    def __mul__(self, v):
        '''produit scalaire
        v : Vector3
        return : float
        '''
        return self.__coord[0] * v.__coord[0] + self.__coord[1] * v.__coord[1] + self.__coord[2] + v.__coord[2]
    
    def __xor__(self, v):
        '''produit vectoriel
        v : Vector3
        return : Vector3
        '''
        x = self.__coord[1] * v[2] - v[1] * self.__coord[2]
        y = self.__coord[2] * v[0] - v[2] * self.__coord[0]
        z = self.__coord[0] * v[1] - v[0] * self.__coord[1]
        return Vector3(x,y,z)

    def normalize(self):
        '''normalisation de l'instance
        return : Vector3 (vecteur unitaire de self)
        '''
        return Vector3(self.__coord[0] / self.__norm, self.__coord[1] / self.__norm, self.__coord[2] / self.__norm)
    
    @staticmethod
    def middle(p1, p2):
        '''Calcul le milieu de [p1, p2]
        p1 : Vector3
        p2 : Vector3
        return : Vector3
        '''
        p = p1+p2
        return Vector3(p[0]/2, p[1]/2, p[2]/2)
    
    def __str__(self):
        '''print
        '''
        return "(%.3f, %.3f, %.3f)"%(self.__coord[0], self.__coord[1], self.__coord[2])
    
class Triangle3:
    '''Triangle 3D personnalisé composé de 3 Vector3
    '''
    def __init__(self, v1, v2, v3):
        '''Constructeur
        v1, v2, v3 : Vector3

        calcul le barycentre, l'aire, la normale et l'élévation de la normale
        '''
        self.__vertices = (v1, v2, v3)
        
        self.__barycenter = (v1+v2+v3).scale(1/3)

        side1 = v2-v1
        side2 = v3-v2
        side3 = v1-v3
        v12 = side1^side2
        n = v12+(side2^side3)+(side3^side1)
        self.__normal = n.normalize()
        self.__area = v12.norm * 0.5

        # recupère l'élévation de la normale
        e = math.acos(abs(self.__normal[2]))
        # passe en degré
        e *= 180/math.pi
        # on reste entre 0 et 90
        if e > 90 and e < 180 : e = 90-(e%90)
        else : e = e%90
        
        self.__elevation = e
    
    # Getters/Setters 
    @property
    def barycenter(self):
        return self.__barycenter
    
    @property
    def normal(self):
        return self.__normal

    @property
    def area(self):
        return self.__area

    @property
    def elevation(self):
        return self.__elevation

    def set_id(self, id):
        self.__id = id

    @property
    def id(self):
        return self.__id

    def translate(self, t):
        '''translation de self d'un vecteur t
        t : Vector3
        Màj du barycentre
        '''
        newvert=[]
        for s in self.__vertices:
            newvert.append(s+t)
        self.__vertices = tuple(newvert)

        self.__barycenter = (self.__vertices[0]+self.__vertices[1]+self.__vertices[2]).scale(1/3)
    
    def rescale(self, h):
        '''Mise à l'échelle de self d'un float h
        h : float
        Màj du barycentre et de l'aire
        Màj de l'aire
        '''
        newvert=[]
        for s in self.__vertices:
            newvert.append(s.scale(h))
        self.__vertices = tuple(newvert)
        self.__barycenter = (self.__vertices[0]+self.__vertices[1]+self.__vertices[2]).scale(1/3)
        
        side1 = self.__vertices[1]-self.__vertices[0]
        side2 = self.__vertices[2]-self.__vertices[1]
        self.__area = (side1^side2).norm * 0.5

    def zrotate(self, omega):
        '''Rotation du triangle autour de l'axe z
        omega : angle en degré

        met à jour __vertices et __barycenter
        '''
        omegarad = omega * math.pi/180
        newvert=[]
        for s in self.__vertices:
            newvert.append(Vector3(math.cos(omegarad)*s[0] - math.sin(omegarad)*s[1], 
                                    math.sin(omegarad)*s[0] + math.cos(omegarad)*s[1],
                                    s[2]))
        self.__vertices = tuple(newvert)
        self.__barycenter = (self.__vertices[0]+self.__vertices[1]+self.__vertices[2]).scale(1/3)

    def __getitem__(self, id):
        '''renvoit un sommet
        id : int
        '''
        return self.__vertices[id]

    def __str__(self):
        '''print les info de self
        '''
        return "S1 : "+str(self.__vertices[0])+"\n"+\
                "S2 : "+str(self.__vertices[1])+"\n"+\
                "S3 : "+str(self.__vertices[2])+"\n"+\
                "n : "+str(self.__normal)+"\n"+\
                "Area = "+str(self.__area)+"\n"
