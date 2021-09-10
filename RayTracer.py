"""
project A:
Ray tracing package
Monte Wren

Follows Cartesian sign convention:
Optical axis is the z-axis.
Must propagate rays from the left of the lens, i.e left to right.
"""
#%% Import packages

import numpy as np
import matplotlib.pyplot as plt

#%% Cell for defining useful functions

def mod(arr):
    """
    Function that returns modulus of numpy array of shape (3,)
    """
    if type(arr) == np.ndarray and np.shape(arr) == (3,):
        modulus = np.sqrt((arr[0]**2)+(arr[1]**2)+(arr[2]**2))
    else:
        raise Exception("Arguments must be numpy.ndarray of shape (3,)")
    
    return modulus

#%% Ray class

class Ray:
    """
    Ray class that initialises at source point p and initial direction vector k
    """
    def __init__(self, p, k): #position p and direction k
        if type(p) == np.ndarray and type(k) == np.ndarray and np.shape(p) == (3,) and np.shape(k) == (3,): #Must be 3D vectors
            self._p = np.array([p])
            self._k = k
        else:
            raise Exception("Arguments must be numpy.ndarray of shape (3,)")
                        
    def __repr__(self):
        return "%s(p=%r, k=%r)" % ("Ray", self._p, self._k)
    
    def __str__(self):
        return "Ray at source {}, direction {}".format(self._p[0], self._k)
            
    def p(self):
        return self._p
    
    def k(self):
        return self._k
    
    def append(self, new_p, new_k):
        """
        method that appends new position pand updates direction vector k
        """
        if type(new_p) == np.ndarray and len(new_p) == 3 and type(new_k) == np.ndarray and len(new_k) == 3:
            self._p = np.vstack((self._p, new_p))
            self._k = new_k

        else:
            raise Exception("Arguments must be numpy.ndarray of shape (3,)")
        
    def vertices(self):
        #I haven't used this method, I did the tasks a slightly different way, 
        #but this would return an array with all the points on the ray path.
        #I've used the p() method instead with indexing
        vertices = []
        for i in range(0, len(self._p)):
            vertices.append(self._p[i])
            return vertices

#%% Optical element base class, this is really just a placeholder
        
class OpticalElement:
     def propagate_ray(self, ray, z_lim, dl = 0.05):
         #In case propagate method is not implemented in inherited class
         raise NotImplementedError()

#%% Output plane class

class OutputPlane(OpticalElement):
    def __init__(self, z):
        self._z = z
    
    def propagate_ray(self, ray, dl = 0): 
        """
        propagates ray from lens to output plane at z
        appends point of intesection with output plane
        OPTIONAL: evaluate each point on the path of the ray with resolution dl
        """
        if ray.p()[-1, 2] > self._z:
            raise Exception('Invalid output plane, must be to the right of lens')
        
        else:
            if ray.k() is None: #Ray is terminated 
                ray._k = None

            else:
                if dl == 0: #This will only calculate and append the point of intersection with the output plane
                    l_refract = ((ray.p()[-1][2]-(self._z))/(-ray.k()[2]))
                    surface_point = ray.p()[-1]
                    new_p = surface_point + l_refract*ray.k()
                    ray.append(new_p, ray.k())
                    
                elif dl < 0:
                    raise Exception('dl must be positive or zero')
                
                else:    
                    l_refract = ((ray.p()[-1][2]-(self._z))/(-ray.k()[2]))
                    surface_point = ray.p()[-1]
                    for i in range(1, int(l_refract/dl)+1):
                        new_p = surface_point + i*dl*ray.k()
                        ray.append(new_p, ray.k())
            
    def propagate_bundle(self, bundle, dl = 0):
        """
        loops through and propagates entire bundle to output plane at z
        So I dont have to write out endless for loops
        """
        for i in range(0, len(bundle._bundle)):
            self.propagate_ray(bundle._bundle[i], dl = dl)

    def bundle_rms(self, bundle):
        """
        calculates rms deviation of ray points from optical axis at the output plane
        Ray bundle must be parallel and centred on optical axis
        """
        sum_residuals_sq = 0
        for i in range(0, len(bundle.get_rays())):
            if bundle._bundle[i].k() is None: #Dont want to include rays that dont make it through lens
                pass
            else:
                #loops through bundle and calculates square of distance of rays from the optical axis
                residual_sq = (bundle.get_rays()[i].p()[-1, 0])**2 + (bundle.get_rays()[i].p()[-1, 1])**2
                sum_residuals_sq = sum_residuals_sq + residual_sq
                #then calculates and returns rmsd
            rmsd = np.sqrt(sum_residuals_sq/len(bundle.get_rays()))
        return rmsd
                    

#%% Spherical refraction class, this cell includes everything to do with the spherical lenses
        
class SphericalRefraction(OpticalElement):
     def __init__(self, z0, curvature, n1, n2, ap_radius):
        """
        z0 - the intercept of the surface with the z-axis, origin of lens is on z-axis
        curvature - the curvature of the surface
        n1, n2 - the refractive indices either side of the surface
        ap_radius - the maximum extent of the surface from the optical axis
        """
        self._z0 = z0
        self._curvature = curvature
        self._n1 = n1
        self._n2 = n2
        self._ap_radius = ap_radius

     def __repr__(self):
         return "%s(z0=%g, curvature=%g, n1=%g, n2=%g, ap_radius=%g)" % ("Spherical lens", self._z0, self._curvature, self._n1, self._n2, self._ap_radius)

     def __str__(self):
         return "Spherical lens, z0 = %g, curvature = %g, n1 = %g, n2 = %g and aperture radius = %g" % (self._z0, self._curvature, self._n1, self._n2, self._ap_radius)
     
     def ref_indices(self):
         return self._n1, self._n2
     
     def plot_lens(self, dx = 0.01):
         """
         Option to plot lens, will be either a plane or arc in xz/yz-plane
         """
         if self._curvature == 0:
             #plots a plane at z0
             y = np.arange(-self._ap_radius, self._ap_radius, dx)
             x = np.full(shape = len(y), fill_value = self._z0)
             plt.plot(x, y, 'k')
        
         else:
             #plots a circle in plane that intersects the optical axis at z0
             x = np.arange(self._z0, self._z0+1/self._curvature, dx)
             y_plus = []
             y_minus = []
             for i in range(0, len(x)):
                 if abs(np.sqrt((1/self._curvature)**2-(x[i]-(self._z0+(1/self._curvature)))**2)) > self._ap_radius:
                     #If y is larger than aperture radius then return None
                     y_plus.append(None)
                     y_minus.append(None)
                 else:
                     y_plus.append(abs(np.sqrt((1/self._curvature)**2-(x[i]-(self._z0+1/self._curvature))**2)))
                     y_minus.append(-1*abs(np.sqrt((1/self._curvature)**2-(x[i]-(self._z0+1/self._curvature))**2)))
             plt.plot(x, y_plus, 'k')
             plt.plot(x, y_minus, 'k')
             
             y = np.arange(-self._ap_radius, self._ap_radius, dx)
             x = np.full(shape = len(y), fill_value = self._z0)
             plt.plot(x, y, 'k')

     def intercept(self, ray):
        """
        Finds the length l of the source of the ray and the point of intersection with the lens
        """
        if ray.p()[-1, 2] > self._z0: #Ray must start to the left of z0 
             return None
        else:
            unit_k = (1/mod(ray.k()))*ray.k()
            if self._curvature == 0: #special case where curvature = 0
                if ray.k()[2] == 0: #This checks if ray will ever hit lens
                    return None
                else:
                    l = (self._z0 - ray.p()[-1][2])/unit_k[2]
                    if ray.p()[-1, 0]**2+ray.p()[-1, 1]**2 > self._ap_radius**2: #Checks if ray misses aperture
                        return None
                    else:
                        return l
            else:
                r = ray.p()[-1] - np.array([0, 0 , self._z0 + (1/self._curvature)]) #vector from centre of curvature to ray source point
                if np.dot(r, unit_k)**2 < (mod(r)**2)-(1/abs(self._curvature))**2:
                    #no solution condition, does not intersect with sphere 
                    return None
                
                else:
                    r = ray.p()[-1] - np.array([0, 0 , self._z0 + (1/self._curvature)])
                    l_1 = -np.dot(r, unit_k) + np.sqrt(((np.dot(r, unit_k))**2)-((mod(r)**2)-(1/abs(self._curvature))**2)) #Length from point to point of intersection
                    l_2 = -np.dot(r, unit_k) - np.sqrt(((np.dot(r, unit_k))**2)-((mod(r)**2)-(1/abs(self._curvature))**2))
                    p1 = ray.p()[-1] + l_1*unit_k #Points of intersection
                    p2 = ray.p()[-1] + l_2*unit_k
                    z0_vect = np.array([0, 0, self._z0])
                    dist1 = mod(p1 - z0_vect) #distances between points of intersection and z0
                    dist2 = mod(p2 - z0_vect)
                    if dist1 < dist2:
                        if (((ray.p()[-1]+l_1*unit_k)[0])**2)+(((ray.p()[-1]+l_1*unit_k)[1])**2) > self._ap_radius**2:
                            #misses aperture
                            return None
                        else:
                            return l_1
                            
                    elif dist1 > dist2:
                        if (((ray.p()[-1]+l_2*unit_k)[0])**2)+(((ray.p()[-1]+l_2*unit_k)[1])**2) > self._ap_radius**2:
                            return None
                        else:
                            return l_2
                        
                    elif dist1 == dist2:
                        return l_1
                    
                    else:
                        pass
        
     def refract(self, ray):
         #returns refracted direction unit vector unit_k2 for incident ray on spherical lens
         unit_k1 = (1/mod(ray.k()))*ray.k() #incident unit direction vector
         if self.intercept(ray) is None: #Does ray have a valid intercept?
             return None #I did have it return unit_k1 as it would carry on in the same direction, but lets assume it cannot propagate past lens
         
         else:
             l = self.intercept(ray)
             n1_2 = (self._n1/self._n2)
             if self._curvature == 0: #special case for zero curvature
                 unit_n = np.array([0, 0, -1])
                 if mod(np.cross(unit_n, unit_k1)) > n1_2: #TIR condition (sin(theta) > n1_2)
                     return None
                 else:
                     c = -np.dot(unit_n, unit_k1)
                     k_refracted = n1_2*unit_k1 + (n1_2*c - np.sqrt(1 - ((n1_2**2)*(1-c**2))))*unit_n #snells law with no trig
                     unit_k_refracted = (1/mod(k_refracted))*k_refracted
                     return unit_k_refracted #returns refracted unit direction vector
             
             elif self._curvature > 0:
                 r = ray.p()[-1] - np.array([0, 0 , self._z0 + (1/self._curvature)])
                 unit_n = 1/(1/abs(self._curvature))*(r + l*unit_k1) #uses geometry of problem, normal vector must point towards the incoming ray
                 if mod(np.cross(unit_n, unit_k1)) > n1_2:
                     return None
                 else:
                     c = -np.dot(unit_n, unit_k1)
                     k_refracted = n1_2*unit_k1 + ((n1_2*c) - np.sqrt(1-((n1_2**2)*(1-c**2))))*unit_n
                     unit_k_refracted = (1/mod(k_refracted))*k_refracted
                     return unit_k_refracted
             
             else: #exactly the same process as for positive curvature except negative unit vector as normal vector must point towards incoming ray
                 r = ray.p()[-1] - np.array([0, 0 , self._z0 + (1/self._curvature)])
                 unit_n = -1/(1/abs(self._curvature))*(r + l*unit_k1)
                 if mod(np.cross(unit_n, unit_k1)) > n1_2:
                     return None
                 else:
                     c = -np.dot(unit_n, unit_k1)
                     k_refracted = n1_2*unit_k1 + ((n1_2*c) - np.sqrt(1 - ((n1_2**2)*(1-c**2))))*unit_n
                     unit_k_refracted = (1/mod(k_refracted))*k_refracted
                     return unit_k_refracted


     def propagate_ray(self, ray, dl = 0): 
         #Evaluates each point on the path of the ray with resolution dl, or just calculate points of intersection with lenses
         unit_k1 = (1/mod(ray.k()))*ray.k() #incident unit vector
         unit_k2 = self.refract(ray) #refracted unit vector
         if unit_k2 is None:
             ray._k = None #Ray is terminated as it is of no observational interest
         
         else:
             if dl == 0: #This will only calculate and append the point of intersection with the lens
                 l = self.intercept(ray)
                 surface_point = ray.p()[-1]
                 new_p = surface_point + l*unit_k1
                 ray.append(new_p, unit_k2)
             
             elif dl < 0:
                 raise Exception('dl must be positive or zero')
                 
             else:
                l = self.intercept(ray)
                surface_point = ray.p()[-1]
                for i in range(1, int(l/dl)+1):
                     new_p = surface_point + i*dl*unit_k1
                     ray.append(new_p, unit_k2)
             return ray
#%% Ray bundle class 

class Bundle(OpticalElement):
    #Cylindrical bundle of rays of uniform flux density
    def __init__(self, radius, initial_direction, radial_seperation = 1, angular_seperation = (2*np.pi)/6, initial_position = np.array([0, 0, 0])):
        self._rad = radius
        self._rs = radial_seperation #radial distance between each point
        self._as = angular_seperation #angle turned before point is placed
        self._k = initial_direction 
        self._initial = initial_position #position of central ray in bundle
        
        bundle = [Ray(self._initial, self._k)] #collimated so all rays have direction k
        rad_seperation = np.arange(0, self._rad, self._rs)
        for i in rad_seperation:
            angle_seperation = np.arange(0, 2*np.pi, self._as/(int(i/self._rs)+1))
            #angular seperation decreases by a factor of i as the radius increases so that the bundle is uniform
            for j in angle_seperation: #turn through given angle then increase the radius by given radial separation until it reaches given radius
                x = self._initial[0]+((i+self._rs)*np.cos(j))
                y = self._initial[1]+((i+self._rs)*np.sin(j))
                z = self._initial[2]
                #ray source points arranged in a circle about initial position with increasing radius
                ray = Ray(np.array([x, y, z]), self._k) 
                bundle.append(ray)
        self._bundle = bundle
    
    def __repr__(self):
        return "%s(radius=%r, direction=%r, radial_seperation = %r, angular_seperation = %r, initial_position = %r" % ("Bundle", self._rad, self._rs, self._as, self._k, self._initial)
    
    def __str__(self):
        return "Bundle of rays centred at {}, radius {} and direction {}".format(self._initial, self._rad, self._k)     
    
    def get_rays(self):
        return self._bundle
            
    def propagate_bundle(self, lens, dl = 0):
        #exactly the same as the propagate method in sphericalrefraction class
        for i in range(0, len(self._bundle)):#Loops through and propagates each ray in bundle
            unit_k1 = (1/mod(self._bundle[i].k()))*self._bundle[i].k()
            unit_k2 = lens.refract(self._bundle[i])
            if unit_k2 is None:
                self._bundle[i]._k = None
        
            else:
                if dl == 0:
                    l = lens.intercept(self._bundle[i])
                    surface_point = self._bundle[i].p()[-1]
                    new_p = surface_point + l*unit_k1
                    self._bundle[i].append(new_p, unit_k2)
                    
                elif dl < 0:
                    raise Exception('dl must be positive or zero')
                    
                else:
                    l = lens.intercept(self._bundle[i])
                    surface_point = self._bundle[i].p()[-1]
                    for j in range(1, int(l/dl)+1):
                        new_p = surface_point + j*dl*unit_k1
                        self._bundle[i].append(new_p, unit_k2)
        return self
    
    def plot_bundle(self, loc = 'op', col = 'r'):
        #Plots a ray bundle at a given location
        #'s' = Source spot diagram
        #'op' = Output plane spot diagram
        #'xz' = Plots ray diagram in xz plane 
        #'yz' = Plots ray diagram in yz plane
        if loc == 's':
            for i in range(0, len(self._bundle)):
                plt.title('Source')
                plt.axis('equal')
                plt.xlabel('/mm')
                plt.ylabel('/mm')
                plt.plot(self._bundle[i].p()[0, 0], self._bundle[i].p()[0, 1], '%s.' %col)
        
        elif loc == 'op':
            for i in range(0, len(self._bundle)):
                if self._bundle[i].k() is None: #Dont want to plot rays that dont go through lens
                    pass
                else:
                    plt.title('Output Plane')
                    plt.axis('equal')
                    plt.xlabel('/mm')
                    plt.ylabel('/mm')
                    plt.plot(self._bundle[i].p()[-1, 0], self._bundle[i].p()[-1, 1], '%s.' %col)
        
        elif loc == 'xz':
            for i in range(0, len(self._bundle)):
                plt.title('xz - plane')
                plt.xlabel('/mm')
                plt.ylabel('/mm')
                plt.plot(self._bundle[i].p()[:, 2], self._bundle[i].p()[:, 0], col)
            
        elif loc == 'yz':
            for i in range(0, len(self._bundle)):
                plt.title('yz - plane')
                plt.xlabel('/mm')
                plt.ylabel('/mm')
                plt.plot(self._bundle[i].p()[:, 2], self._bundle[i].p()[:, 1], col)

        else:
            raise Exception("Give a valid location: 's', 'op', 'xz' or 'yz'")
            
#%% Useful lenses and that

T9_lens = SphericalRefraction(20, 0.03, 1, 1.5, 2)
T9_op = OutputPlane(120)

triangle_lens1 = SphericalRefraction(20, 0.03, 1, 1.5, 0.5)
triangle_lens2 = SphericalRefraction(220, -0.03, 1.5, 1, 0.5)

T12_lens = SphericalRefraction(20, 0.03, 1, 1.5, 10)
T12_op = OutputPlane(120)

T15_lens1 = SphericalRefraction(20, 0, 1, 1.5168, 7.5)
T15_lens2 = SphericalRefraction(25, 0.02, 1, 1.5168, 7.5)
T15_lens3 = SphericalRefraction(20, 0.02, 1, 1.5168, 7.5)
T15_lens4 = SphericalRefraction(25, 0, 1, 1.5168, 7.5)

T15_op1 = OutputPlane(171.75)
T15_op2 = OutputPlane(240)
