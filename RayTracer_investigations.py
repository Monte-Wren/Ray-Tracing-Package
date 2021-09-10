import RayTracer as rt
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import fmin_tnc

#%% Task 9 lens

T9_lens = rt.SphericalRefraction(20, 0.03, 1, 1.5, 20.000001) #slightly larger than 20 due to floating point errors when bundle radius = aperture radius (some rays dont go through when they should)

#%% function that finds the focus for a given lens, with error

def find_focus(lens, ray_position = 0.1, resolution = 0.01, decimal_places = 3): #can only round values to within resolution
    ray = rt.Ray(np.array([0, ray_position, 0]), np.array([0, 0, 1]))
    lens.propagate_ray(ray)
    op = rt.OutputPlane(140)
    op.propagate_ray(ray, dl = resolution)
    foci = []
    for i in range(0, len(ray.p())):
        if abs(np.around(ray.p()[i, 0], decimal_places)) == 0 and abs(np.around(ray.p()[i, 1], decimal_places)) == 0:
            focus = ray.p()[i, 2] - lens._z0
            foci.append(focus)
        else:
            pass
    
    foci = np.array(foci)
    mean = np.mean(foci)
    error = np.std(foci)
    print('The focus is at', mean, '+/-', error,'/mm')
    return(mean, error)

#agrees with theoretically calculated value

#%% Finds the root mean square deviation from the optical axis for the lens in task 9

bundle1 = rt.Bundle(10, np.array([0,0,1]), radial_seperation = 1)

bundle1.propagate_bundle(T9_lens)
rt.T9_op.propagate_bundle(bundle1)
T9_lens.plot_lens()
bundle1.plot_bundle(loc = 'yz')
plt.show()
bundle1.plot_bundle(loc = 's')
plt.show()
bundle1.plot_bundle(loc = 'op')
plt.show()

print('Root mean square deviation for radius', bundle1._rad, 'mm is:', rt.T9_op.bundle_rms(bundle1), 'mm')

#%%Function that uses scipy.optimize.fmin_tnc to find the optimum lens curvatures to minimies rms deviation of rays from optical axis
#(May take a while)

def find_optimum_lens(focus, z0 = 20, lens_seperation = 5, n1 = 1, n2 = 1.5, aperture_radius = 20, bundle_rad_sep = 1, initial_guess = [0,0]):
    def rms_curvature(curvatures):
        output_plane = rt.OutputPlane(focus)
        bundle1 = rt.Bundle(0.5*aperture_radius, np.array([0,0,1]), radial_seperation = bundle_rad_sep) #makes sure bundle is well within aperture radius as radius of lens could be below aperture radius
        lens1 = rt.SphericalRefraction(z0, curvatures[0], n1, n2, aperture_radius)
        lens2 = rt.SphericalRefraction(z0+lens_seperation, curvatures[1], 1, 1.5, aperture_radius)
        
        bundle1.propagate_bundle(lens1, dl = 0).propagate_bundle(lens2, dl = 0)
        output_plane.propagate_bundle(bundle1, dl = 0)
        
        return output_plane.bundle_rms(bundle1)
    
    opt_curv = opt_curv = fmin_tnc(rms_curvature, x0 = initial_guess, approx_grad = True)
    
    return opt_curv[0]

optimum = find_optimum_lens(100)

print('Optimum lens curvatures are', optimum[0], 'and', optimum[1])

lens1 = rt.SphericalRefraction(20, optimum[0], 1, 1.5, 20)
lens2 = rt.SphericalRefraction(25, optimum[1], 1, 1.5, 20)
OP = rt.OutputPlane(100)

bundle1 = rt.Bundle(10, np.array([0,0,1]), radial_seperation = 1)
bundle1.propagate_bundle(lens1, dl =0).propagate_bundle(lens2, dl = 0)
OP.propagate_bundle(bundle1, dl = 0)
lens1.plot_lens()
lens2.plot_lens()
bundle1.plot_bundle('yz')
plt.show()
bundle1.plot_bundle('op')
plt.show()