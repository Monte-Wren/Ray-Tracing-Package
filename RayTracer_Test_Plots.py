"""
Here is the code the generates all the test plots
"""


#%%

import RayTracer as rt #import package
import numpy as np
import matplotlib.pyplot as plt

#%% Task-9 plots

T9_op = rt.OutputPlane(120)

T9_ray1 = rt.Ray(np.array([0,1,0]), np.array([0,0,1]))
T9_ray2 = rt.Ray(np.array([0,0,0]), np.array([0,0,1]))
T9_ray3 = rt.Ray(np.array([0,0,0]), np.array([0,0.01,1]))
rt.T9_lens.propagate_ray(T9_ray1)
rt.T9_lens.propagate_ray(T9_ray2)
rt.T9_lens.propagate_ray(T9_ray3)
rt.T9_op.propagate_ray(T9_ray1)
rt.T9_op.propagate_ray(T9_ray2)
rt.T9_op.propagate_ray(T9_ray3)

plt.plot(T9_ray1.p()[:, 2], T9_ray1.p()[:, 1])
plt.plot(T9_ray2.p()[:, 2], T9_ray2.p()[:, 1])
plt.plot(T9_ray3.p()[:, 2], T9_ray3.p()[:, 1])
plt.legend(['Ray 1', 'Ray 2', 'Ray 3'])
plt.ylabel('/mm')
plt.xlabel('/mm')
plt.title('Few example rays through Task 9 lens')
rt.T9_lens.plot_lens()
plt.show()

#%% Task 11 misc plots

#Image of triangle sent through lens

T11_ray1 = rt.Ray(np.array([0.1, 0, 0]), np.array([0, 0, 1]))
T11_ray2 = rt.Ray(np.array([-1/20, np.sqrt(3)/20, 0]), np.array([0, 0, 1]))
T11_ray3 = rt.Ray(np.array([-1/20, -np.sqrt(3)/20, 0]), np.array([0, 0, 1]))

triangle_op = rt.OutputPlane(250)

rt.triangle_lens1.propagate_ray(T11_ray1)
rt.triangle_lens1.propagate_ray(T11_ray2)
rt.triangle_lens1.propagate_ray(T11_ray3)
rt.triangle_lens2.propagate_ray(T11_ray1)
rt.triangle_lens2.propagate_ray(T11_ray2)
rt.triangle_lens2.propagate_ray(T11_ray3)

triangle_op.propagate_ray(T11_ray1)
triangle_op.propagate_ray(T11_ray2)
triangle_op.propagate_ray(T11_ray3)

plt.plot(T11_ray1.p()[:, 2], T11_ray1.p()[:, 1])
plt.plot(T11_ray2.p()[:, 2], T11_ray2.p()[:, 1])
plt.plot(T11_ray3.p()[:, 2], T11_ray3.p()[:, 1])
plt.ylabel('/mm')
plt.xlabel('/mm')
plt.title('Triangle image')
plt.legend(['Ray 1', 'Ray 2', 'Ray 3'])
rt.triangle_lens1.plot_lens(dx = 0.001)
rt.triangle_lens2.plot_lens(dx = 0.001)
plt.show()

x = [T11_ray1.p()[0, 0], T11_ray2.p()[0, 0], T11_ray3.p()[0, 0], T11_ray1.p()[-1, 0], T11_ray2.p()[-1, 0], T11_ray3.p()[-1, 0]]
y = [T11_ray1.p()[0, 1], T11_ray2.p()[0, 1], T11_ray3.p()[0, 1], T11_ray1.p()[-1, 1], T11_ray2.p()[-1, 1], T11_ray3.p()[-1, 1]]

for i in range(0, 2):
    plt.plot([x[i],x[i+1]], [y[i], y[i+1]], 'bo-')

for i in range(3, 5):
    plt.plot([x[i],x[i+1]], [y[i], y[i+1]], 'ro-')

plt.plot([x[0], x[2]], [y[0], y[2]], 'bo-', label = 'Incident triangle')
plt.plot([x[3], x[5]], [y[3], y[5]], 'ro-', label = 'Output plane')

plt.legend(loc = 1)
plt.ylabel('/mm')
plt.xlabel('/mm')
plt.title('Lens reverses images')

plt.show()

#%% Task 12 plot

T12_bundle1 = rt.Bundle(5, np.array([0 ,0, 1]))
T12_bundle2 = rt.Bundle(5, np.array([0, 0.1, 1]), initial_position = np.array([0 , -2, 0]))

T12_bundle1.propagate_bundle(rt.T12_lens)
T12_bundle2.propagate_bundle(rt.T12_lens)
rt.T12_op.propagate_bundle(T12_bundle1)
rt.T12_op.propagate_bundle(T12_bundle2)

T12_bundle1.plot_bundle(loc = 'yz', col = 'b')
T12_bundle2.plot_bundle(loc = 'yz', col = 'r')
plt.title('Task 12: Two ray bundles through lens from task 9')
plt.show()
T12_bundle1.plot_bundle(loc = 's', col = 'b')
plt.title('Source (Task 12)')
plt.show()
T12_bundle2.plot_bundle(loc = 's', col = 'r')
plt.title('Source (Task 12)')
plt.show()
T12_bundle1.plot_bundle(loc = 'op', col = 'b')
plt.title('Output plane (Task 12) - Spherical aberration can be clearly seen')
plt.show()
T12_bundle2.plot_bundle(loc = 'op', col = 'r')
plt.title('Output Plane - Comatic aberration can be clearly seen')
plt.show()

#%% Plano-convex lens

#Planar then convex lens
T15_bundle1 = rt.Bundle(5, np.array([0 ,0, 1]))
T15_bundle2 = rt.Bundle(5, np.array([0, 0.1, 1]), initial_position = np.array([0 , -2, 0]))

T15_bundle1.propagate_bundle(rt.T15_lens1).propagate_bundle(rt.T15_lens2)
T15_bundle2.propagate_bundle(rt.T15_lens1).propagate_bundle(rt.T15_lens2)
rt.T15_op1.propagate_bundle(T15_bundle1)
rt.T15_op1.propagate_bundle(T15_bundle2)

T15_bundle1.plot_bundle(loc = 'yz', col = 'b')
T15_bundle2.plot_bundle(loc = 'yz', col = 'r')
rt.T15_lens1.plot_lens()
rt.T15_lens2.plot_lens()
plt.title('Task 15: Two ray bundles through plano-convex lens (plane then convex)')
plt.show()
T15_bundle1.plot_bundle(loc = 'op', col = 'b')
plt.title('Source (Task 15): Plane then convex')
plt.show()
T15_bundle2.plot_bundle(loc = 'op', col = 'r')
plt.title('Output plane (Task 15): Plane then convex')
plt.show()


#convex then planar lens
T15_bundle1 = rt.Bundle(5, np.array([0 ,0, 1]))
T15_bundle2 = rt.Bundle(5, np.array([0, 0.1, 1]), initial_position = np.array([0 , -2, 0]))

T15_bundle1.propagate_bundle(rt.T15_lens3, dl = 0).propagate_bundle(rt.T15_lens4)
T15_bundle2.propagate_bundle(rt.T15_lens3, dl = 0).propagate_bundle(rt.T15_lens4)
rt.T15_op2.propagate_bundle(T15_bundle1)
rt.T15_op2.propagate_bundle(T15_bundle2)

T15_bundle1.plot_bundle(loc = 'yz', col = 'b')
T15_bundle2.plot_bundle(loc = 'yz', col = 'r')
rt.T15_lens3.plot_lens()
rt.T15_lens4.plot_lens()
plt.title('Task 15: Two ray bundles through plano-convex lens (convex then plane)')
plt.show()


