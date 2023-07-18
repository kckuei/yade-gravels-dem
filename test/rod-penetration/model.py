# -*- encoding=utf-8 -*-

from __future__ import print_function
from builtins import range
import random
from yade import ymport

## PhysicalParameters

## Variant of mesh
mesh = 'coarse'
#mesh = 'fine'
#mesh = 'tiny'

## Import geometry
rod = O.bodies.append(ymport.stl('rod-' + mesh + '.stl', wire=True))

# Spheres
sphereRadius = 0.01
nbSpheres = (32, 11, 32)
print("Creating %d spheres..." % (nbSpheres[0] * nbSpheres[1] * nbSpheres[2]), end=' ')
for i in range(nbSpheres[0]):
	for j in range(nbSpheres[1]):
		for k in range(nbSpheres[2]):
			x = (i * 2 - nbSpheres[0]) * sphereRadius * 1.1 + sphereRadius * random.uniform(-0.1, 0.1)
			y = -j * sphereRadius * 2.2 - 0.01
			z = (k * 2 - nbSpheres[2]) * sphereRadius * 1.1 + sphereRadius * random.uniform(-0.1, 0.1)
			r = random.uniform(sphereRadius, sphereRadius * 0.9)
			fixed = False
			color = [0.51, 0.52, 0.4]
			if (i == 0 or i == nbSpheres[0] - 1 or j == nbSpheres[1] - 1 or k == 0 or k == nbSpheres[2] - 1):
				fixed = True
				color = [0.21, 0.22, 0.1]
			O.bodies.append(sphere([x, y, z], r, color=color, fixed=fixed))
print("done\n")

## Estimate time step
#O.dt=PWaveTimeStep()
O.dt = 0.0001

## Engines
O.engines = [
        ## Resets forces and momenta the act on bodies
        ForceResetter(),

        ## Using bounding boxes find possible body collisions.
        InsertionSortCollider([
                Bo1_Sphere_Aabb(),
                Bo1_Facet_Aabb(),
        ]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        ## Apply gravity
        ## Motion equation
        NewtonIntegrator(damping=0.3, gravity=[0, -9.81, 0]),
        ## Apply kinematics to rod
        TranslationEngine(ids=rod, translationAxis=[0, -1, 0], velocity=0.1),
        ## Save force on rod
        #ForceRecorder(ids=rod,file='force-'+mesh+'.dat',iterPeriod=50),
]

import sys, time

print("Start simulation: " + mesh)
nbIter = 10000

from yade import qt
qt.View()

O.stopAtIter = nbIter
O.run()

#for t in xrange(2):
#	start=time.time();O.run(nbIter);O.wait();finish=time.time()
#	speed=nbIter/(finish-start); print '%g iter/sec\n'%speed
#print "FINISH"
#quit()
