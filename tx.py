from __future__ import print_function

# generate loose packing
from yade import pack, qt, plot

sigmaIso = -1e5
initialFric = 0
finalFric = 30

## define materials
O.materials.append(CohFrictMat(
        density=2650,                       # Density [kg/m3]
        young=1e7,                          # Particle modulus [Pa]
        poisson=.3,                         # Ks/Kn ratio
        frictionAngle=radians(initialFric), # Local friction [rad]
        isCohesive=False,                   # Turn off adhesion
        momentRotationLaw=True,             # Turn on rotational stiffness
        etaRoll=radians(initialFric),       # Rotational friction [rad]
        etaTwist=0,                         # Turn off twisting
        label="granr"                       # Material label
        ))


O.periodic = True
sp = pack.SpherePack()
## uniform distribution
sp.makeCloud((0, 0, 0), (4, 4, 4), rMean=.1, rRelFuzz=.3, periodic=True)

# setup periodic boundary, insert the packing
sp.toSimulation(material="granr")

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom()],
                        [Ip2_FrictMat_FrictMat_FrictPhys(),
                         Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()],
                        [Law2_ScGeom_FrictPhys_CundallStrack(),
                         Law2_ScGeom6D_CohFrictPhys_CohesionMoment()],
                        label='interactionLoop'),
        GlobalStiffnessTimeStepper(
                # use adaptive stiffness-based timestepper
                timestepSafetyCoefficient=0.8,
                timeStepUpdateInterval=100,
                # set to True for domain decomp
                parallelMode=False,
                label='timeStepper'
        ),
        NewtonIntegrator(damping=.2),
        PeriTriaxController(
                label='triax',
                # specify target values and whether they are strains or stresses
                goal=(sigmaIso, sigmaIso, sigmaIso),
                stressMask=7,
                # type of servo-control
                dynCell=True,
                maxStrainRate=(0.01, 0.01, 0.01),
                # wait until the unbalanced force goes below this value
                maxUnbalanced=1e-5,
                relStressTol=1e-5,
                # call this function when goal is reached and the packing is stable
                doneHook='compactionFinished()'
        ),
        PyRunner(command='addPlotData()', iterPeriod=100),
]


def addPlotData():
	plot.addData(
	        unbalanced=unbalancedForce(),
	        i=O.iter,
	        sxx=triax.stress[0],
	        syy=triax.stress[1],
	        szz=triax.stress[2],
	        exx=triax.strain[0],
	        eyy=triax.strain[1],
	        ezz=triax.strain[2],
	        # save all available energy data
	        Etot=O.energy.total(),
	        **O.energy
	)


# enable energy tracking in the code
O.trackEnergy = True

# define what to plot
plot.plots = {
        'i': ('unbalanced',),
        'i ': ('sxx', 'syy', 'szz'),
        ' i': ('exx', 'eyy', 'ezz'),
        # energy plot
        ' i ': (O.energy.keys, None, 'Etot'),
}
# show the plot
plot.plot()


def compactionFinished():
    ## reassign friction values for shear
    # for future contacts change material
	O.materials[0].frictionAngle = radians(finalFric) # radians
	O.materials[0].etaRoll = radians(finalFric)
	# for existing contacts, set contact friction directly
	for i in O.interactions:
	    i.phys.tangensOfFrictionAngle = tan(radians(finalFric))
	    i.phys.maxRollPl = tan(radians(finalFric))

    # set the current cell configuration to be the reference one
	O.cell.trsf = Matrix3.Identity
	# change control type: keep constant confinement in x,y, 20% compression in z
	triax.goal = (sigmaIso, sigmaIso, -.2)
	triax.stressMask = 3
	# allow faster deformation along x,y to better maintain stresses
	triax.maxStrainRate = (1., 1., .01)
	# next time, call triaxFinished instead of compactionFinished
	triax.doneHook = 'triaxFinished()'
	# do not wait for stabilization before calling triaxFinished
	triax.maxUnbalanced = 10


def triaxFinished():
	print('Finished')
	O.pause()

