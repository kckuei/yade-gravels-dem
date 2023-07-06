# periodic triaxial test simulation
#
# The initial packing is either
#
# 1. random cloud with uniform distribution, or
# 2. cloud with specified granulometry (radii and percentages), or
# 3. cloud of clumps, i.e. rigid aggregates of several particles
#
# The triaxial consists of 2 stages:
#
# 1. isotropic compaction, until sigmaIso is reached in all directions;
#    this stage is ended by calling compactionFinished()
# 2. constant-strain deformation along the z-axis, while maintaining
#    constant stress (sigmaIso) laterally; this stage is ended by calling
#    triaxFinished()
#
# Controlling of strain and stresses is performed via PeriTriaxController,
# of which parameters determine type of control and also stability
# condition (maxUnbalanced) so that the packing is considered stabilized
# and the stage is done.
#
from __future__ import print_function


## Inputs
sigmaIso = -1e5         # Target stress [Pa]
initialFric = 0.0       # Initial friction at specimen genesis
finalFric = 0.5         # Final friction before shearing
strain_inc = 0.01 #0.0001     # Strain increment to take before every force rebalance
targetStrain = 0.4      # Target strain (+val comp), must be a multiple of the strain incr



#import matplotlib
#matplotlib.use('Agg')

# generate loose packing
from yade import pack, qt, plot

O.periodic = True
sp = pack.SpherePack()
if 1:
	## uniform distribution
	sp.makeCloud((0, 0, 0), (4, 4, 4), rMean=.1, rRelFuzz=.3, periodic=True)
else:
	## create packing from clumps
	# configuration of one clump
	c1 = pack.SpherePack([((0, 0, 0), .03333), ((.03, 0, 0), .017), ((0, .03, 0), .017)])
	# make cloud using the configuration c1 (there could c2, c3, ...; selection between them would be random)
	sp.makeClumpCloud((0, 0, 0), (16, 16, 16), [c1], periodic=True, num=2000)

# setup periodic boundary, insert the packing
sp.toSimulation()

# assign initial friction
O.materials[0].frictionAngle = initialFric

# define strain target levels and initialize target index
#strain_levels = numpy.arange(strain_inc, targetStrain + 10*strain_inc, strain_inc)

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()],
                              label='collider'),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom()], 
                        [Ip2_FrictMat_FrictMat_FrictPhys()], 
                        [Law2_ScGeom_FrictPhys_CundallStrack()],
                        label='interactionLoop'),
        GlobalStiffnessTimeStepper(
                # use adaptive stiffness-based timestepper
                timestepSafetyCoefficient=0.8,
                timeStepUpdateInterval=100,
                # set to True for domain decomp
                parallelMode=False,
                label='timeStepper'
        ),
        NewtonIntegrator(
                # non-viscous newton damping
                damping=.2,
                # create label 
                label='newton'
        ),
        PeriTriaxController(
                # create label/var name for triax engine
                label='triax',
                # specify target values and whether they are strains or stresses
                goal=(sigmaIso, sigmaIso, sigmaIso),
                stressMask=7,
                # type of servo-control
                dynCell=True,
                maxStrainRate=(.01, .01, .01),
                # wait until the unbalanced force goes below this value
                maxUnbalanced=1e-5,
                relStressTol=1e-5,
                # call this function when goal is reached and the packing is stable
                doneHook='compactionFinished()'
        ),
        PyRunner(command='addPlotData()', iterPeriod=100),
]
#O.dt = .5 * PWaveTimeStep()


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
            poros=yade._utils.porosity(),
            porosvox=yade._utils.voxelPorosity(200,*yade._utils.aabbExtrema()),
            # mechanical coordination number Zm
            Zm=utils.avgNumInteractions(skipFree=True),
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
        '  i': ('Zm', None, 'poros', 'porosvox'),
        'ezz': ('sxx','syy','szz')
}
# show the plot
plot.plot()



def compactionFinished():
    ## set friction angle back to non-zero value and prepare for shear
	# tangensOfFrictionAngle is computed by the Ip2_* functor from material
	# for future contacts change material (there is only one material for all particles)
	O.materials[0].frictionAngle = finalFric  # radians
	# for existing contacts, set contact friction directly
	for i in O.interactions:
	    i.phys.tangensOfFrictionAngle = tan(finalFric)

    # set the current cell configuration to be the reference one
	O.cell.trsf = Matrix3.Identity
    # next time, start shearing specimen
	triax.doneHook = 'startShear()'
	# init strain targets, and save strain inc and final target
	triax.strain_inc = strain_inc
	triax.strain_target_stage = strain_inc
	triax.strain_target_final = targetStrain
	# set grow damping
	triax.growDamping = 0.02
	

def startShear():
	# change control type: keep constant confinement in x,y, compression in z
	triax.goal = (sigmaIso, sigmaIso, -triax.strain_target_stage)
	# print stage info
	print("Target strain stage: ", triax.goal)
	triax.stressMask = 3
	# allow faster deformation along x,y to better maintain stresses
	triax.maxStrainRate = (1., 1., .01)
	# set the done hook based on stage of test
	if (triax.strain[2] <= -triax.strain_target_final):
	    # reached target strain, next time, call triaxFinished instead of startShear
	    triax.doneHook = 'triaxFinished()'
	else:
	    # otherwise update strain target level and recursively call startShear
	    triax.strain_target_stage += triax.strain_inc
	    triax.doneHook = 'startShear()'
	 
	# wait for stabilization before calling calling next doneHook
	triax.maxUnbalanced= 10 #0.001
    


def triaxFinished():
	print('Finished')
	O.pause()

