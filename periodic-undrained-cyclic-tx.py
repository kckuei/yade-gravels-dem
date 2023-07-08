# periodic cyclic triaxial test simulation
#
# The initial packing is either
#
# 1. random cloud with uniform distribution, or
# 2. cloud with specified granulometry (radii and percentages), or
#
# The triaxial consists of 2 stages:
#
# 1. isotropic compaction, until sigmaIso is reached in all directions;
#    this stage is ended by calling compactionFinished()
# 2. constant-strain deformation along the xyz axes, while maintaining
#    constant volume; this stage is ended by calling triaxFinished()
#
# Controlling of strain and stresses is performed via PeriTriaxController,
# of which parameters determine type of control and also stability
# condition (maxUnbalanced) so that the packing is considered stabilized
# and the stage is done.


## Imports
from __future__ import print_function
from yade import pack, qt, plot

#import matplotlib
#matplotlib.use('Agg')


## Inputs
sigmaIso = -1e5             # Target stress [Pa]
initialFric = 30            # Initial friction at specimen genesis [deg]
finalFric = 30              # Final friction before shearing [deg]
strain_inc = 0.0001 #0.01   # Strain increment to take before every force rebalance
stageMaxUnbalanced = 10.    # Max unbalanced force for each strain inc stage
targetStrain = 0.15         # Target strain, e.g. 3%, 5%
targetCSR = 0.2             # Target cyclic stress ratio, CSR_txc = qd/(2*sigma'co)
triaxGrowDamping = 0.02     # Triax grow damping

pk_maxStrainRate = 0.01     # Max isotropic strain rate
pk_maxUnbalanced = 1e-5     # Max unbalanced force for iso pack
pk_relStressTol = 1e-5      # Stress tolerance for iso pack

nmin = (0, 0, 0)            # Min cell boundary
nmax = (4, 4, 4)            # Max Cell boundary

## Define material
O.materials.append(CohFrictMat(
        density=2650,               # Density [kg/m3]
        young=1e7,                  # Particle modulus [Pa]
        poisson=.3,                 # Ks/Kn ratio
        frictionAngle=radians(30),  # Local friction [rad]
        isCohesive=False,           # Turn off adhesion
        momentRotationLaw=True,     # Turn on rotational stiffness
        etaRoll=radians(30),        # Rotational friction [rad]
        etaTwist=0,                 # Turn off twisting
        label="granr"               # Material label
        ))

## Setup periodic boundary and generate loose packing
# Periodic boundary
O.periodic = True
# Define a spherical packing
sp = pack.SpherePack()
# Uniform distributiomn
sp.makeCloud(nmin, nmax, rMean=.1, rRelFuzz=.3, periodic=True)
# Insert the packing
sp.toSimulation(material="granr")

# Overwrite with initial friction for compaction
O.materials[0].frictionAngle = radians(initialFric)
O.materials[0].etaRoll = radians(initialFric)

## Define engines
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()],
                              label='collider'),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom()], 
                        [Ip2_FrictMat_FrictMat_FrictPhys(),
                         Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()
                        ], 
                        [Law2_ScGeom_FrictPhys_CundallStrack(),
                         Law2_ScGeom6D_CohFrictPhys_CohesionMoment()
                        ],
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
                maxStrainRate=(pk_maxStrainRate, pk_maxStrainRate, pk_maxStrainRate),
                # wait until the unbalanced force goes below this value
                maxUnbalanced=pk_maxUnbalanced,
                relStressTol=pk_relStressTol,
                # call this function when goal is reached and the packing is stable
                doneHook='compactionFinished()'
        ),
        PyRunner(command='addPlotData()', iterPeriod=100),
]
#O.dt = .5 * PWaveTimeStep()


def mean_stress():
    pm=triax.stressTensor.trace()/3
    return pm

def deviatoric_stress():
    qd=triax.stress[2]-((triax.stress[0]+triax.stress[1])/2)
    return qd

def porosity(factor=1.0, resolution=200):
    vmin, vmax = yade._utils.aabbExtrema()
    return yade._utils.voxelPorosity(resolution, vmin*factor, vmax*factor)

def addPlotData():
	plot.addData(
            i=O.iter, 
            # Unbalanced force
	        unbalanced=unbalancedForce(),
            # stress strain
	        sxx=triax.stress[0],
	        syy=triax.stress[1],
	        szz=triax.stress[2],
	        exx=triax.strain[0],
	        eyy=triax.strain[1],
	        ezz=triax.strain[2],
	        # mean and deviatoric stress
	        pm=mean_stress(),
	        qd=deviatoric_stress(),
	        # porosity
	        poros=yade._utils.porosity(),
	        porosvox0=yade._utils.voxelPorosity(200,*yade._utils.aabbExtrema()),
	        porosvox1=porosity(factor=0.5),
	        porosvox2=porosity(factor=0.9),
	        porosvox3=porosity(factor=1.1),
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
        # unbalanced force plot
        'i': ('unbalanced',),
        # stress plot
        'i ': ('sxx', 'syy', 'szz'),
        # strain plot
        ' i': ('exx', 'eyy', 'ezz'),
        # energy plot
        ' i ': (O.energy.keys, None, 'Etot'),
        # coordination porosity plot
        '  i': ('Zm', None, 'poros', 'porosvox0', 'porosvox1', 'porosvox2', 'porosvox3'),
        # stress-strain
        'ezz': ('qd'),
        # q-p
        'pm': ('qd')
}
# show the plot
plot.plot()



def compactionFinished():
    ## set friction angle back to non-zero value and prepare for shear
	# for future contacts change material
	O.materials[0].frictionAngle = radians(finalFric) # radians
	O.materials[0].etaRoll = radians(finalFric)          
    # for existing contacts, set contact friction directly
	for i in O.interactions:
	    i.phys.tangensOfFrictionAngle = tan(radians(finalFric))
	    i.phys.maxRollPl = tan(radians(finalFric))

	# set the current cell configuration to be the reference one
	O.cell.trsf = Matrix3.Identity
	# next time, start shearing specimen
	triax.doneHook = 'startShear()'
	# init strain and stress targets for cyclic loading
	triax.stress_target = targetCSR*(2.*abs(sigmaIso))
	triax.strain_sign = -1
	triax.strain_inc = strain_inc
	triax.strain_target_stage = strain_inc*triax.strain_sign
	triax.strain_target_final = targetStrain
	# set grow damping
	triax.growDamping = triaxGrowDamping  
	# save compaction data and reset 
	triax.growDamping = triaxGrowDamping  
	plot.saveDataTxt(O.tags['d.id'] + '_compaction.txt')


def startShear():
	## change control type for  constant volume, undrained loading
    # strain controlled in xyz
	triax.stressMask = 0
    # avoid hang ups on near-zero strain targets
	if (numpy.isclose(0,triax.strain_target_stage)):
	    triax.strain_target_stage += triax.strain_sign*triax.strain_inc
	triax.goal = (-triax.strain_target_stage*0.5,
                  -triax.strain_target_stage*0.5,
                   triax.strain_target_stage)
	print("Target strain stage: ", triax.goal)
	#triax.maxStrainRate = (.1, .1, .1)
	
	## compute mean and deviatoric stress
	pm = triax.stressTensor.trace()/3.
	qd = triax.stress[2] - (triax.stress[0]+triax.stress[1])/2.
	
	## set the done hook based on stage of test
	# terminate if reached target strain
	if (abs(triax.strain[2]) >= triax.strain_target_final):
	    triax.doneHook = 'triaxFinished()'
	    return
	# handle stress reversals
	if (qd >= triax.stress_target):
	    if (triax.strain_sign > 0):
	        triax.strain_sign = -1.*triax.strain_sign
	elif (qd <= -triax.stress_target):
	    if (triax.strain_sign < 0):
	        triax.strain_sign = -1.*triax.strain_sign
    # increment strain and recurse on startShear hook
	triax.strain_target_stage += triax.strain_sign*triax.strain_inc
	triax.doneHook = 'startShear()'
	 
	# don't wait for stabilization before calling calling next doneHook
	triax.maxUnbalanced = stageMaxUnbalanced
    

def triaxFinished():
	print('Finished')
	# save the plot data
	plot.saveDataTxt(O.tags['d.id'] + '_shear.txt')
	O.pause()

