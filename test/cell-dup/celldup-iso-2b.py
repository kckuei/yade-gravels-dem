from __future__ import print_function
from yade import pack, qt, plot, export


sigmaIso = -1e5

# define materials
O.materials.append(CohFrictMat(density=2650,young=9.6e7, poisson=.28,frictionAngle=0,
    isCohesive=False,momentRotationLaw=True, etaRoll=0,etaTwist=0,label="granr"))
O.materials.append(CohFrictMat(young=1e8, poisson=0.28, frictionAngle=0, 
    density=0, isCohesive=False, momentRotationLaw=False, etaTwist=0, label='walls'))

# set periodic boundaries
O.periodic = False

# define log-linear psd
dmin, dmax = 0.07, 0.13
start, end = numpy.log10(dmin), numpy.log10(dmax)
psdSizes = numpy.logspace(start, end, num=25, endpoint=True) 
psdCumm = numpy.linspace(0, 1., len(psdSizes))

# make packing
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), psdSizes=psdSizes, psdCumm=psdCumm, periodic=True)
sp.toSimulation(material="granr")

O.engines = [
        ForceResetter(label='resetter'),
        InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()],
                    label='collider'),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(),Ig2_Wall_Sphere_ScGeom()],
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
        PeriTriaxController(
            label='triax',
            dead=False,
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
            doneHook='changeBC()'
        ),
        TriaxialStressController(label='triaxw',dead=True),
        NewtonIntegrator(damping=.2),
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
         pmw=triaxw.meanStress,
         Etot=O.energy.total(),
         **O.energy
 )

# enable energy tracking in the code
O.trackEnergy = True

# define what to plot
plot.plots = {
        'i': ('unbalanced',),
        'i ': ('sxx', 'syy', 'szz', 'pmw'),
        ' i': ('exx', 'eyy', 'ezz'),
        # energy plot
        ' i ': (O.energy.keys, None, 'Etot'),
}
# show the plot
plot.plot()

def changeBC():
    # turn off periodic boundaries
    O.periodic=False
    
    # deactivate PeriTriaxController
    triax.dead=True

    # Update engines for walls
    collider.boundDispatcher.functors += [Bo1_Box_Aabb()]
    interactionLoop.geomDispatcher.functors += [Ig2_Wall_Sphere_ScGeom()]

    ## switch to wall TriaxialStressController
    # create and append 4 walls of a cube sized to our mn, mx parameters
    mn, mx = yade._utils.aabbExtrema()  # corners of the final aperiodic packing
    walls = aabbWalls([mn, mx], thickness=0, material='walls')
    wallIds = O.bodies.append(walls)
    # update wall triax goals, assign aabbWall Ids to triax..
    triaxw.internalCompaction=False
    triaxw.max_vel=1
    triaxw.stressMask=7
    triaxw.goal1=triaxw.goal2=triaxw.goal3=sigmaIso
    triaxw.wall_left_id=wallIds[0]
    triaxw.wall_right_id=wallIds[1]
    triaxw.wall_bottom_id=wallIds[2]
    triaxw.wall_top_id=wallIds[3]
    triaxw.wall_back_id=wallIds[4]
    triaxw.wall_front_id=wallIds[5]
    # activate wall triax
    triaxw.dead=False    



