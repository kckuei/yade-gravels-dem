from __future__ import print_function
from yade import pack, qt, plot, export
import pandas as pd

# define output path folder
output_path = 'output'

# define some inputs
sigmaIso = -1e5
initialFric = 0
finalFric = 36.50
initialEtaRoll = 0.0
finalEtaRoll = 0.12

# define materials
O.materials.append(CohFrictMat(density=2650,young=9.6e7,poisson=.28,
    frictionAngle=radians(initialFric),isCohesive=False,momentRotationLaw=True,
    etaRoll=initialEtaRoll,etaTwist=0,label='spheres'))
O.materials.append(CohFrictMat(young=1e8, poisson=0.28, frictionAngle=0, 
    density=0, isCohesive=False, momentRotationLaw=False, etaTwist=0, label='walls'))

# define periodic boundaries
O.periodic = True

# define log-linear psd
dmin, dmax = 0.07, 0.13
start, end = numpy.log10(dmin), numpy.log10(dmax)
psdSizes = numpy.logspace(start, end, num=25, endpoint=True) 
psdCumm = numpy.linspace(0, 1., len(psdSizes))

# make cloud and insert packing
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), psdSizes=psdSizes, psdCumm=psdCumm, periodic=True)
sp.toSimulation(material="spheres")

triaxw = TriaxialStressController(dead=True)

# define engines
O.engines = [
        ForceResetter(label='resetter'),
        InsertionSortCollider([Bo1_Sphere_Aabb()],label='collider'),
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
            doneHook='duplicateCell()'
        ),
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
        pmw=triaxw.meanStress,
        exx=triax.strain[0],
        eyy=triax.strain[1],
        ezz=triax.strain[2],
        n0=yade._utils.porosity(),
        n1=yade._utils.voxelPorosity(200,*yade._utils.aabbExtrema()),
        Zm=utils.avgNumInteractions(skipFree=True),
        Cn=utils.avgNumInteractions(skipFree=False),
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

def duplicateCell():
    # Duplicate in x-dim
    # save iso packing and reload to dataframe
    export.text(output_path + "/sphere_coordinates")
    df = pd.read_csv("output/sphere_coordinates",sep='\t',names=['x','y','z','r'],skiprows=1)
    # get base cell size
    dx, dy, dz = O.cell.hSize[0][0], O.cell.hSize[1][1], O.cell.hSize[2][2]
    # update dimensions of box by 1 rve unit
    O.cell.setBox((dx*2,dy,dz)) 
    # duplicate along x dimension
    for _, (x, y, z, r) in df.iterrows():
        O.bodies.append(utils.sphere([x+dx,y,z], r))

    # Duplicate in y-dim
    # save iso packing and reload to dataframe
    export.text(output_path + "/sphere_coordinates")
    df = pd.read_csv("output/sphere_coordinates",sep='\t',names=['x','y','z','r'],skiprows=1)
    # get base cell size
    dx, dy, dz = O.cell.hSize[0][0], O.cell.hSize[1][1], O.cell.hSize[2][2]
    # update dimensions of box by 1 rve unit
    O.cell.setBox((dx,dy*2,dz)) 
    # duplicate along y dimension
    for _, (x, y, z, r) in df.iterrows():
        O.bodies.append(utils.sphere([x,y+dy,z], r))

    # Duplicate in z-dim
    # save iso packing and reload to dataframe
    export.text(output_path + "/sphere_coordinates")
    df = pd.read_csv("output/sphere_coordinates",sep='\t',names=['x','y','z','r'],skiprows=1)
    # get base cell size
    dx, dy, dz = O.cell.hSize[0][0], O.cell.hSize[1][1], O.cell.hSize[2][2]
    # update dimensions of box by 1 rve unit
    O.cell.setBox((dx,dy,dz*2)) 
    # duplicate along z dim
    for _, (x, y, z, r) in df.iterrows():
        O.bodies.append(utils.sphere([x,y,z+dz], r))

    triax.doneHook = 'stage2Compaction'
    O.pause()

def stage2Compaction():
    triax.maxUnbalanced=1e-5
    triax.relStressTol=1e-5
    triax.goal=(sigmaIso*1.1, sigmaIso*1.1, sigmaIso*1.1)
    triax.doneHook = 'finishedCompaction'

def finishedCompaction():
    plot.saveDataTxt(output_path + '/histories')
    vtkExporter = export.VTKExporter(output_path + '/vtkExporterTesting')
    vtkExporter.exportSpheres(what={'dist':'b.state.pos.norm()'})
    vtkExporter.exportInteractions(what=dict(kn='i.phys.kn'))
    vtkExporter.exportContactPoints(what={'nn': 'i.geom.normal'})
    vtkExporter.exportPeriodicCell()
    triax.doneHook = 'changeBC()'
    #O.pause()

def changeBC():
    print("Change boundary condition")
    
    # save the cell dim which gets wiped when made aperiodic
    cellDim = O.cell.hSize[0][0]
    
    # turn off periodic boundaries
    O.periodic=False
    
    # turn off/deactivate PeriTriaxController
    triax.dead=True

    # Update engines for walls
    collider.boundDispatcher.functors += [Bo1_Box_Aabb()]
    interactionLoop.geomDispatcher.functors += [Ig2_Box_Sphere_ScGeom(),Ig2_Wall_Sphere_ScGeom()]
    
    ## switch to wall TriaxialStressController
    # create and append 4 walls of a cube sized to our mn, mx parameters
    mn, mx = yade._utils.aabbExtrema()  # corners of the final aperiodic packing
    walls = aabbWalls([mn, mx], thickness=0, material='walls')
    wallIds = O.bodies.append(walls)
    # activate triax controller, update goals, assign aabbWall Ids to triax..
    triaxw.internalCompaction=False
    triaxw.max_vel=1
    triaxw.stressDamping=0.25
    triaxw.stressMask=7
    triaxw.goal1=triaxw.goal2=triaxw.goal3=-sigmaIso*1.1
    triaxw.wall_left_id=wallIds[0]
    triaxw.wall_right_id=wallIds[1]
    triaxw.wall_bottom_id=wallIds[2]
    triaxw.wall_top_id=wallIds[3]
    triaxw.wall_back_id=wallIds[4]
    triaxw.wall_front_id=wallIds[5]
    triaxw.dead=False

    O.pause()



