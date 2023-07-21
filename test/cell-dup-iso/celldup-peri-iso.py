# Cell duplication, two-stage compaction, aperiodic only.

from __future__ import print_function
from yade import pack, qt, plot, export
import pandas as pd


# output path
output_path = 'output'

# define inputs
sigmaIso = -1e5
initialFric = 0
finalFric = 36.50
initialEtaRoll = 0.0
finalEtaRoll = 0.12


# define materials
O.materials.append(CohFrictMat(density=2650, young=9.6e7, poisson=.28,
	frictionAngle=radians(initialFric), isCohesive=False, momentRotationLaw=True,
	etaRoll=initialEtaRoll, etaTwist=0, label='granr'))

# define periodic boundaries
O.periodic = True

# define log-linear psd
dmin, dmax = 0.07, 0.13
start, end = numpy.log10(dmin), numpy.log10(dmax)
psdSizes = numpy.logspace(start, end, num=25, endpoint=True) 
psdCumm = numpy.linspace(0, 1., len(psdSizes))

# define cloud and insert packing
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), psdSizes=psdSizes, psdCumm=psdCumm, periodic=True)
sp.toSimulation(material="granr")

# echo cell quantities
print('\nsp.dim(): ', sp.dim())
print('sp.aabb(): ', sp.aabb())
print('sp.cellSize: ', sp.cellSize)
print('sp.center(): ', sp.center())
print("O.cell.refSize: ", O.cell.refSize)
print("O.cell.hSize: ",O.cell.hSize)

# define engines
# compact to 90% goal stress with PeriTriaxController
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
        NewtonIntegrator(damping=.2),
        PeriTriaxController(
                label='triax',
                # specify target values and whether they are strains or stresses
                goal=(sigmaIso*0.90, sigmaIso*0.90, sigmaIso*0.90),
                stressMask=7,
                # type of servo-control
                dynCell=True,
                maxStrainRate=(0.01, 0.01, 0.01),
                # wait until the unbalanced force goes below this value
                maxUnbalanced=1e-6,
                relStressTol=1e-6,
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


def compactionFinished():
    print('\nCompaction Finished')
    print('sp.dim(): ', sp.dim())
    print('sp.aabb(): ', sp.aabb())
    print('sp.cellSize: ', sp.cellSize)
    print('sp.center(): ', sp.center())
    print("O.cell.refSize: ", O.cell.refSize)
    print("O.cell.hSize: ",O.cell.hSize)
    triax.doneHook  = 'dupCell()'

def dupCell():
    print('\nDuplicating Cell')
    print('sp.dim(): ', sp.dim())
    print('sp.aabb(): ', sp.aabb())
    print('sp.cellSize: ', sp.cellSize)
    print('sp.center(): ', sp.center())
    print("O.cell.refSize: ", O.cell.refSize)
    print("O.cell.hSize: ",O.cell.hSize)
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

    triax.doneHook = 'compactionStage2()'
    #O.pause()

def compactionStage2():
    print('\nCompact to full target stress')
    triax.maxUnbalanced=1e-5
    triax.relStressTol=1e-5
    triax.goal=(sigmaIso, sigmaIso, sigmaIso)
    triax.doneHook = 'finished()'

def finished():
    print('Finished')
    plot.saveDataTxt(output_path + '/histories')
    vtkExporter = export.VTKExporter(output_path + '/vtkExporterTesting')
    vtkExporter.exportSpheres(what={'dist':'b.state.pos.norm()'})
    vtkExporter.exportInteractions(what=dict(kn='i.phys.kn'))
    vtkExporter.exportContactPoints(what={'nn': 'i.geom.normal'})
    vtkExporter.exportPeriodicCell()
    triax.doneHook = 'changeBC()'
    O.pause()

