from __future__ import print_function
from yade import pack, qt, plot, export
import pandas as pd


output_path = 'output'

sigmaIso = -1e5
initialFric = 0
finalFric = 36.50
initialEtaRoll = 0.0
finalEtaRoll = 0.12

## define materials
O.materials.append(CohFrictMat(
        density=2650, # Density [kg/m3]
        young=9.6e7, # Particle modulus [Pa]
        poisson=.28, # Ks/Kn ratio
        frictionAngle=radians(initialFric), # Local friction [rad]
        isCohesive=False, # Turn off adhesion
        momentRotationLaw=True, # Turn on rotational stiffness
        etaRoll=initialEtaRoll, # Rotational friction [rad]
        etaTwist=0, # Turn off twisting
        label="granr" # Material label
        ))

O.periodic = True
sp = pack.SpherePack()
if 0:
    ## uniform distribution
    sp.makeCloud((0, 0, 0), (4, 4, 4), rMean=.1, rRelFuzz=.3, periodic=True)
else:
    ## log-linear distribution
    dmin, dmax = 0.07, 0.13
    numbins = 25
    start, end = numpy.log10(dmin), numpy.log10(dmax)
    psdSizes = numpy.logspace(start, end, num=numbins, endpoint=True) 
    psdCumm = numpy.linspace(0, 1., len(psdSizes))
    sp.makeCloud((0, 0, 0), (1, 1, 1), psdSizes=psdSizes, psdCumm=psdCumm, periodic=True)
    print(sp.aabb()) # Print axes aligned bounding box coordinates


print('sp.dim(): ', sp.dim())
print('sp.aabb(): ', sp.aabb())
print('sp.cellSize: ', sp.cellSize)
print('sp.center(): ', sp.center())
print("O.cell.refSize: ", O.cell.refSize)
print("O.cell.hSize: ",O.cell.hSize)


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
    print('Compaction Finished')
    doneHook='dupCell()'
    print('sp.dim(): ', sp.dim())
    print('sp.aabb(): ', sp.aabb())
    print('sp.cellSize: ', sp.cellSize)
    print('sp.center(): ', sp.center())
    print("O.cell.refSize: ", O.cell.refSize)
    print("O.cell.hSize: ",O.cell.hSize)
    triax.doneHook  = 'dupCell()'

def dupCell():
    print('Duplicating Cell')

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
    print('Recompacting second time after duplication')
    triax.maxUnbalanced=1e-6
    triax.relStressTol=1e-6
    triax.goal=(sigmaIso*1.1, sigmaIso*1.1, sigmaIso*1.1)
    triax.doneHook = 'finished()'

def finished():
    print('Finished')
    plot.saveDataTxt(output_path + '/histories')
    vtkExporter = export.VTKExporter(output_path + '/vtkExporterTesting')
    vtkExporter.exportSpheres(what={'dist':'b.state.pos.norm()'})
    vtkExporter.exportInteractions(what=dict(kn='i.phys.kn'))
    vtkExporter.exportContactPoints(what={'nn': 'i.geom.normal'})
    vtkExporter.exportPeriodicCell()
    O.pause()

