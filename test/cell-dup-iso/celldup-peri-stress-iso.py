# Performs cell duplication with two-stage periodic isotropic compaction
# with PeriTriaxialController followed by aperiodic walled isotropic compaction
# with StressTriaxController

from __future__ import print_function
from yade import pack, qt, plot, export
import pandas as pd


# output path
output_path = 'output'

# packing stress
pk_sigmaIso1 = -1e5
pk_sigmaIso2 = -1e5

# cell
pk_mn = (0, 0, 0)
pk_mx = (1, 1, 1)

# sphere props
dmin, dmax = 0.07, 0.13
sp_density = 2650
sp_young = 9.67e7
sp_poisson = 0.28
sp_fric_init = 0.0		# deg
sp_fric_final = 36.50		# deg
sp_fric_rr_init = 0.0		# rad
sp_fric_rr_final = 0.12		# rad

# wall props
wl_fric = 0.0			# deg
wl_fric_rr = 0.0		# rad

# create a dict to store vars for global access
myvars = {}


# define materials
O.materials.append(CohFrictMat(density=sp_density, young=sp_young, poisson=sp_poisson,
	frictionAngle=radians(sp_fric_init), isCohesive=False, momentRotationLaw=True,
	etaRoll=sp_fric_rr_init, etaTwist=0, label='granr'))

# define periodic boundaries
O.periodic = True

# define log-linear psd
start, end = numpy.log10(dmin), numpy.log10(dmax)
psdSizes = numpy.logspace(start, end, num=25, endpoint=True) 
psdCumm = numpy.linspace(0, 1., len(psdSizes))

# define cloud and insert packing
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), psdSizes=psdSizes, psdCumm=psdCumm, periodic=True, seed=1)
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
        NewtonIntegrator(damping=.2, label='newton'),
        TranslationEngine(label='translation', dead=True),
        PeriTriaxController(
                label='triax',
                # specify target values and whether they are strains or stresses
                goal=(pk_sigmaIso1, pk_sigmaIso1, pk_sigmaIso1),
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
        PyRunner(command='checkTriaxW()', iterPeriod=1000, 
        	label='pyrnr_checkWallTriax', dead=True),
        PyRunner(command='calm()', iterPeriod=10,
        	label='pyrnr_calm', dead=True),
]

# init wall triax controller
triaxw = TriaxialStressController()


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
    	nw=yade.utils.voxelPorosityTriaxial(triaxw,200),
    	orMax=yade._utils.maxOverlapRatio(),
    	Zm=utils.avgNumInteractions(skipFree=True),
    	Cn=utils.avgNumInteractions(skipFree=False),
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
        '  i': ('n0','n1'),
        '  i ': ('Zm','Cn'),
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
    print('\nDuplicate Cell and Compact to Full Stress')
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

    triax.maxUnbalanced=1e-5
    triax.relStressTol=1e-5
    triax.goal=(pk_sigmaIso2, pk_sigmaIso2, pk_sigmaIso2)
    triax.doneHook = 'finished()'
    #O.pause()

def finished():
    print('Finished')
    plot.saveDataTxt(output_path + '/histories')
    vtkExporter = export.VTKExporter(output_path + '/vtkExporterTesting')
    vtkExporter.exportSpheres(what={'dist':'b.state.pos.norm()'})
    vtkExporter.exportInteractions(what=dict(kn='i.phys.kn'))
    vtkExporter.exportContactPoints(what={'nn': 'i.geom.normal'})
    vtkExporter.exportPeriodicCell()
    triax.doneHook = 'switchController()'
    #O.pause()
    
def switchController():
    print("Change boundary conditions to wall triax")
    
    # save the cell dim which gets wiped when made aperiodic
    cellDim = O.cell.hSize[0][0]
    
    # turn off periodic boundaries
    O.periodic=False
    
    # turn off/deactivate PeriTriaxController
    triax.dead=True
    
    # add material for wall
    O.materials.append(CohFrictMat(young=1e8, poisson=0.28, frictionAngle=0, 
    	density=0, isCohesive=False, momentRotationLaw=False, etaTwist=0, label='walls'))
    
    # update engines for walls
    collider.boundDispatcher.functors += [Bo1_Box_Aabb()]
    interactionLoop.geomDispatcher.functors += [Ig2_Box_Sphere_ScGeom(),Ig2_Wall_Sphere_ScGeom()]

    # create and append 4 walls of a cube sized to our mn, mx parameters
    mn, mx = yade._utils.aabbExtrema()  # corners of the final aperiodic packing
    walls = aabbWalls([mn, mx], thickness=0, material='walls')
    wallIds = O.bodies.append(walls)
    
    # update wall triax settings, goals, and assign aabbWall Ids
    triaxw.internalCompaction=False
    triaxw.max_vel=0.005		# needs to be sufficiently small to avoid osc
    triaxw.stressDamping=0.6		# Use high damping for now, change later
    triaxw.stressMask=7
    triaxw.goal1=triaxw.goal2=triaxw.goal3=pk_sigmaIso2
    triaxw.wall_left_id=wallIds[0]
    triaxw.wall_right_id=wallIds[1]
    triaxw.wall_bottom_id=wallIds[2]
    triaxw.wall_top_id=wallIds[3]
    triaxw.wall_back_id=wallIds[4]
    triaxw.wall_front_id=wallIds[5]
    
    # add wall triax to engine
    O.engines +=[triaxw]
    
    # activate triaxw
    triaxw.dead=False
    
    ## reassign final friction values
    # for future contacts change material
    O.materials[0].frictionAngle = radians(sp_fric_final)
    O.materials[0].etaRoll = sp_fric_rr_final
    # for existing contacts, set contact friction directly
    for i in O.interactions:
    	i.phys.tangensOfFrictionAngle = radians(sp_fric_final)
    	i.phys.maxRollPl = sp_fric_rr_final
    
    	
    ## Add the cone penetrometer
    # make the cylinder body
    cellToProbeRatio = 12.
    center = Vector3(cellDim/2, cellDim*1.8, cellDim/2.)
    radius = cellDim/cellToProbeRatio
    height = cellDim*1.2
    rodIds = geom.facetCylinder(center, radius, height, Quaternion((1, 0, 0), .5*pi))
    # make the cone
    coneHeight = radius/(numpy.tan(numpy.radians(30)))
    coneIds = geom.facetCone(Vector3(cellDim/2, center[1]-height/2.-coneHeight/2, cellDim/2.), 0, radius, coneHeight, Quaternion((1, 0, 0), .5*pi))
    # import the geoms
    rodIds = O.bodies.append(rodIds)
    coneIds = O.bodies.append(coneIds)
    
    myvars['rodIds'] = rodIds
    myvars['coneIds'] = coneIds
    
    
    # re-establish wall stresses with high damping and calming sphere velocities
    # to preserve state
    newton.damping = 0.8		# Increasing damping temp
    pyrnr_calm.dead = False
    pyrnr_checkWallTriax.dead = False
    
    
def checkTriaxW():
    fac = abs((pk_sigmaIso2-triaxw.meanStress)/pk_sigmaIso2)
    if unbalancedForce() < 0.0018 and fac < 0.001:
    	# deactivate pyrunners
        pyrnr_calm.dead = True
        pyrnr_checkWallTriax.dead = True
        
        # reset damping
        newton.damping = 0.2
        triaxw.stressDamping = 0.2
        
        # add material for rod
        O.materials.append(CohFrictMat(young=1e8, poisson=0.28, frictionAngle=0, density=0, 
        	isCohesive=False, momentRotationLaw=False, etaTwist=0, label='rod'))
        	
        # update engines for facets
        collider.boundDispatcher.functors += [Bo1_Facet_Aabb()]
        interactionLoop.geomDispatcher.functors += [Ig2_Facet_Sphere_ScGeom()]

	    
       	# activate translation engine for rod
        translation.dead = False
        translation.ids = myvars['coneIds'] + myvars['rodIds']
        translation.translationAxis = [0,-1,0]
        translation.velocity = 0.01
        
        O.pause()
    else:
        # Matybe adjust iterPeriod of calm as we approach target stress?
        # or reduce the damping?
        newton.damping = 0.6*fac + 0.2
        triaxw.stressDamping = 0.6*fac + 0.2

def calm():
    # zero linear and angular velocoties
    for b in O.bodies:
    	b.state.vel = Vector3(0.,0.,0.)
    	b.state.angVel = Vector3(0.,0.,0.)



