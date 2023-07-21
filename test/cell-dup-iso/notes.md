### Blocks
---
```python
sp.dim()        # returns dimensions of packing in terms of aabb(), as 3-tuple
sp.aabb()       # axis-aligned bounding box coordinates, 2 3-tuples
sp.cellSize()   # size of periodic cell
sp.center()     # bounding box center

O.cell.setBox( )       # To update box size
O.cell.refSize    # Base cell vectors (columns of the matrix)
O.cell.hSize      # Reference size of the cell (initial cell vectors)

# save the iso packing
import pandas as pd
df = pd.read_csv("output/sphere_coordinates",sep='\t',names=['x','y','z','r'],skiprows=1)

# get base cell size
dx, dy, dz = O.cell.hSize[0][0], O.cell.hSize[1][1], O.cell.hSize[2][2]

## This works
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


## This doesn't quite work, blows up
# update dimensions of box by 2 rve units in xyz directions
nx, ny, nz = 2, 2, 2
O.cell.setBox((dx*nx,dy*ny,dz*nz))

# duplicate rves in xyz directions
for i in range(nx):
  for j in range(ny):
    for k in range(nz):
	  for _, (x, y, z, r) in df.iterrows():
	    if !(i == 0 and j == 0 and k == 0):
	      O.bodies.append(utils.sphere([x+dx*i,y+dy*j,z+dz*k], r))

```


### Switching Boundaries
---
```python
# wall material (define earlier)
O.materials.append(FrictMat(young=1e8, poisson=0.5, frictionAngle=0, density=0, label='walls'))
O.materials.append(FrictMat(young=1e8, poisson=0.5, frictionAngle=0, density=0, label='penetrometer'))

# interaction loop (define earlier)
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
                        [Ip2_FrictMat_FrictMat_FrictPhys(),
                         Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()],
                        [Law2_ScGeom_FrictPhys_CundallStrack(),
                         Law2_ScGeom6D_CohFrictPhys_CohesionMoment()],
                        label='interactionLoop'),
        ...
]
...

# turn off periodic boundaries
O.periodic=False
# turn off PeriTriaxController
triax.dead=True

## switch to TriaxialStressController
# create and append 4 walls of a cube sized to our mn, mx parameters
mn, mx = Vector3(0, 0, 0), Vector3(1, 1, 1)  # corners of the initial packing
walls = aabbWalls([mn, mx], thickness=0, material='walls')
wallIds = O.bodies.append(walls)
# overwrite triax controller
triax = TriaxialStressController(
        internalCompaction=False,
        stressMask=7,
        goal1=sigmaIso*1.1,
        goal2=sigmaIso*1.1,
        goal3=sigmaIso*1.1,
)
O.engines += [triax]


## Need to update stress-strain histories


# IMPORTANT, periodic needs to be off before generating facets!
# save the cell dim before clearing
cellDim = O.cell.hSize[0][0]
O.periodic=False

## Add the cone penetrometer
# make the cylinder body
cellToProbeRatio = 12.
center = Vector3(cellDim/2, cellDim*1.8, cellDim/2.)
radius = cellDim/cellToProbeRatio
height = cellDim*1.2
rod = geom.facetCylinder(center, radius, height, Quaternion((1, 0, 0), .5*pi))
# make the cone
coneHeight = radius/(numpy.tan(numpy.radians(30)))
cone = geom.facetCone(Vector3(cellDim/2, center[1]-height/2.-coneHeight/2, cellDim/2.), 0, radius, coneHeight, Quaternion((1, 0, 0), .5*pi))
# import the geoms
rodIds = O.bodies.append(rod)
coneIds = O.bodies.append(cone)


## Add penetromer kinematics
trans = TranslationEngine(ids=rod, translationAxis=[0, -1, 0], velocity=0.1)
O.engines += [trans]

## Also update the triaxController to fix vertical walls, and only control lateral stress



```