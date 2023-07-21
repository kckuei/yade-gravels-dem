from yade import pack
r = 0.005
num = 20
dimSmall = (0.04, 0.04, 0.04)
dimFull = (1, 1, 1)
sp = pack.SpherePack()
sp.makeCloud((0,0,0),dimSmall,r,num=num,periodic=True,seed=1)
print(len(sp),sp.dim()) # 20 particles in this case
sp.cellFill(dimFull)
print(len(sp),sp.dim()) # 312500 particles in this case (is "instantaneous")
#filterSpherePack(...) # if needed
#sp.toSimulation() # takes 30 s on my laptop
