from yade import pack
dimSmall = (0.04, 0.04, 0.04)
dimFull = (0.2, 0.2, 0.2)
sp = pack.SpherePack()

dmin, dmax = 0.001, 0.005
numbins = 25
start, end = numpy.log10(dmin), numpy.log10(dmax)
psdSizes = numpy.logspace(start, end, num=numbins, endpoint=True)
psdCumm = numpy.linspace(0, 1., len(psdSizes))
sp.makeCloud((0, 0, 0), dimSmall, psdSizes=psdSizes, psdCumm=psdCumm, periodic=True)

print(len(sp),sp.dim()) # 1015 particles in this case
sp.cellFill(dimFull)
print(len(sp),sp.dim()) # 126875 particles in this case (is "instantaneous")
#filterSpherePack(...) # if needed
#sp.toSimulation() # takes 30 s on my laptop

