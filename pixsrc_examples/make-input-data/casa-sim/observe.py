
default('simobserve')
projname='merger.uv'
project = projname
skymodel = '../../data/merger.uv.fits'
indirection = 'J2000 03h30m00 -28d00m00'
inbright = '0.0001'
incell = '0.01arcsec'
incenter = '245GHz'
inwidth = '2GHz'
integration = '1s'
antennalist =  "alma.out20.cfg"
totaltime = '20000s'
overwrite = T
simobserve()

default('simanalyze')
projname='merger.uv'
project = projname
image = T
vis = project+'.alma.out20.noisy.ms'
print vis
imsize = [600,600]
cell = '0.01arcsec'
#imsize = [192, 192]
threshold = '1e-9Jy'
#niter = 10000
niter = 30000
weighting = 'natural'
analyze = F
overwrite = T
verbose = T
simanalyze()
