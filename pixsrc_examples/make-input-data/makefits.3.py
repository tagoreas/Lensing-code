#! /usr/bin/env python

import pyfits
from pyraf import iraf
from subprocess import call

iraf.set(imclobber="yes")
iraf.hedit("cross.fits","CTYPE1",'RA---TAN',add=1,verify=0)
iraf.hedit("cross.fits","CTYPE2",'DEC--TAN',add=1,verify=0)
iraf.hedit("cross.fits","CD1_1",-3.3333333333E-5,add=1,verify=0)
iraf.hedit("cross.fits","CD1_2",0,add=1,verify=0)
iraf.hedit("cross.fits","CD2_1",0,add=1,verify=0)
iraf.hedit("cross.fits","CD2_2",3.3333333333E-5,add=1,verify=0)
iraf.hedit("cross.fits","CRPIX1",26,add=1,verify=0)
iraf.hedit("cross.fits","CRPIX2",26,add=1,verify=0)
iraf.hedit("cross.fits","CRVAL1",161.25,add=1,verify=0)
iraf.hedit("cross.fits","CRVAL2",58,add=1,verify=0)

iraf.noao()
iraf.artdata()

iraf.gauss(input="cross.fits",output="cross.fits",sigma=2., ratio=0.7, theta=35)

iraf.mknoise(input="cross.fits",background=0,gain=10,rdnoise=1)
