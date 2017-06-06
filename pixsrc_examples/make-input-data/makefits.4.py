#! /usr/bin/env python

import pyfits
from pyraf import iraf
from subprocess import call

iraf.set(imclobber="yes")
iraf.hedit("ring.fits","CTYPE1",'RA---TAN',add=1,verify=0)
iraf.hedit("ring.fits","CTYPE2",'DEC--TAN',add=1,verify=0)
iraf.hedit("ring.fits","CD1_1",-3.3333333333E-5,add=1,verify=0)
iraf.hedit("ring.fits","CD1_2",0,add=1,verify=0)
iraf.hedit("ring.fits","CD2_1",0,add=1,verify=0)
iraf.hedit("ring.fits","CD2_2",3.3333333333E-5,add=1,verify=0)
iraf.hedit("ring.fits","CRPIX1",26,add=1,verify=0)
iraf.hedit("ring.fits","CRPIX2",26,add=1,verify=0)
iraf.hedit("ring.fits","CRVAL1",161.25,add=1,verify=0)
iraf.hedit("ring.fits","CRVAL2",58,add=1,verify=0)

iraf.noao()
iraf.artdata()

#iraf.gauss(input="ring.fits",output="ring.fits",sigma=0.5, ratio=1, theta=0)

iraf.mknoise(input="ring.fits",background=0,gain=10,rdnoise=0.5)
