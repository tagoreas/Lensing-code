#! /usr/bin/env python

import pyfits
from pyraf import iraf
from subprocess import call

iraf.set(imclobber="yes")
iraf.hedit("gauss.fits","CTYPE1",'RA---TAN',add=1,verify=0)
iraf.hedit("gauss.fits","CTYPE2",'DEC--TAN',add=1,verify=0)
iraf.hedit("gauss.fits","CD1_1",-5.55555555555E-5,add=1,verify=0)
iraf.hedit("gauss.fits","CD1_2",0,add=1,verify=0)
iraf.hedit("gauss.fits","CD2_1",0,add=1,verify=0)
iraf.hedit("gauss.fits","CD2_2",5.55555555555E-5,add=1,verify=0)
iraf.hedit("gauss.fits","CRPIX1",16,add=1,verify=0)
iraf.hedit("gauss.fits","CRPIX2",16,add=1,verify=0)
iraf.hedit("gauss.fits","CRVAL1",161.25,add=1,verify=0)
iraf.hedit("gauss.fits","CRVAL2",58,add=1,verify=0)

iraf.noao()
iraf.artdata()

iraf.gauss(input="gauss.fits",output="gauss.fits",sigma=2., ratio=0.7, theta=35)

iraf.mknoise(input="gauss.fits",background=0,gain=10,rdnoise=1)
