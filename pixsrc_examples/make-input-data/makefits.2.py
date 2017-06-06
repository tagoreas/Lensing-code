#! /usr/bin/env python

import pyfits
from pyraf import iraf
from subprocess import call

iraf.set(imclobber="yes")
iraf.hedit("merger.fits","CTYPE1",'RA---TAN',add=1,verify=0)
iraf.hedit("merger.fits","CTYPE2",'DEC--TAN',add=1,verify=0)
iraf.hedit("merger.fits","CD1_1",-1.66666666666E-5,add=1,verify=0)
iraf.hedit("merger.fits","CD1_2",0,add=1,verify=0)
iraf.hedit("merger.fits","CD2_1",0,add=1,verify=0)
iraf.hedit("merger.fits","CD2_2",1.66666666666E-5,add=1,verify=0)
iraf.hedit("merger.fits","CRPIX1",51,add=1,verify=0)
iraf.hedit("merger.fits","CRPIX2",51,add=1,verify=0)
iraf.hedit("merger.fits","CRVAL1",161.25,add=1,verify=0)
iraf.hedit("merger.fits","CRVAL2",58,add=1,verify=0)

iraf.noao()
iraf.artdata()

iraf.gauss(input="merger.fits",output="merger.fits",sigma=2., ratio=0.7, theta=35)

iraf.mknoise(input="merger.fits",background=0,gain=10,rdnoise=1)
