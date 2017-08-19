
**What does this code do?**

The code here can be compiled into *three* distinct 
libraries (See the INSTALL file).
They do various things:
- Modeling quasar time delays using triaxial lenses.
- Modeling central velocity dispersions of
  triaxial lenses.
- Modeling lensed, extended source emission.

One of these (*libstple_lenscalc.so*) is used for
calculating deflections, time delays, solving the
lens equations, and calculating a chi^2 value.

The second (*libeval-tps-avd-hmpr.so*) is used for
calculating aperture velocity dispersions and two-
dimensional half light radii.

The third library (*libpixsrc.so*) is used for modeling
the extended source emission of a lensed object.

**How do I use it?**

"fit-lenses.py" does MCMC sampling of the lens
parameter space. It requires the libstple_lenscalc.so
and libeval-tps-avd-hmpr.so libraries, but libpixsrc.so
is (sort of) optional.

Additionally, libpixsrc.so can be used by
"fit-lenses.py" or by Chuck Keeton's gravlens/lensmodel
code.

Examples of using pixsrc can be found in the
"pixsrc_examples" directory.
pixsrc documentation can be found there as well.
The examples require lensmodel to run. 
