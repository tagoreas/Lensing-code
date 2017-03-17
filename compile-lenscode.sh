#!/bin/bash

BASEDIR=`pwd`

# external libraries
libs='-lgsl -L'$BASEDIR'/lib -lpixsrc'

# attempt to get number of threads to be used during lens modelling (parallel MCMC sampling)
numlibs=$(cat /proc/cpuinfo | grep -i processor | awk '{print $NF+1}' | tail -n1)

# cleanup
mkdir -p lib
rm -f lib/libeval-tps-avd-hmpr.so lib/libstple_lenscalc.*.so

# compile code and create library
names='eval-tps-avd-hmpr stple_lenscalc'
cd src/lenscode
rm -f *.o
for name in $names; do
    g++ -Wno-unused-variable -Wno-unused-but-set-variable -O3 -g -fPIC -Wall -c $name.cpp
    g++ -Wno-unused-variable -Wno-unused-but-set-variable -O3 -g -fPIC -shared -o lib$name.so -Wl,-Bdynamic $name.o $libs
    mv lib$name.so ../../lib/
done
cd ../../

# copy library so multiple threads can have their own copy (a hack, but it works with python's ctypes)
let numlibs=numlibs-1
cd lib
for i in $(seq 0 $numlibs); do 
    ln -s libstple_lenscalc.so libstple_lenscalc.$i.so
done
cd ..
