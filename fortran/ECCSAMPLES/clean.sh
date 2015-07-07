#!/bin/bash
#

modulename='ECCSAMPLES'
progname='example'

# clean
rm -f ECCSAMPLES.o f2.o	hypdrvf1.o rkf45.o check_parms.o f21r.o	hypgeo.o \
specfun.o d1mach.o f2gam.o hypgeof1.o writef1.o f1bnl.o	gamma.o	hypstartf1.o \
f1conv.o get_transformation2.o inversebeta.o f1v3.o hyp.o isneg.o
rm -f eccsamplesmod.mod	inversebetamod.mod
rm -f *~
rm -f `echo $modulename".o" "ECCSAMPLESmod.mod" $progname`
