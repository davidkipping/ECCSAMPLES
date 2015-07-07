#!/bin/bash
#

compiler='gfortran'
flags='-O3'

cflag='-c'
oflag='-o'
f1dir='../f1/'			# directory of the f1 routines
f1depdir='../f1/webroutines/'	# directory of the f1 dependencies
extrasdir='../extras/'		# directory of the extra dependencies
extraflag='-fno-range-check'
modulename='ECCSAMPLES'
progname='example'

# suppress legacy warnings from gfortran
if [ "$compiler" == "gfortran" ]; then
  flags=`echo $flags "-std=legacy -w"`
fi

# clean
rm -f ECCSAMPLES.o f2.o	hypdrvf1.o rkf45.o check_parms.o f21r.o	hypgeo.o \
specfun.o d1mach.o f2gam.o hypgeof1.o writef1.o f1bnl.o	gamma.o	hypstartf1.o \
f1conv.o get_transformation2.o inversebeta.o f1v3.o hyp.o isneg.o
rm -f eccsamplesmod.mod	inversebetamod.mod
rm -f `echo $modulename".o" "ECCSAMPLESmod.mod" $progname`

# f1
`echo $compiler $flags $cflag $f1dir"gamma.f90"`
`echo $compiler $flags $cflag $f1dir"f2.f90"`
`echo $compiler $flags $cflag $f1dir"f2gam.f90"`
`echo $compiler $flags $cflag $f1dir"check_parms.f90"`
`echo $compiler $flags $cflag $f1dir"f1bnl.f90"`
`echo $compiler $flags $cflag $f1dir"f1conv.f90"`
`echo $compiler $flags $cflag $f1dir"f1v3.f90"`
`echo $compiler $flags $cflag $f1dir"hypgeof1.f90"`
`echo $compiler $flags $cflag $f1dir"f21r.f"`
`echo $compiler $flags $cflag $f1dir"get_transformation2.f90"`
`echo $compiler $flags $cflag $f1dir"hypdrvf1.f90"`
`echo $compiler $flags $cflag $f1dir"hypstartf1.f90"`
`echo $compiler $flags $cflag $f1dir"isneg.f90"`
`echo $compiler $flags $cflag $f1dir"writef1.f90"`
# f1 dependencies
if [ "$compiler" == "gfortran" ]; then
  `echo $compiler $flags $extraflag $cflag $f1depdir"d1mach.for"`
else
  `echo $compiler $flags $cflag $f1depdir"d1mach.for"`
fi
`echo $compiler $flags $cflag $f1depdir"hyp.f"`
`echo $compiler $flags $cflag $f1depdir"hypgeo.f"`
`echo $compiler $flags $cflag $f1depdir"rkf45.for"`
# extras
`echo $compiler $flags $cflag $extrasdir"specfun.f90"`
`echo $compiler $flags $cflag $extrasdir"inversebeta.f"`
# esamples
`echo $compiler $flags $cflag $modulename".f90"`
`echo $compiler $flags $oflag $progname "example.f90 gamma.o f2.o \
f2gam.o check_parms.o f1bnl.o f1conv.o f1v3.o hypgeof1.o f21r.o \
get_transformation2.o hypdrvf1.o hypstartf1.o isneg.o writef1.o \
d1mach.o hyp.o hypgeo.o rkf45.o specfun.o inversebeta.o ECCSAMPLES.o"`

# result
found=`ls $progname | wc -l`
if [ "$found" -eq 1 ]; then
  echo 'Compilation successful'
else
  echo 'Compilation failed, suggest turning off warning suppression flags'
fi
