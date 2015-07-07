COMPILING

In the ECCSAMPLES subdirectory, compile using...
./make.sh
then execute using...
./example_call

INPUTS

- Change alpha and beta shape parameters describing the eccentricity's
  underlying distribution (assumed to be a Beta distribution) in the
  example_call.f90 file.

- Switch on/off whether the planet is known to be transiting via the
  transit logical in the example_call.f90 file.

- Change the number of sample of output using the n integer parameter in
  example_call.f90

OUTPUTS

- ECCSAMPLES subroutine returns an eccentricity and an argument of periastron
  back to example_call.f90. These are then added to a vector for each of length 
  n within example_call.f90. Finally, they are outputted in a plain text file, 
  ECCSAMPLES.dat by the example_call.f90 code.
