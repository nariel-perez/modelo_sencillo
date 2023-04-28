

outx="microgels-mc.x"
#outx="otro.x"

rm -f -r $outx*


all_files="parameters.f95 minpack.f90 fvn_common.f90 fvn_interpol.f90 simple_microgel.f95 nvt_mc.f95 microgels-mc.f95"



read -r -p "Debugging? [y/N] " response

if [[ "$response" =~ ^(yes|y)$ ]];then

#flags="-O -Wall -g -fcheck=all -fbacktrace -ffpe-trap=overflow"
flags="-O -Wall -g -fcheck=all -fbacktrace"

else

#flags="-O3 -Wall"
flags="-O3 -Wall -Wno-unused-dummy-argument -Wno-unused-variable"

fi



gfortran $flags $all_files -o $outx

rm -f *.mod
