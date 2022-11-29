# Setup Env Var
export INPUT='triggered_convec.inp'
export EXE='../../bin/polar_vortex_threshold.ex'
export MYDIR='./examples/2d.polar_vortex_threshold'

./patch.py
./configure.py $(head $MYDIR/$INPUT | grep configure | cut -d ' ' -f3-)
make clean
make -j8
mkdir -p $MYDIR/data
cd $MYDIR && mpiexec -n 16 $EXE -i $INPUT
../../combine.py -d data -o test
