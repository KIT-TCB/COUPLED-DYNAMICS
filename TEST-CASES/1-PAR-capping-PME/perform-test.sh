# tests the FO hamiltonian along a given trajectory (trajectories are not numerically identical even in normal MD simulations with and without MPI)
# compiling without MPI gives slightly different TB_HAMILTONAINA.xvg than with.


# failure of the test can be caused by:
#   different slater koster files
#   FO construction routine
#   side-chain capping routine
#   calculation of ESP of MM environment (PME)  


progpath=/home/kit/ipc/yc6285/COUPLED-DYNAMICS-build/src/kernel/



#check which program we can test
dir=`pwd`
tests=""
if [ -e $progpath/mdrun_mpi ] ;then tests="mpi-multinode mpi-singlenode"; fi
if [ -e $progpath/mdrun ] ;then tests=$tests" non-mpi"; fi



echo "Testing program in $progpath"
for run in $tests
do
  rm $dir/$run -rf 
  mkdir $dir/$run
  cd $dir/$run
  cp $dir/rerun.trr $dir/charge-transfer.dat $dir/mol.spec $dir/ct.tpr .

  if [ $run == "non-mpi" ]
    then $progpath/mdrun -pd -rerun rerun.trr -s ct.tpr > out.dat 2> /dev/null
  elif [ $run == "mpi-multinode" ]
    then mpirun -np 3 "$progpath/mdrun_mpi -pd -rerun rerun.trr -s ct.tpr" > out.dat 2> /dev/null
  else 
    mpirun -np 1 "$progpath/mdrun_mpi -pd -rerun rerun.trr -s ct.tpr" > out.dat 2> /dev/null
  fi

  if [ -e $dir/$run/TB_HAMILTONIAN.xvg ] 
  then
    awk '{for (i=1; i<=NF; i++ ){ printf "%7.3f ", $i} ; print ""}' TB_HAMILTONIAN.xvg > significant_TB_HAMILTONIAN.xvg
    differences=`diff $dir/expected_TB_HAMILTONIAN.xvg $dir/$run/significant_TB_HAMILTONIAN.xvg | wc -l` 
    if [ $differences != 0 ]; then echo "FAILED TEST   $run" ;else echo "PASSED TEST   $run" ; fi
  else 
    echo "FAILED TEST   $run"
  fi
done



