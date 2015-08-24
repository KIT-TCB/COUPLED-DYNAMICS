# tests the adaptive QM zone, but only the initial building of an optimized QM zone around the lowest energy site. too few steps to see a moving QM zone.
# note: mpi version doesn't work yet on multiple nodes

# failure of the test (if previous test were passed) can be caused by:
#   adapt_QMzone


progpath=/home/kit/ipc/yc6285/COUPLED-DYNAMICS-build/src/kernel/



#check which program we can test
dir=`pwd`
tests=""
if [ -e $progpath/mdrun_mpi ] ;then tests="mpi-singlenode"; fi
if [ -e $progpath/mdrun ] ;then tests=$tests" non-mpi"; fi



echo "Testing program in $progpath"
for run in $tests
do
  rm $dir/$run/* -f
  cd $dir/$run
  cp $dir/charge-transfer.dat $dir/mol.spec $dir/ct.tpr .

  if [ $run == "non-mpi" ]
    then $progpath/mdrun -pd -s ct.tpr -nsteps 1 > out.dat 2> /dev/null
  elif [ $run == "mpi-multinode" ]
    then mpirun -np 3 "$progpath/mdrun_mpi -pd -s ct.tpr -nsteps 1" > out.dat 2> /dev/null
  else 
    mpirun -np 1 "$progpath/mdrun_mpi -pd -s ct.tpr -nsteps 1" > out.dat 2> /dev/null
  fi
  
  awk '{for (i=1; i<=NF; i++ ){ printf "%7.3f ", $i} ; print ""}' TB_HAMILTONIAN.xvg > significant_TB_HAMILTONIAN.xvg
  sed -i 's/-0.000/ 0.000/g' significant_TB_HAMILTONIAN.xvg
  differences=`diff $dir/expected_TB_HAMILTONIAN.xvg $dir/$run/significant_TB_HAMILTONIAN.xvg | wc -l` 
  if [ $differences != 0 ]; then echo "FAILED TEST   $run" ;else echo "PASSED TEST   $run" ; fi
done



