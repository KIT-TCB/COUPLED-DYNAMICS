# tests the bridge mediated coupling (TDA) 

# failure of the test (if previous test were passed) can be caused by:
#   TDA routine


progpath=/home/kit/ipc/yc6285/COUPLED-DYNAMICS-build/src/kernel/



#check which program we can test
dir=`pwd`
tests=""
if [ -e $progpath/mdrun_mpi ] ;then tests="mpi-multinode mpi-singlenode"; fi
if [ -e $progpath/mdrun ] ;then tests=$tests" non-mpi"; fi



echo "Testing program in $progpath"
for run in $tests
do
  rm $dir/$run/* -f 
  cd $dir/$run
  cp $dir/charge-transfer.dat $dir/et.spec  $dir/etf2.spec $dir/ct.tpr .

  if [ $run == "non-mpi" ]
    then $progpath/mdrun -pd -s ct.tpr > out.dat 2> /dev/null
  elif [ $run == "mpi-multinode" ]
    then mpirun -np 3 "$progpath/mdrun_mpi -pd -s ct.tpr " > out.dat 2> /dev/null
  else 
    mpirun -np 1 "$progpath/mdrun_mpi -pd -s ct.tpr " > out.dat 2> /dev/null
  fi
  
  differences=`diff $dir/expected_CT_TDA.xvg $dir/$run/CT_TDA.xvg | wc -l` 
  if [ $differences != 0 ]; then echo "FAILED TEST   $run" ;else echo "PASSED TEST   $run" ; fi
done



