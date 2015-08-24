for test in */
do
echo "testing case $test" 
cd $test
./perform-test.sh
cd ../
done
