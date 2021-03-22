#!/bin/bash
# a bash script to run testing and bechmark between boost graph library and my sssp library
#both use the same input format, and output only result

#build boost-graph first

cd boost-graph
mkdir build
cd build 
cmake ..
make -j
cp boost_graph ../../
cd ../../
echo "build boost-graph done"

#build RandGraph

cd RandGraph
mkdir build 
cd build
cmake ..
make -j
cp randGraph ../../
cd ../../
echo "build RandGraph done"

#build my sssp 
cd ../
make -j
cp sssp testing/
cd testing
echo "build sssp done"
echo "Now start testing for 100 times"




errors=0
for i in {1..100}
do
./randGraph

#call boost
boost_result=$(./boost_graph < "tempMatrix.data")

#call sssp
sssp_result=$(./sssp < "tempMatrix.data")

echo "boost output: ""$boost_result" "vs" "sssp output:" "$sssp_result"
 
# compare
if [ "$boost_result" != "$sssp_result" ];then
echo "Found a different value " "boost:""$boost_result"  "sssp" "$sssp_result"

errorfile="$errors.errorlog"
errors=$[errors+1]
cp tempMatrix.data $errorfile

fi

#first generate
done





echo "Now compare time"
echo "boost"
time for i in {1..100}
do
./boost_graph < "tempMatrix.data"
done
echo "boost time"

echo "sssp"
time for i in {1..100}
do
./sssp < "tempMatrix.data"
done
echo "sssp time"
