#!/bin/bash
#compiles all progs

compiler="g++"
arguments="-O3 -frounding-math -mfpmath=387 -finline-functions -lgmp -fopenmp -std=c++11"
#codePath="/project/fh2-project-atpc18/ufpoa/code"
#runPath="/project/fh2-project-atpc18/ufpoa/run"
codePath="../code"
runPath="../run"

eval $compiler $codePath/sievesimp1.cpp  -o $runPath/sievesimp1 $arguments &

eval $compiler $codePath/sievesimp1B.cpp -o $runPath/sievesimp1B $arguments &

eval $compiler $codePath/sievesimp2.cpp  -o $runPath/sievesimp2 $arguments &

eval $compiler $codePath/sievesimp2B.cpp -o $runPath/sievesimp2B $arguments &

eval $compiler $codePath/sievesimp3.cpp  -o $runPath/sievesimp3 $arguments &

eval $compiler $codePath/sievesimp4.cpp  -o $runPath/sievesimp4 $arguments &

eval $compiler $codePath/sievesimp5.cpp  -o $runPath/sievesimp5 $arguments &

eval $compiler $codePath/sievesimp6.cpp  -o $runPath/sievesimp6 $arguments &

eval $compiler $codePath/sievesimp6B.cpp -o $runPath/sievesimp6B $arguments 

echo "Job Done"
