#!/bin/bash
doTest=0
echo "First round, single threaded:"
eval './sievesimp1 5000000000000000000 40000000 1 $doTest'
eval './sievesimp2 5000000000000000000 40000000 1 $doTest'
eval './sievesimp3 5000000000000000000 40000000 1 $doTest'
eval './sievesimp4 5000000000000000000 40000000 1 $doTest'
eval './sievesimp6 5000000000000000000 40000000 1 $doTest'
echo "=========================================="

echo ""
echo "Second round, two threads:"
eval './sievesimp3 5000000000000000000 40000000 2 $doTest'
eval './sievesimp4 5000000000000000000 40000000 2 $doTest'
eval './sievesimp6 5000000000000000000 40000000 2 $doTest'
echo "=========================================="

echo ""
echo "Third round, four threads:"
eval './sievesimp3 5000000000000000000 40000000 4 $doTest'
eval './sievesimp4 5000000000000000000 40000000 4 $doTest'
eval './sievesimp6 5000000000000000000 40000000 4 $doTest'
echo "=========================================="
