CXX:=icpc

TARGET = sievesimp

LDFLAGS=-g -O3 -frounding-math -mfpmath=387 -finline-functions -qopenmp -std=c++11
LDLIBS=-lgmp 
sievesimp: sievesimp.cpp diophappr.h solvemodz.h
	$(CXX) $(LDFLAGS) sievesimp.cpp -o sievesimp $(LDLIBS)
#sievesimp: sievesimp.o
#	$(CXX) $(LDFLAGS) -o sievesimp $(LDLIBS)
#sievesimp.o: sievesimp.cpp diophappr.h solvemodz.o
#	$(CXX) $(LDFLAGS) -c sievesimp.cpp diophappr.h solvmodz.h $(LDLIBS)
#solvemodz.o: solvemodz.h
#	$(CXX) $(LDFLAGS) -c solvemodz.h $(LDLIBS)


