CXX:=g++

TARGET = simplesieve

LDFLAGS=-g -O3 -frounding-math -mfpmath=387 -finline-functions  -std=c++11
LDLIBS=-lgmp
simplesieve:   diophappr.h solvemodz.h
	$(CXX) $(LDFLAGS) simplesieve.cpp -o simplesieve $(LDLIBS)
sievesimp: sievesimp.o
	$(CXX) $(LDFLAGS) -o sievesimp $(LDLIBS)
sievesimp.o: sievesimp.cpp diophappr.h solvemodz.o
	$(CXX) $(LDFLAGS) -c sievesimp.cpp diophappr.h solvmodz.h $(LDLIBS)
solvemodz.o: solvemodz.h
	$(CXX) $(LDFLAGS) -c solvemodz.h $(LDLIBS)
