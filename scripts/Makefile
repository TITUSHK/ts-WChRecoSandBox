CXX = `root-config --cxx`


CXXFLAGS  = -g -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -O3
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lCore  -lRGL -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer -lTMVA
INCPATH   = -I./../include/

all : 	objects

clean :
	rm -f *.o gen

objects : 	SandFit.o  

SandFit.o : SandFit.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATt) -c $< -o $@

SandFit : SandFit.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -L./../lib -lWCLAnalysis
