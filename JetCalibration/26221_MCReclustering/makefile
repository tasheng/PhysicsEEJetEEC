

default: TestRun

TestRun: Execute
	./Execute --Input Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root --Output Samples/ALEPHMCRecluster/Jet-014.root \
		--Reco t --Gen tgen --GenBefore --JetR 0.2,0.4,0.6,0.8,1.0

Execute: JetCluster.cpp
	g++ JetCluster.cpp -o Execute \
		`root-config --libs --cflags` \
		`$(FASTJET_BASE)/bin/fastjet-config --libs --cxxflags` \
		-I$(ProjectBase)/CommonCode/include $(ProjectBase)/CommonCode/library/*.o

