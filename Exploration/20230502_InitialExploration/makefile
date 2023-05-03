JECBase     = $(ProjectBase)/CommonCode/text/
MCJEC       = $(JECBase)/JEC_EEAK4_MC_20200625.txt
ResidualJEC = $(JECBase)/JEC_EEAK4_DataL2_20200625.txt,$(JECBase)/JEC_EEAK4_DataL3_20200625.txt
DataJEC     = $(MCJEC),$(ResidualJEC)

default: Prepare TestRun

Prepare:
	ln -sf $(ProjectBase)/CommonCode/source/Dictionary_rdict.pcm .

TestRun: Execute
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root \
		--Jet akR4ESchemeJetTree --JetR 0.4 \
		--JEC $(DataJEC)

Execute: FirstExploration.cpp
	g++ FirstExploration.cpp -o Execute \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
