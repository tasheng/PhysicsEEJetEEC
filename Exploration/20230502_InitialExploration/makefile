JECBase     = $(ProjectBase)/CommonCode/text/
MCJEC       = $(JECBase)/JEC_EEAK4_MC_20200625.txt
ResidualJEC = $(JECBase)/JEC_EEAK4_DataL2_20200625.txt,$(JECBase)/JEC_EEAK4_DataL3_20200625.txt
DataJEC     = $(MCJEC),$(ResidualJEC)

default: Prepare TestRun

Prepare:
	ln -sf $(ProjectBase)/CommonCode/source/Dictionary_rdict.pcm .

TestRun: Execute
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root \
		--Output PlotReco.root \
		--Particle t --IsReco true \
		--Jet akR4ESchemeJetTree --JetR 0.4 --JEC $(MCJEC)
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root \
		--Output PlotGen.root \
		--Particle tgen --IsReco false \
		--Jet akR4ESchemeGenJetTree --JetR 0.4 --JEC ""
	./Execute --Input $(ProjectBase)/Samples/ALEPH/LEP1Data1994P2_recons_aftercut-MERGED.root \
		--Output PlotData.root \
		--Particle t --IsReco true \
		--Jet akR4ESchemeJetTree --JetR 0.4 --JEC $(DataJEC)

Execute: FirstExploration.cpp
	g++ FirstExploration.cpp -o Execute \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
