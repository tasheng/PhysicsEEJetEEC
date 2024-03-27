

TestRun: Execute
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root \
		--Output Output/TestLEP1MC_1994_014_Reco.root
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root \
		--Output Output/TestLEP1MC_1994_014_Gen.root --GenLevel true --Tree tgen

RunAll: Execute
	bash RunAll.sh

RunMC: Execute
	./Execute --Input $(ProjectBase)/Samples/HERWIG7/26436_AustinPilot/LEP-Matchbox-S1000-1_0_0.root \
		--Output Output/Herwig.root --GenLevel true
	./Execute --Input $(ProjectBase)/Samples/SHERPA/26436_AustinPilot/Sherpa_RNG100_0_0.root \
		--Output Output/Sherpa.root --GenLevel true
	./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240206_Pilot_Hannah/LEP1_pythia8_MC_withSphericity.root \
		--Output Output/Pythia8.root --GenLevel true --Tree tgen
	./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240211_Vincia_Hannah/LEP1_pythia8_MC_VINCIA.root\
		--Output Output/Pythia8Vincia.root --GenLevel true --Tree tgen
	./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240211_Dire_Hannah/LEP1_pythia8_MC_DIRE.root \
		--Output Output/Pythia8Dire.root --GenLevel true --Tree tgen

Execute: ReduceTree.cpp
	g++ ReduceTree.cpp -o Execute \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o \
		`root-config --cflags --libs`
