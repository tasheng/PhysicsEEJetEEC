default: TestRun

TestRun: Execute
	./Execute --Input ../../Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
			     --EvtSelEffName EvtSelEff.root --MakeEvtSelEffCorrFactor true
	./Execute --Input ../../Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
			     --EvtSelEffName EvtSelEff.root

Execute: evtSelEffCorr.cpp
	g++ evtSelEffCorr.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o