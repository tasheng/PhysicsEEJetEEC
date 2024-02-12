default: FullPlot 

Prepare:
	ln -sf $(ProjectBase)/CommonCode/source/Dictionary_rdict.pcm .

TestRun: Execute
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root \
		--Output TestReco.root \
		--Particle t --IsReco true
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root \
		--Output TestGen.root \
		--Particle tgen --IsReco false
	./Execute --Input $(ProjectBase)/Samples/ALEPH/LEP1Data1994P2_recons_aftercut-MERGED.root \
		--Output TestData.root \
		--Particle t --IsReco true

RunPythia8: Execute
	# ./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240206_Pilot_Hannah/LEP1_pythia8_MC.root \
	# 	--Output PlotPythia8.root \
	# 	--Particle tgen --IsReco false --DoEENormalize true --DoWeight false
	./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240206_Pilot_Hannah/LEP1_pythia8_MC_withSphericity.root \
		--Output PlotPythia8Sphericity.root \
		--Particle tgen --IsReco false --DoEENormalize true --DoWeight false --CheckSphericity true
	./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240206_Pilot_Hannah/LEP1_pythia8_MC_withSphericity_v2.root \
		--Output PlotPythia8SphericityV2.root \
		--Particle tgen --IsReco false --DoEENormalize true --DoWeight false --CheckSphericity true

RunPythia8Vincia: Execute
	./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240211_Vincia_Hannah/LEP1_pythia8_MC_VINCIA.root \
		--Output PlotPythia8Vincia.root \
		--Particle tgen --IsReco false --DoEENormalize true --DoWeight false --CheckSphericity true

RunPythia8Dire: Execute
	./Execute --Input $(ProjectBase)/Samples/PYTHIA8/20240211_Dire_Hannah/LEP1_pythia8_MC_DIRE.root \
		--Output PlotPythia8Dire.root \
		--Particle tgen --IsReco false --DoEENormalize true --DoWeight false --CheckSphericity true

RunSherpa: Execute
	./Execute --Input $(ProjectBase)/Samples/SHERPA/26436_AustinPilot/Sherpa_RNG100_0_0.root \
		--Output PlotSherpa.root \
		--Particle t --IsReco false --DoEENormalize true --DoWeight false --DoCheck true

RunHerwig: Execute
	./Execute --Input $(ProjectBase)/Samples/HERWIG7/26436_AustinPilot/LEP-Matchbox-S1000-1_0_0.root \
		--Output PlotHerwig.root \
		--Particle t --IsReco false --DoEENormalize true --DoWeight false --DoCheck true

RunHerwigPartial: Execute
	./Execute --Input $(ProjectBase)/Samples/HERWIG7/26436_AustinPilot/LEP-Matchbox-S1000-1_0_0.root \
		--Output PlotHerwigPartial.root --Fraction 0.05 \
		--Particle t --IsReco false --DoEENormalize true --DoWeight false

RunAll: Execute
	bash RunAll.sh

RunGenReject3Jet: Execute
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
		--Output PlotGenNo3Jet_01.root \
      --Particle tgen --IsReco false --DoEENormalize true --DoWeight false \
		--Reject3Jet true --Jet akR4ESchemeGenJetTree

TestPlot: ExecutePlot
	./ExecutePlot

FullPlot: ExecutePlot
	./ExecutePlot --Input PlotGenSTheta.root,PlotReco.root \
		--Label "Archived MC","Archived MC + detector" \
		--Prefix "GenReco"
	./ExecutePlot --Input PlotGenSTheta.root,PlotPythia8.root,PlotPythia8Sphericity.root,PlotPythia8SphericityV2.root \
		--Label "Archived MC","Pythia8 + no cut","Pythia8 + sphericity cut 1","Pythia8 + sphericity cut 2" \
		--Prefix "Pythia"
	./ExecutePlot --Input PlotGenSTheta.root,PlotPythia8Sphericity.root,PlotHerwig.root,PlotSherpa.root,PlotPythia8Vincia.root,PlotPythia8Dire.root \
		--Label "Archived MC","Pythia8","Herwig","Sherpa","Pythia8 Vincia","Pythia8 Dire" \
		--Prefix "MC"
	./ExecutePlot --Input PlotReco.root,PlotData.root \
		--Label "Archived MC + detector","Data" \
		--Prefix "RecoData"
	./ExecutePlot --Input PlotGenSTheta.root,PlotGen.root \
		--Label "Archived MC","Archived MC (no sphericity flag)" \
		--Prefix "GenSphericity"
	./ExecutePlot --Input PlotData.root,PlotDataLEP2Cut.root \
		--Label "Data","Data with LEP2 Cut" \
		--Prefix "LEP1LEP2"
	./ExecutePlot --Input PlotGenSTheta.root,PlotGenNo3Jet_01.root \
		--Label "Archived MC","Archived MC; E_{jet 3} < 5 GeV" \
		--Prefix "ThirdJet"
	./ExecutePlot --Input TestGen.root \
		--Label "Example" \
		--Prefix "Example" --DoRatio false

StatPlot: ExecuteStat
	./ExecuteStat

Convert:
	bash ConvertToPNG.sh

Execute: FirstExploration.cpp
	g++ FirstExploration.cpp -o Execute \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExecutePlot: MakePlot.cpp
	g++ MakePlot.cpp -o ExecutePlot \
		-I$(ProjectBase)/CommonCode/include \
		`root-config --libs --cflags` -g

ExecuteStat: MakeStatisticsPlot.cpp
	g++ MakeStatisticsPlot.cpp -o ExecuteStat \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o


