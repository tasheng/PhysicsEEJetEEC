default: TestPlot

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

RunAll: Execute
	bash RunAll.sh

TestPlot: ExecutePlot
	./ExecutePlot

Execute: FirstExploration.cpp
	g++ FirstExploration.cpp -o Execute \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExecutePlot: MakePlot.cpp
	g++ MakePlot.cpp -o ExecutePlot \
		-I$(ProjectBase)/CommonCode/include \
		`root-config --libs --cflags`
