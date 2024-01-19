

default: TestRun

TestRun: Execute
	./Execute --Input Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
		--Output LEP1MC1994_recons_aftercut-001_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00

FullRun: Execute
	bash MatchAll.sh

Execute: MatchTracks.cpp
	g++ MatchTracks.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

