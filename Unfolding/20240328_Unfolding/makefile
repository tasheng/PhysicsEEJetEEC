default: TestRun

TestRun: Execute
	./Execute --Input Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
		--Output LEP1MC1994_recons_aftercut-001_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00

FullRun: Execute
	bash MatchSome.sh

Execute: MatchEEC.cpp
	g++ MatchEEC.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

