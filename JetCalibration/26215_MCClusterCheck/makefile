

default: TestRun

TestRun: Execute
	# ./Execute --Input Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root --Output Meow.root \
	# 	--Jet akR4ESchemeJetTree --P t --JetR 0.4
	# ./Execute --Input Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root --Output Meow.root \
	# 	--Jet akR4ESchemeGenJetTree --P tgen --JetR 0.4
	./Execute --Input Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root --Output Meow.root \
		--Jet akR4ESchemeGenJetBeforeTree --P tgenBefore --JetR 0.4

Execute: JetCluster.cpp
	g++ JetCluster.cpp -o Execute \
		`root-config --libs --cflags` \
		`$(FASTJET_BASE)/bin/fastjet-config --libs --cxxflags` \
		-I$(ProjectBase)/CommonCode/include $(ProjectBase)/CommonCode/library/*.o

