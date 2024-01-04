default: TestRun

Preapre:
	ln -sf $(ProjectBase)/Samples .

TestRun: Execute
	./Execute --Input LEP1MC1994_recons_aftercut-014.root \
		--JetInput LEP1MC1994_recons_aftercut-014-Recluster.root \
		--Output LEP1MC1994_recons_aftercut-014-Combined.root \
		--InputTree RecoR4 \
		--OutputTree akR4ESchemeReclusteredJetTree \
		--Particle t

Execute: PrepareTree.cpp
	CompileRootMacro PrepareTree.cpp Execute
