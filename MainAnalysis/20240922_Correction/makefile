

default: TestRun

TestRun: Execute
# 	./Execute --InputData /data/hbossi/PhysicsEEJetEEC/Unfolding/20240922_UnfoldingThetaZ/unfoldingE2C_DoubleLogBinning_FrozenMCBinning_09232024.root
	./Execute --InputData /data/hbossi/PhysicsEEJetEEC/Unfolding/20240922_UnfoldingThetaZ/unfoldingE2C_DataUnfolding_DoubleLogBinning_09232024.root

Execute: Correction.cpp
	g++ Correction.cpp -o Execute `root-config --cflags --libs` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
