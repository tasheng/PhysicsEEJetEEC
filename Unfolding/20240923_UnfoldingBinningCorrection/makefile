default: TestRun

TestRun: Execute
	./Execute --Input /data/hbossi/PhysicsEEJetEEC/Unfolding/20240922_UnfoldingThetaZ/unfoldingE2C_DoubleLogBinning_FrozenBinning_v2_09232024.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_z.root \
			     --UnfoldingBinCorrArgName z --MakeUnfoldingBinCorrFactor true
	./Execute --Input /data/hbossi/PhysicsEEJetEEC/Unfolding/20240922_UnfoldingThetaZ/unfoldingE2C_DoubleLogBinning_FrozenBinning_v2_09232024.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_z.root \
			     --UnfoldingBinCorrArgName z
	./Execute --Input /data/hbossi/PhysicsEEJetEEC/Unfolding/20240922_UnfoldingThetaZ/unfoldingE2C_DoubleLogBinning_FrozenBinning_v2_09232024.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_theta.root \
			     --UnfoldingBinCorrArgName theta --MakeUnfoldingBinCorrFactor true
	./Execute --Input /data/hbossi/PhysicsEEJetEEC/Unfolding/20240922_UnfoldingThetaZ/unfoldingE2C_DoubleLogBinning_FrozenBinning_v2_09232024.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_theta.root \
			     --UnfoldingBinCorrArgName theta

Execute: unfoldingBinningCorrection.cpp
	g++ unfoldingBinningCorrection.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o