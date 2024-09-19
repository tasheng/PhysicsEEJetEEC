	# ./ExecuteHistFillerGen --Input "$ProjectBase/Samples/MCSamples/LEP1_pythia8_MC_VINCIA.root" --Output "$ProjectBase/Samples/MCSamples/PYTHIA8_VINCIA.root" \
	# 	--Gen tgen --Fraction 1.0


    # ./ExecuteHistFillerGen --Input "$ProjectBase/Samples/MCSamples/Sherpa_RNG100_0_0.root" --Output "$ProjectBase/Samples/MCSamples/Sherpa.root" \
    #     --Gen t --Fraction 1.0 --IsSherpa true

    ./ExecuteHistFillerGen --Input "$ProjectBase/Samples/MCSamples/LEP1_HERWIG.root" --Output "$ProjectBase/Samples/MCSamples/HERWIG.root" \
        --Gen t --Fraction 1.0 --IsSherpa true