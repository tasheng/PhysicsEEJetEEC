inputFile=unfoldingE2C_DataUnfolding_DoubleLogBinning_BinningOption3_01032025.root
echo "Running Plotting scripts for input file: $inputFile"
echo "-------> Starting to plot unfolding stability"
./ExecuteUnfoldingStab --Input $inputFile \
        --Output StabilityCheckDataUnfolding_BinningOption3_RatioToMC \
        --Label "Iteration 1","Iteration 2","Iteration 3","Iteration 4","Iteration 5","Archived MC"\
        --Prefix "01092025" \
        --DoRatio true \
        --DoWeight false

echo "-------> Starting to plot unfolding closure"
./ExecuteClosure --Input $inputFile \
        --Output ClosureCheckDataUnfolding_BinningOption3_ \
        --Label "MC Gen 1D","Iteration 3"\
        --Prefix "01092025" \
        --DoRatio true \
        --Iter 3 \
        --DoWeight false