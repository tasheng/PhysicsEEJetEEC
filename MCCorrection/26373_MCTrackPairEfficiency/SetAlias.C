{
   MatchedTree->SetAlias("GenTheta", "atan(GenZ/GenP)");
   MatchedTree->SetAlias("GenP", "sqrt(GenX*GenX+GenY*GenY+GenZ*GenZ)");

   PairTree->SetAlias("GenPT1", "sqrt(GenX1*GenX1+GenY1*GenY1)");
   PairTree->SetAlias("GenPT2", "sqrt(GenX2*GenX2+GenY2*GenY2)");
   PairTree->SetAlias("GenP1", "sqrt(GenX1*GenX1+GenY1*GenY1+GenZ1*GenZ1)");
   PairTree->SetAlias("GenP2", "sqrt(GenX2*GenX2+GenY2*GenY2+GenZ2*GenZ2)");
   PairTree->SetAlias("GenTheta1", "atan(GenZ1/GenP1)");
   PairTree->SetAlias("GenTheta2", "atan(GenZ2/GenP2)");
   PairTree->SetAlias("GenTheta", "((GenTheta1+GenTheta2)/2)");

   TLatex Latex;
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetNDC();
}
