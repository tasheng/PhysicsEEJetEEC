{
   RecoR4->SetAlias("jtp", "(jtpt*cosh(jteta))");
   RecoR4->SetAlias("jtE", "sqrt(jtp*jtp+jtm*jtm)");
   RecoR4->SetAlias("jtpz", "(jtpt*sinh(jteta))");
   RecoR4->SetAlias("jttheta", "acos(jtpz/jtp)");
   RecoR4->SetAlias("jtthetagap", "(3.14159/2-abs(jttheta-3.14159/2))");
   
   RecoR10->SetAlias("jtp", "(jtpt*cosh(jteta))");
   RecoR10->SetAlias("jtE", "sqrt(jtp*jtp+jtm*jtm)");
   RecoR10->SetAlias("jtpz", "(jtpt*sinh(jteta))");
   RecoR10->SetAlias("jttheta", "acos(jtpz/jtp)");
   RecoR10->SetAlias("jtthetagap", "(3.14159/2-abs(jttheta-3.14159/2))");
   
   GenR10->SetAlias("jtp", "(jtpt*cosh(jteta))");
   GenR10->SetAlias("jtE", "sqrt(jtp*jtp+jtm*jtm)");
   GenR10->SetAlias("jtpz", "(jtpt*sinh(jteta))");
   GenR10->SetAlias("jttheta", "acos(jtpz/jtp)");
   GenR10->SetAlias("jtthetagap", "(3.14159/2-abs(jttheta-3.14159/2))");
   
   GenBeforeR10->SetAlias("jtp", "(jtpt*cosh(jteta))");
   GenBeforeR10->SetAlias("jtE", "sqrt(jtp*jtp+jtm*jtm)");
   GenBeforeR10->SetAlias("jtpz", "(jtpt*sinh(jteta))");
   GenBeforeR10->SetAlias("jttheta", "acos(jtpz/jtp)");
   GenBeforeR10->SetAlias("jtthetagap", "(3.14159/2-abs(jttheta-3.14159/2))");
}
