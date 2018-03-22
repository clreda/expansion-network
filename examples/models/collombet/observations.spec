// LMPP to CLP stable state
#LineageOne[0] |=  $LymphoidMyeloidPP;
#LineageOne[0] |=  $Il7Csf1NoGeneOverExpression;
#LineageOne[10] |=  $FinalStateCLP;
#LineageOne[11] |=  $FinalStateCLP;
#LineageOne[12] |=  $FinalStateCLP;
#LineageOne[20] |=  $FinalStateCLP;

// LMPP to GMP stable state
#LineageTwo[0] |=  $LymphoidMyeloidPP;
#LineageTwo[0] |=  $Il7Csf1NoGeneOverExpression;
#LineageTwo[5] |=  $PreGM;
#LineageTwo[10] |=  $Gmp;
#LineageTwo[19] |=  $FinalStateGMP;
#LineageTwo[20] |=  $FinalStateGMP;

// CLP and stimulation with Csf1 and Il7 to B Cell stable state
#LineageThree[0] |=  $FinalStateCLP;
#LineageThree[0] |=  $Il7Csf1GeneOverExpression;
#LineageThree[15] |=  $FinalStateBCell;
#LineageThree[16] |=  $FinalStateBCell;
#LineageThree[17] |=  $FinalStateBCell;
#LineageThree[20] |=  $FinalStateBCell;

// GMP and stimulation with Csf1 and Il7 to Mac stable state
#LineageFour[0] |= $FinalStateGMP;
#LineageFour[0] |=  $Il7Csf1GeneOverExpression;
#LineageFour[13] |= $FinalStateMac;
fixpoint(#LineageFour[13]);

// Perturbations

$Il7Csf1GeneOverExpression :=
{
 FE(Il7)=1 and 
 FE(Csf1)=1
};

$Il7Csf1NoGeneOverExpression :=
{
 FE(Il7)=0 and 
 FE(Csf1)=0
};

// Gene expression patterns

$LymphoidMyeloidPP:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 0 and
 Foxo1 = 0 and
 E2a = 1 and
 Il7ra = 0 and
 Ets1 = 0 and
 Mef2c = 0 and
 Flt3 = 1 and
 Ikzf1 = 1 and
 Gfi1 = 0 and
 Runx1 = 1 and
 Spi11 = 1 and
 Csf1r = 0 and
 Cebpa = 0 and
 Cebpb = 0 and
 Egr2 = 0 and
 Id2 = 0 and
 Mac1 = 0
};

$FinalStateCLP:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 0 and
 Foxo1 = 1 and
 E2a = 1 and
 Il7ra = 1 and
 Ets1 = 1 and
 Mef2c = 1 and
 Flt3 = 1 and
 Ikzf1 = 1 and
 Gfi1 = 1 and
 Runx1 = 1 and
 Spi11 = 0 and
 Spi12 = 0 and
 Csf1r = 1 and
 Cebpa = 0 and
 Cebpb = 0 and
 Egr2 = 0 and
 Id2 = 0 and
 Mac1 = 0
};

$Clp:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 1 and
 Foxo1 = 0 and
 E2a = 1 and
 Il7ra = 1 and
 Ets1 = 0 and
 Mef2c = 0 and
 Flt3 = 1 and
 Ikzf1 = 1 and
 Gfi1 = 0 and
 Runx1 = 1 and
 Spi11 = 0 and
 Spi12 = 0 and
 Csf1r = 0 and
 Cebpa = 0 and
 Cebpb = 0 and
 Egr2 = 0 and
 Id2 = 0 and
 Mac1 = 0
};

$ProB:=
{
 Cd19 = 1 and
 Pax5 = 1 and
 Ebf1 = 1 and
 Foxo1 = 1 and
 E2a = 0 and
 Il7ra = 1 and
 Ets1 = 1 and
 Mef2c = 1 and
 Flt3 = 0 and
 Ikzf1 = 1 and
 Gfi1 = 1 and
 Runx1 = 0 and
 Spi11 = 0 and
 Spi12 = 0 and
 Csf1r = 0 and
 Cebpa = 0 and
 Cebpb = 1 and
 Egr2 = 1 and
 Id2 = 1 and
 Mac1 = 0
};

$PreB:=
{
 Cd19 = 1 and
 Pax5 = 1 and
 Ebf1 = 1 and
 Foxo1 = 1 and
 E2a = 0 and
 Il7ra = 1 and
 Ets1 = 1 and
 Mef2c = 1 and
 Flt3 = 0 and
 Ikzf1 = 1 and
 Gfi1 = 0 and
 Runx1 = 0 and
 Spi11 = 0 and
 Spi12 = 0 and
 Csf1r = 0 and
 Cebpa = 0 and
 Cebpb = 1 and
 Egr2 = 1 and
 Id2 = 0 and
 Mac1 = 0
};

$PreGM:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 0 and
 Foxo1 = 0 and
 E2a = 1 and
 Il7ra = 0 and
 Ets1 = 0 and
 Mef2c = 0 and
 Flt3 = 1 and
 Ikzf1 = 1 and
 Gfi1 = 0 and
 Runx1 = 1 and
 Spi11 = 1 and
 Csf1r = 0 and
 Cebpa = 1 and
 Cebpb = 0 and
 Egr2 = 0 and
 Id2 = 0 and
 Mac1 = 0
};

$Gmp:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 0 and
 Foxo1 = 0 and
 E2a = 1 and
 Il7ra = 0 and
 Ets1 = 0 and
 Mef2c = 0 and
 Flt3 = 0 and
 Ikzf1 = 1 and
 Gfi1 = 1 and
 Runx1 = 1 and
 Spi11 = 0 and
 Spi12 = 0 and
 Csf1r = 0 and
 Cebpa = 1 and
 Cebpb = 0 and
 Egr2 = 0 and
 Id2 = 0 and
 Mac1 = 0
};

$Mac:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 0 and
 Foxo1 = 0 and
 E2a = 0 and
 Il7ra = 0 and
 Ets1 = 0 and
 Mef2c = 0 and
 Flt3 = 0 and
 Ikzf1 = 0 and
 Gfi1 = 0 and
 Runx1 = 1 and
 Spi12 = 1 and
 Csf1r = 1 and
 Cebpa = 0 and
 Cebpb = 1 and
 Egr2 = 0 and
 Id2 = 1 and
 Mac1 = 1
};

$FinalStateBCell:=
{
 Cd19 = 1 and
 Pax5 = 1 and
 Ebf1 = 1 and
 Foxo1 = 1 and
 E2a = 1 and
 Il7ra = 1 and
 Ets1 = 1 and
 Mef2c = 1 and
 Flt3 = 1 and
 Ikzf1 = 1 and
 Gfi1 = 1 and
 Runx1 = 1 and
 Spi11 = 0 and
 Spi12 = 0 and
 Csf1r = 0 and
 Cebpa = 0 and
 Cebpb = 0 and
 Egr2 = 1 and
 Id2 = 0 and
 Mac1 = 0
};

$FinalStateGMP:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 0 and
 Foxo1 = 0 and
 E2a = 0 and
 Il7ra = 0 and
 Ets1 = 0 and
 Mef2c = 0 and
 Flt3 = 0 and
 Ikzf1 = 0 and
 Gfi1 = 1 and
 Runx1 = 1 and
 Spi12 = 0 and
 Csf1r = 1 and
 Cebpa = 1 and
 Cebpb = 0 and
 Egr2 = 0 and
 Id2 = 0 and
 Mac1 = 0
};

$FinalStateMac:=
{
 Cd19 = 0 and
 Pax5 = 0 and
 Ebf1 = 0 and
 Foxo1 = 0 and
 E2a = 0 and
 Il7ra = 0 and
 Ets1 = 0 and
 Mef2c = 0 and
 Flt3 = 0 and
 Ikzf1 = 0 and
 Gfi1 = 0 and
 Runx1 = 1 and
 Spi12 = 1 and
 Csf1r = 1 and
 Cebpa = 1 and
 Cebpb = 1 and
 Egr2 = 1 and
 Id2 = 1 and
 Mac1 = 1
};

// From SI (Collombet et al 2017)

$SignatureMac:=
{
 Mac1 = 1
};

$SignaturePreProB:=
{
 E2a = 1 and
 Foxo1 = 1 and
 Ebf1 = 1 and
 Cd19 = 0
};

$SignatureBCell:=
{
 Cd19 = 1
};

$SignatureGMP:=
{
 Spi11 = 1 and
 Cebpa = 1 and
 Gfi1 = 1 and
 Cebpb = 0
};

$SignatureCLP:=
{
 E2a = 1 and
 Foxo1 = 1 and
 Cebpa = 0 and
 Ebf1 = 0
};

$SignatureMP:=
{
 Spi11 = 1 and
 Ikzf1 = 1 and
 Flt3 = 1 and
 Foxo1 = 0
};

