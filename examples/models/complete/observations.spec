// First run
#LineageOne[0] |=  $OverExpression;
#LineageOne[0] |=  $Initial1;
#LineageOne[20] |=  $Final1;
fixpoint(#LineageOne[20]);

// Second run
#LineageTwo[0] |=  $OverExpression;
#LineageTwo[0] |=  $Initial2;
#LineageTwo[20] |=  $Final2;
fixpoint(#LineageTwo[20]);

// Perturbations

$OverExpression :=
{
 FE(S1)=1 and 
 FE(S2)=0
};

// Gene expression patterns

$Initial1:=
{
 S1 = 1 and
 S2 = 0 and
 G2 = 1
};

$Initial2:=
{
 S1 = 1 and
 S2 = 0
};

$Final1:=
{
 G1 = 1
};

$Final2:=
{
 G1 = 1 and
 G2 = 1
};
