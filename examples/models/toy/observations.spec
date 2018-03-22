// Observation predicates
$Conditions1 := { S1 = 0 and S2 = 1};
$Conditions2 := { S1 = 1 and S2 = 1};
$Expression1 := {A = 1 and B = 1 and C = 1};
$Expression2 := {A = 0 and B = 1 and C = 1};

// Observations
#Experiment1[0] |= $Conditions1 and
#Experiment1[0] |= $Expression1 and
#Experiment1[18] |= $Expression2 and
fixpoint(#Experiment1[18]);

#Experiment2[0] |= $Conditions2 and
#Experiment2[0] |= $Expression2 and
#Experiment2[18] |= $Expression1 and
fixpoint(#Experiment2[18]);