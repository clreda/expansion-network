#!/bin/bash

if [ "$1" == "--clean" ];
then
echo "MSG: Cleaning"
cd Python/
rm -rf *.png results/ *.pyc *.dot *.tmp *.txt ; mkdir results/
cd ../examples/models/complete/
rm -f model_expanded.net complete-full+_expanded.net
cd ../../
rm -f *.pyc
cd ../
else
echo "MSG: Expansion procedure"
R -e 'setwd("R/"); source("expansion.R"); setwd("../examples/models/complete/"); source("complete_crms.R"); model_regular <- readModel("model.net"); writeModel(model_regular, title="model_expanded.net", format="expanded"); model_expanded <- expansion(crms, "model.net", title="complete-full+_expanded.net", type="full", method="positive", inferTF=F, format="expanded")'
echo "MSG: Graph for the regular model"
cd Python/
python solve.py launch complete igraph
echo "MSG: Graph for the expanded model"
python solve.py launch complete igraph --model "complete-full+_expanded"
echo "MSG: Model synthesis for the non-expanded model"
python solve.py run --simplify complete/model_expanded.net complete/observations.spec > results/report-complete-model_expanded.txt
echo "MSG: Model synthesis for the expanded model"
python solve.py run --simplify --visualize complete/complete-full+_expanded.net complete/observations.spec > results/report-complete-full+_expanded.txt
echo "MSG: Generate trajectories for the first solution model found (expanded model)"
python solve.py launch complete --model complete-full+_expanded --q0 Initial2 --nstep 5 --solmax 2 --modelID 1 --steadyStates 0 --expnames Final1 > results/trajectories-model1-5steps-Initial2-Final1-noSS.txt
cd ../
fi
echo "MSG: End of vignette script"
exit
