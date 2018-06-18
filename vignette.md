This is an example of using the code for model expansion and synthesis, on a dummy model, called "complete". Related files, called "model.net" (the non-expanded model file), "observations.spec" (the experimental data file) and "complete_crms.R" (file which contains the TF bindings), are located at examples/models/complete. Associated script is vignette.sh (start the script with the following command: "./vignette.sh". To clean the resulting files, type "./vignette.sh --clean").

**Note:** There should be a **white space** at the end of the 6th line (line for node declaration) of "model.net".

**Note:** For perturbation patterns (in observations.spec), the word "KnockDown" should appear in the name of the KO perturbation (ex. $KnockDownNanog), resp. "OverExpression" in the FE perturbation (ex. $OverExpressionOct4). It is advised to divide the observations into initial conditions, perturbations, and intermediary/final states.

### Getting started
Open a terminal in the root folder.

### Model expansion
R
```R
setwd("R/") ; source("expansion.R") ; setwd("../examples/models/complete/") ; source("complete_crms.R")

// Convert the non-expanded model into a file readable by the Python program
model_regular <- readModel("model.net") ; writeModel(model_regular, title="model_expanded.net", format="expanded")

// Expand the model file
model_expanded <- expansion(crms, "model.net", title="complete-full+_expanded.net", type="full", method="positive", inferTF=F, format="expanded") ; q(save="no")

### Get the graph associated to the regular model
(The image file is automatically created in the Python folder)
```bash
cd Python/ ; python solve.py launch complete igraph
```

### Get the graph associated to the expanded model
(The image file is automatically created in the Python folder)
```bash
python solve.py launch complete igraph --model "complete-full+_expanded"
```

### Model synthesis for regular model
Here all solutions have the exact same topology (because all interactions are definite in the abstract model)
The following command stores the output (giving the values of selected interaction vector and gene regulatory function types for each node) in a text file
Solution GRFs are written in files in Python/results/ (one by solution model)
```bash
python solve.py run --simplify complete/model_expanded.net complete/observations.spec > results/report-complete-model_expanded.txt
```

### Model synthesis for expanded model
```bash
python solve.py run --simplify complete/complete-full+_expanded.net complete/observations.spec > results/report-complete-full+_expanded.txt
```

### Generate (at most) 2 trajectories of length 5 for the first solution found for the expanded model
Starting in state Initial2 (defined in observations.spec file located in the same folder as the model file) (with no fixpoint condition on the 5th state reached) and look for condition Final1 (defined in observations.spec)
```bash
python solve.py launch complete --model complete-full+_expanded --q0 Initial2 --nstep 5 --solmax 2 --modelID 1 --steadyStates 0 --expnames Final1 > results/trajectories-model1-5steps-Initial2-Final1-noSS.txt ; cd ../
```
