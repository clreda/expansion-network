## Automated Inference of Gene Regulatory Networks Using Explicit Regulatory Modules

This GitHub project provides the code for the network inference and model expansion procedures, along with all the related scripts used in the paper *Automated Inference of Gene Regulatory Networks Using Explicit Regulatory Modules* (non published).

The folder "docs" contains supplementary material for the tests performed in the core paper.

## Requirements & Installation

**Python:** Packages **Z3** (SMT solver) and **igraph** (GRN visualization).

For Debian Linux:

`sudo apt install python-pip`

`sudo python -m pip install z3-solver`

`sudo python -m pip install python-igraph`


**R:** Packages **QCA** (for boolean function simplification), **XML** (parsing of SBML/XML files) and **igraph** (GRN visualization).

For Debian Linux:

`R`

`> install.packages("QCA")`

`> install.packages("XML")`

`> install.packages("igraph")`


## Description



- *examples* Example models and test files.

- *Python* Contains the code for the network inference procedure.

- *R* Contains code related to the model expansion procedure, boolean function simplification and decomposition, conversion from SBML model to IGRAPH object/REIN model, and regulatory module analysis.


## Network Inference Procedure

Method similar to **RE:IN** (**reference:** *Yordanov, B., Dunn, S. J., Kugler, H., Smith, A., Martello, G., & Emmott, S. (2016). A method to identify and analyze biological programs through automated reasoning. NPJ systems biology and applications, 2, 16010.*).

This code can be used for *regular* or *expanded* Boolean models, with some slight modifications in the format of the model file compared to **RE:IN**-formatted files:

In the line that describes the nodes, instead of:

`Name[perturbations](regulation functions)`

 
it should be read:



- `Name[perturbations]{}(regulation functions)` (if the node `Name` is NOT a regulatory module).

- `Name[perturbations]{gene}(regulation functions)` (otherwise: `gene` is the gene regulated by module `Name`).


Compared to **RE:IN**, two types of constraints have been added:



- If an interaction of type "TF to RM" is selected, then the interaction "RM -> gene" must be selected.

- If an interaction of type "RM to gene" is not selected, then no associated interaction of type "TF to RM" can be selected.


Results from the solver are stored in the folder "results/" located in the folder related to the solver, and models should be present in the subdirectory "models" of the folder "examples". This can be modified in file `global_paths.py`.

### Usage

#### Test files
To test functions from *filename* in {*shortcuts* | *utils* | *grn_inference* | *launch_model*}, type in the terminal (in the "Python" folder):

`python tests.py filename`


#### Network inference from an abstract model and an experiments file

`python solve.py run [--simplify] [--visualize] model experiments`


- Option *--simplify* uses the boolean reducer to simplify the expression of the resulting GRFs.

- Option *--visualize* build each graph image corresponding to each model candidate found.


"model" and "experiments" are the respective names of the model and the experiments files, that must contain the file extension.

**Example:** To test the **Toy** model from RE:IN (located at "../examples/models/toy"), and visualize all the GRNs corresponding to the model solutions, type the following command:

`python solve.py run --visualize toy/model_expanded.net toy/observations.spec`


#### Get trajectory from a candidate model and an initial state

By default, the model file and the experiments files are respectively named "model.net" and "observations.spec". This can be modified in file `global_paths.py`. They are stored in the models folder in a subdirectory called "model_name". 

`python solve.py launch model_name [igraph] [--model (default:model_expanded)] [--experiments (default:observations)] [--KO (default:"")] [--FE (default:"")] [--modelID (default:0)] [--q0 (default:111...11)] [--nstep (default:20)] [--solmax (default:10)] [--steadyStates (default:0)] [--expnames condition1 condition2 ...]`


- if *igraph* is present then it returns the igraph associated with model solution number *modelID* (or the full abstract model if no non-empty model solution list is provided, see `solve.py`)

- *model* is the model file name (without the extension ".net").

- *experiments* is the experiments file name (without the extension ".spec").

- *KO* is a set of knocked-out perturbations (as a condition name that appears in the experiments file, with "KnockDown" in it).

- *FE* is a set of forcibly-expressed perturbations (as a condition name that appears in the experiments file, with "OverExpression" in it).

- *modelID* is the index of the model candidate on which the function should be applied.

- *q0* is the initial state (a sequence of 0's and 1's of size #nodes or a condition name that appears in the experiments file).

- *nstep* is the length of the trajectories to generate.

- *solmax* is the maximum number of trajectories to generate.

- *steadyStates*, if equal to 1, adds a fix point constraint at step *nstep* to find steady states.

- *expnames* checks at each step of the trajectories if the conditions (which should appear in the experiments file) appear.


All the arguments between brackets are optional. The optional arguments following a double dash can be given in any order, **except for** argument list associated with "--expnames" that must be last.

**Example:** 

To use the regular **Collombet** model (located at "../examples/models/collombet/model_expanded.net"), and get the GRN corresponding to the abstract model (in PNG format), type the following command:

`python solve.py launch collombet igraph`


To use an expanded **Collombet** model (located at "../examples/models/collombet/expanded.net"), and get the GRN corresponding to the abstract model (in PNG format), type the following command:

`python solve.py launch collombet igraph --model expanded`


To use the regular **Collombet** model (located at "../examples/models/collombet/model_expanded.net"), and generate at most 5 trajectories of length 40 using the first model solution found by the solver, from initial state *LymphoidMyeloidPP* (available in the "observations.spec" file associated with the Collombet model), with the fix point condition applied at states of step 40, and look for condition *FinalStateMac*, type the following command:

`python solve.py launch collombet --q0 LymphoidMyeloidPP --nstep 40 --solmax 5 --steadyStates 1 --expnames FinalStateMac`


## Boolean Reducer

Uses Quine-McCluskey minimization. Reduces functions with less than 8 variables.

### Specifications

**Generic Boolean Formula OBJECT**


- Obtained from a text file (such as GRF output of RE:IN) with function "readModel" in *boolean_reducer.R*

- Can be used to generate a REIN-like formatted text file with function "writeModel" in *boolean_reducer.R*

Named nested list such as:
generics = list(conditions: ((string character) R vector) R named list (names = outcomes), outcomes: (string character) R list, formula: (boolean R function) R named list (names = outcomes))

### Usage

To minimize the functions (one per line: Gene' = ( AND | NOT | OR | ...)-) written in file "filename" located in the same folder as file "boolean_reducer.R":

`̀R -e "source('boolean_reducer.R'); minimizedFunctionsList <- boolean_reducer(filename)"` 

The previous command will write in the same folder a text file named "simplified--filename" the minimized functions.

## Expansion procedure

### Specifications

**Generic Model OBJECT**

- Obtained from a REIN-like formatted text file

- with function "readModel" in *io_expansion.R*

- Can be used to generate a REIN-like formatted text file

- with function "writeModel" in *io_expansion.R*

Named nested list such as:
model = list(directives: list(directive_updates: c("async", "sync"), directive_length: positive integer cast to string character, directive_uniqueness: c("interactions", "full", "paths"), directive_limit: nonnegative integer cast to string character, directive_regulation: c("threshold", "legacy", "default")), nodes: (string character) R named list (names = nodes), conditions: (positive integer list cast to single string character) R named list (names = nodes), perturbations: (c("", "-", "-", "--", "!")) R named list (names = nodes), definites: (3-sized string character R vector) R named list, optionals: (3-sized string character R vector) R named list)

e.g. model has 3 nodes A, B, C, with no perturbed genes, all conditions set to 
function template 0, and 2 positive definite interactions A -> B and A -> C
and 1 negative positive interaction B -> C
model = list(directives=list(directive_update="sync", directive_length="20", directive_uniqueness="interactions", directive_limit="0", directive_regulation="legacy"), nodes=list(A="A", B="B", C="C"), conditions=list(A="0", B="0", C="0"), perturbations=list(A="", B="", C=""), definites=list(definite1=c("A", "B", "positive"), definite2=c("A", "C", "positive")), optionals=list(optional1=c("B", "C", "negative")))

**Regulatory Module Set OBJECT**

Named nested list such as:
crms = list(nodes: (string character) R vector, isTF: (boolean) R named list (names = nodes), crms: (string character R list) R named list (names = nodes), tfs: ((string character R list) R named list (names = crms of the node)) R named list (names = nodes))

e.g. non-TF node G having 2 CRMs C1 with TFs (t11, t12), C2 with TF (t21)
crms = list(nodes=c("G", "t11", "t12", "t21"), isTF=list(G=FALSE, t11=TRUE, t12=TRUE, t21=TRUE), crms=list(G=list("C1", "C2"), t11=NULL, t12=NULL, t21=NULL), tfs=list(G=list(C1=list("t11", "t12"), C2=list("t21")), t11=NULL, t12=NULL, t21=NULL))

### Usage

Write a **Generic Model OBJECT** *model* as a RE:IN file called "model.net":

`R`

`> source("io_expansion.R")`

`> writeModel(model, title="model.net", format="RE:IN")`


Write a **Generic Model OBJECT** as an "expanded model" file called "model.net":

`R`

`> source("io_expansion.R")`

`> writeModel(model, title="model.net", format="expanded")`


Convert a RE:IN file called "model.net" into a **Generic Model OBJECT**:

`R`

`> source("io_expansion.R")`

`> readModel("model.net")`


Expand fully a RE:IN model called "model.net", using a **Regulatory Module Set OBJECT** *crms* and write the result in a RE:IN file called "model_expanded.net":

`R`

`> source("expansion.R")`

`> model <- expansion(crms, filename="model.net", title="model_expanded", type="full", inferTF=FALSE, format="RE:IN")`


Partially:

`R`

`> source("expansion.R")`

`> model <- expansion(crms, filename="model.net", title="model_expanded", type="partial", inferTF=FALSE, format="RE:IN")`


With TF binding inference:

`R`

`> source("expansion.R")`

`> model <- expansion(crms, filename="model.net", title="model_expanded", type="full", inferTF=TRUE, format="RE:IN")`


As an "expanded model" file:

`R`

`> source("expansion.R")`

`> model <- expansion(crms, filename="model.net", title="model_expanded", type="full", inferTF=FALSE, format="expanded")`


## Regulatory Module Analysis

Compute the overlap matrix from a **Regulatory Module Set OBJECT** *crms*:

`R`

`> source("rm_analysis.R")`

`> matrix <- compute_overlap_matrix(crms)`


Build a graph (in DOT) associated with a RE:IN model called "model.net":

`R`

`> source("rm_analysis.R")`

`> build_graph(name="model")`


Build the matrix associated with a RE:IN model called "model.net":

`R`

`> source("rm_analysis.R")`

`> build_matrix(name="model")`


## SBML to IGRAPH/REIN 

To convert a SBML model "model.xml" into a igraph object:

`R`

`> source("sbml2rein.R")`

`> igraph_object <- sbml2rein(filename="model.xml", igraph=TRUE)`


To convert a SBML model "model.xml" into a RE:IN file with "sync" or "async" updates, maximum length of experiments *n*, uniqueness condition "interactions"/"full"/"paths", etc.:

`R`

`> source("sbml2rein.R")`

`> generic_model_object <- sbml2rein(filename="model.xml", updates="async", length=n, ...)`


**Clémence Réda, 2018**
