## Supplementary material

These are additional material for the tests performed in the core paper.

### CISVIEW_Collombet.pdf

This is the file where information and analysis of TF bindings in the Collombet dataset [Collombet et al., 2017] extracted from the CisView database [Sharov et al., 2006] is gathered. 

*collombet_crms_cisview.R* 

This is the R file used to expand the Collombet model in the algorithm.

### CISVIEW_Dunn.pdf 

This is the file where information and analysis of TF bindings in the Dunn dataset [Dunn et al., 2014] extracted from the CisView database [Sharov et al., 2006] is gathered.

*dunn_crms_cisview.R*

This is the R file used to expand the Dunn model in the algorithm.

### "collombet" and "dunn" folders

Contains the first solution models for (resp. non) expanded models (Collombet et al.'s and Dunn et al.'s): resulting (simplified) GRFs and selected interactions in the abstract model, associated with the figures presented at the end of Appendix in the core paper.

### Commands (to get the results presented in the paper)

To be run in the root folder of the GitHub.

*Collombet et al.'s model*

**Expansion procedure**
```R
R -e 'setwd("R/"); source("expansion.R"); setwd("../examples/models/collombet/"); source("collombet_crms_cisview.R"); model_regular <- readModel("model.net"); writeModel(model_regular, title="model_expanded.net", format="expanded"); model_expanded <- expansion(crms, "model.net", title="collombet-partial+_expanded.net", type="partial", method="positive", inferTF=F, format="expanded")'
```

**Network inference (expanded model)**

```python
python solve.py run --simplify --visualize collombet/collombet-partial+_expanded.net collombet/observations.spec > results/report-collombet-partial+_expanded.txt
```

*Dunn et al.'s model*

**Expansion procedure**
```R
R -e 'setwd("R/"); source("expansion.R"); setwd("../examples/models/pluripotency/"); source("pluripotency_crms_cisview.R"); model_regular <- readModel("model.net"); writeModel(model_regular, title="model_expanded.net", format="expanded"); model_expanded <- expansion(crms, "model.net", title="pluripotency-partial+_expanded.net", type="partial", method="positive", inferTF=F, format="expanded")'
```

**Network inference (expanded model)**

```python
python solve.py run --simplify --visualize pluripotency/pluripotency-partial+_expanded.net pluripotency/observations.spec > results/report-pluripotency-partial+_expanded.txt
```
