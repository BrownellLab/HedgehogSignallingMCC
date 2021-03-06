## Exploration of Hedgehog signaling in MCC samples

### Publications

Carroll, T. M. et al. **Hedgehog signaling inhibitors fail to reduce Merkel cell carcinoma viability.** *Journal of Investigative Dermatology* [doi:10.1016/j.jid.2017.01.008](http://dx.doi.org/10.1016/j.jid.2017.01.008)

### Introduction

We use merkel cell carcinoma ([GSE50451](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50451)) basal cell carcinoma ([GSE7553](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7553)) and medulloblastoma ([GSE37418](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37418)) samples as comparison.

See the [reports](reports) directory for the results of analysis.

Note that the datasets are on separate platforms (hgu133plus2 and hgu133a2) - we keep probes in common and then a $z$-score normalization on individual data sets.

### Project Information

This project was created using `ProjectTemplate`.

To load this project, you'll first need to `setwd()` into the directory
where this README.md file is located. Then you need to run the following two
lines of R code:

```
library('ProjectTemplate')
load.project()
```

After you enter the second line of code, you'll see a series of automated
messages as ProjectTemplate goes about doing its work. This work involves:

* Reading in the global configuration file contained in `config`.
* Loading any R packages you listed in he configuration file.
* Reading in any datasets stored in `data` or `cache`.
* Preprocessing your data using the files in the `munge` directory.

For more details about ProjectTemplate, see http://projecttemplate.net
