Using data from our 23 MCC samples (platform HG-U133A v2), a dataset including 15 basal cell carcinoma (GSE7553), and a dataset of 76 medulloblastoma samples, 10 of which are classified as SHH-driven (GSE37418), I created a subset of data for SHH pathway genes using the KEGG pathway (hsa04340) as a reference.
To remove redundant and potentially uninformative probes, I first used all probes to derive expression heatmap plots for the each dataset individually.
The baseline was chosen to be the medulloblastoma samples, since there are both SHH and non-SHH driven samples.
For example, there were multiple probes for PTCH1, but one of them exhibited consistently high expression in all datasets across all samples.
The list of genes chosen after filtering is seen in Table 1.

A summary value for each gene was made by first averaging expression values across probes for each gene.
Then, the average over all samples was taken for each gene.
This value was used as input to a false color gene network made in Cytoscape - red indicates higher expression, blue indicates lower.
The middle of each color scale is set independently for each dataset to the mean of all expression values in the entire dataset (not limited to hedgehog pathway genes).

Symbol and Affymetrix probe ID of genes in SHH pathway from KEGG:

Affy ID | Gene Symbol 
---|---
208570\_at | WNT1 
218629\_at | SMO 
207586\_at | SHH 
209816\_at | PTCH1 
205710\_at | LRP2 
205201\_at | GLI3 
208057\_s\_at | GLI2 
207034\_s\_at | GLI2 
206646\_at | GLI1 
204457\_s\_at | GAS1 
204456\_s\_at | GAS1 
222374\_at | BTRC 
204901\_at | BTRC 
211518\_s\_at | BMP4 

Update figures here:

#+CAPTION: The Homo sapiens hedgehog pathway from KEGG (hsa04340).
#+LABEL: fig:kegghedgehogpathway
[[./hsa04340.png]]

#+CAPTION: Heatmap of expression of genes in hedgehog pathway from basal cell carcinoma samples.
[[../graphs/heatmap_basal.png]]

#+CAPTION: Average expression of genes in hedgehog pathway from basal cell carcinoma samples.
[[./hedgehog_basal.png]]

\newpage

#+CAPTION: Heatmap of expression of genes in hedgehog pathway from medulloblastoma samples.
[[../graphs/heatmap_medullo.png]]

#+CAPTION: Average expression of genes in hedgehog pathway from medulloblastoma, SHH subgroup samples.
[[./hedgehog_medullo_SHH.png]]

\newpage

#+CAPTION: Heatmap of expression of genes in hedgehog pathway from Merkel cell carcinoma samples.
[[../graphs/heatmap_mcc.png]]


#+CAPTION: Average expression of genes in hedgehog pathway from Merkel cell carcinoma samples.
[[./hedgehog_mcc.png]]


