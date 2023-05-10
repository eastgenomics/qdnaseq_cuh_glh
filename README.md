# qdnaseq_cuh_glh

This repo contains several script to run qDNASeq tool. Each script has it's own purpose which will be described below.

## qDNASeq_analysis.Rmd

This script compares the outputs from qDNASeq to a truth WGS dataset. This script only requires three main inputs for the analysis to be modified or added, this is:
* Metadata - Contains rows per sample and columns as case, CUH_number, sample name and tissue type
* CUH number
* Qualimap depth - Contains rows per sample and columns as sample name and qualimap depth.

There is a analysis_done with yes or no input which can be modified to rerun the analysis, this is because running one analysis for a case can take from 1 hour to 5 mins depending on the bin size (the smaller the bin, the longer the analysis takes). The script is usually run from the root project directory. The subfolders should be the qDNASeq outputs at a different depth and bin size. In the sWGS folders should be the truth dataset that starts with the CUH number. An example of the output directory is summarized below:

```bash
├── sWGS
    ├── case5_down_100bin
    ├── case5_down_15bin
    ├── case5_down_500bin
    ├── case5_down_50bin
    ├── case5_down_5bin
    ├── case5_merged_100bin
    ├── case5_merged_15bin
    ├── case5_merged_500bin
    ├── case5_merged_50bin
    ├── case5_merged_5bin
    ├── case5_norm_100bin
    ├── case5_norm_15bin
    ├── case5_norm_500bin
    ├── case5_norm_50bin
    ├── case5_norm_5bin
```