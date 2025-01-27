# Metatranscriptomics reveal the entire diversity and tree-specific differences of bark microbiomes and their microbial food webs in the canopy of a floodplain forest
This repository is a collection of scripts that guide you through the data analyses of the paper of [Freudenthal et al. 2024](https://doi.org/10.1093/ismejo/wrae206). The raw sequencing data are available at the the NCBI BioProject [PRJNA1105877](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1105877), the metadata and the trait databases are provided in [OriginalData](Data/OriginalData). Intermediate files, as well as results and figures, can be obtained by running the data analysis scripts.

We show an overview of our [workflow](Scripts/01RawDataProcessing.md) for creating count tables from raw Illumina data (paired-end sequences, 150bp) as well as of our [preprocessing steps](Scripts/02CountDataPreprocessing.md) for the count tables.

Further, we provide all scripts that were used to calculate [rarefaction curves](Scripts/03RarefactionCurves.md), for visualizing community composition with [sankey diagrams](Scripts/04CommunityComposition.md), for [network inference](Scripts/05Network.md), for assessing [alpha](Scripts/06AlphaDiversity.md) and [beta](Scripts/07BetaDiversity.md) diversity and for visualizing the significant differences in the communities across tree species with [heatmaps](Scripts/08Heatmaps.md).

Self-written Functions that are used in the Scripts are stored at [Functions](Scripts/Functions).
