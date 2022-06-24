# A computationally-enhanced hiCLIP atlas reveals Staufen1 RNA binding features and links 3’ UTR structure to RNA metabolism

The repo contains code for our manuscript: [**"A computationally-enhanced hiCLIP atlas reveals Staufen1 RNA binding features and links 3’ UTR structure to RNA metabolism"**](https://www.biorxiv.org/content/10.1101/2022.06.13.495933v1).

*Tosca*, our Nextflow RNA proximity ligation data analysis pipeline, is available at: https://github.com/amchakra/tosca.

## Abstract

The structure of mRNA molecules plays an important role in its interactions with trans-acting factors, notably RNA binding proteins (RBPs), thus contributing to the functional consequences of this interplay. However, current transcriptome-wide experimental methods to chart these interactions are limited by their poor sensitivity. Here we extend the hiCLIP atlas of duplexes bound by Staufen1 (STAU1) ~10-fold, through careful consideration of experimental assumptions, and the development of bespoke computational methods which we apply to existing data. We present Tosca, a Nextflow computational pipeline for the processing, analysis and visualisation of proximity ligation sequencing data generally. We use our extended duplex atlas to discover insights into the RNA selectivity of STAU1, revealing the importance of structural symmetry and duplex-span-dependent nucleotide composition. Furthermore, we identify heterogeneity in the relationship between STAU1-bound 3' UTRs and metabolism of the associated RNAs that we relate to RNA structure: transcripts with short-range proximal 3' UTR duplexes have high degradation rates, but those with long-range duplexes have low rates. Overall, our work enables the integrative analysis of proximity ligation data delivering insights into specific features and effects of RBP-RNA structure interactions.

## Structure

* `figures` contains code to generate the figures
* `linker` contains code for the linker hybrids
* `no_linker` contains code for the no-linker hybrids
* `no_rnase` contains code for the no RNase hybrids (short-range structures)
* `paris` contains code for the PARIS hybrids
* `ref` contains code to generate the reference files
