## A Model of Episodic Directional Selection (MEDS)
MEDS.bf implements the Model of Episodic Directional Selection described in [Murrell et al](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002507). MEDS uses a codon model with a directional component along a priori selected foreground branches.

## Files
The data and simulations used in the manuscript can be found in MEDSfiles.bz2

## Labeling Foreground Branches
The foreground branches of the phylogeny must be labelled. To do this, "{FG}" is placed after the sequence name (but before the colon) in the newick tree string. Here is an example (sub)tree with "Branch1" labelled foreground:
```
(Branch1{FG}:0.1,Branch2:0.1)
```
Editing the newick files by hand is inefficient and confusing. We suggest using [phylotree.js](http://phylotree.hyphy.org/). Detailed instructions for labeling branches can be found [here](http://phylotree.hyphy.org/).


## Running the Analysis
To run a MEDS analysis, a nucleotide coding sequence alignment and a rooted newick tree with tagged foreground are needed. Executing the batch file will prompt for the tree and data file. The underlying nucleotide substitution model is REV, but any 6 character nucleotide model string (010010 for HKY85, 012345 for REV) can be specified by modifying the batch file. MEDS can take quite some time and should be run in HyPhy's command line mode. A small reverse transcriptase demo alignment (MEDSdemo.fas) and phylogeny (MEDSdemo.nwk) are available.

## The Output
The output for a MEDS analysis is a file containing, for each site, the maximized likelihood values for a null model, an episodic diversifying model (FEEDS), and 21 MEDS analyses: one for each amino acid (including stop codons, although these never occur in coding sequences, so that test is never run). For computational convenience, the directional model is only optimized for amino acids that are actually observed at the site in question, and omission codes (-99) are used as place holders for parameter values for amino acids that don't occur. Don't be alarmed if these appear almost everywhere: only a handful of amino acids occur at any particular site.

Also, a re-parameterized version of the ωT parameters are reported. ωT is the target amino acid rate multiplier. To improve optimization speed, we re-parameterized ωT = ((1/(1-xT))-1), and optimized the xT parameters. The values reported in the output file are actually the maximum likelihood xT values.

To make things simple, a python script for processing the results files is available. The user inputs the output file from a MEDS analysis and selects a p-value threshold with:

```
python MEDSproc.py <resultsFile.csv> <alphaLevel>
```
All sites and substitutions with significant directional or diversifying selection are returned, along with the corresponding maximum likelihood parameter values for the null and alternative models. The xT parameters are also transformed into ωT values.
