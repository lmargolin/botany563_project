# Project Notebook Log 

## 2026-01-30
- Created repository
- Set up folder structure 

## 2026-02-10
- Picked Data Set
- I plan on working with 11 species of pikas (ochotona). This dataset contains sequences from 12 loci that were aligned with ClustalW. Here is the link to the dataset on DRYAD: https://datadryad.org/dataset/doi:10.5061/dryad.bn547#usage 
- Quality Control: QC was already preformed when data was accesses. Raw reads were not provided. No further action needed to be completed. 
    
## 2026-02-16
- Working on Alignment 
- Needed to get raw reads from paper -> GenBank Acc. Nrs. KP292978-KP293227
- Installed EDirect to use NCBI's command-line toos to produce a FASTA file pika_raw_sequences.fasta that contains all unaligned sequences in that GenBank range
- code
    esearch -db nucleotide -query "KP292978:KP293100[ACCN]" \
        | efetch -format fasta > part1.fasta

    esearch -db nucleotide -query "KP293101:KP293227[ACCN]" \
        | efetch -format fasta > part2.fasta

    cat part1.fasta part2.fasta > pika_raw_sequences.fasta
- double checked by opening the file and saw nucleotides, and the file size was 165KB
- Quality Control: could not run FastQC because my raw data is in FASTA format not FASTQ. Instead used seqkit
        conda install -c bioconda seqkit
        seqkit stats pika_raw_sequences.fasta
- The dataset consists of 250 raw sequences in FASTA format, ranging from 234 to 909 bases (average 570 bp). Quality control was performed using seqkit stats to check sequence length distribution and overall integrity, confirming the sequences are suitable for downstream analyses.
    ## Multiple Sequence Alignment(1)
    - Method: MAFFT (--auto)
    - Algorithm: FFT-based progressive alignment with iterative refinement
    - Assumptions: Sequences are homologous, evolution is mostly vertical
    - Limitations: May misalign highly divergent regions; guide-tree errors propagate
    - command: mafft --auto pika_raw_sequences.fasta > pika_aligned.fasta
    ## Multiple Sequence Alignment using ClustalW (2)
    - Algorithm: Progressive alignment based on a guide tree
    - Assumptions: Sequences are homologous; early alignments are correct
    - Limitations: Errors early in the guide tree propagate; less accurate for divergent sequences
    - command: clustalw -INFILE=pika_raw_sequences.fasta -TYPE=DNA -OUTFILE=pika_clustalw_aligned.fasta -OUTPUT=FASTA
    
##2026-03-03
- Working on Distance and Parsimony 
1) installed necessary packages (R)
    install.packages("adegenet", dep=TRUE)
    install.packages("phangorn", dep=TRUE)
2) got in correct working directory 
    setwd("~/Desktop/botany563_project/data/raw")
3) Loaded packages (R)
    library(ape)
    library(adegenet)
    library(phangorn)
4) Loading the sample data and concert to pangorn object
    dna <- read.dna("pika_aligned.fasta", format="fasta")
    dna2 <- as.phyDat(dna)
5) We need a starting tree for the search on tree space and compute the parsimony score of this tree (422)
    tre.ini <- nj(dist.dna(dna, model="raw"))
    parsimony(tre.ini, dna2)
    - result: [1] 4257
6) search for the tree with maximum parsimony 
    tre.pars <- optim.parsimony(tre.ini, dna2)
    - result: Final p-score 3433 after  87 nni operations
7) plot tree
    plot(tre.pars, cex=0.6)
8) labels were too crowded so tried again
    plot(tre.pars, cex=0.5)
    
    
##2026-03-19
- Maximum Likelihood (RAxML)
- Description of RAxML: RAxML-NG (Randomized Axelerated Maximum Likelihood – Next Generation) is a widely used software tool for inferring phylogenetic trees using the maximum likelihood (ML) framework. It estimates the tree topology, branch lengths, and substitution model parameters that maximize the likelihood of observing a given multiple sequence alignment under a specified model of sequence evolution. RAxML-NG is optimized for large datasets and supports parallel computation, making it efficient for analyses with many taxa and sites.
- Assumptions of RAxML: RAxML assumes that the input sequences are correctly aligned and homologous, meaning each column represents a shared evolutionary position. It relies on an explicit substitution model (e.g., LG+G+F for protein data) that assumes a particular pattern of sequence evolution, including stationarity, reversibility, and homogeneity across the tree (unless partitioned models are used). Sites are typically assumed to evolve independently, with rate heterogeneity modeled (e.g., via a Gamma distribution).
- Limitations of RAxML: RAxML’s accuracy depends strongly on the quality of the input alignment; misaligned regions or excessive gaps can bias results. Violations of model assumptions (e.g., compositional bias, heterotachy, or non-independence among sites) can lead to incorrect tree inference. The ML search is heuristic, so it does not guarantee finding the global optimum tree, especially for large or complex datasets. Additionally, duplicated or highly similar sequences can influence branch length estimation and support values if not handled appropriately.
1) check the alignment
    raxml-ng --check --msa pika_clustalw_aligned.fasta --model LG+G8+F
2) find the ML tree
    raxml-ng --msa pika_clustalw_aligned.fasta --model LG+G8+F
3) Open R and dirty plot
    library(ape)
    tre = read.tree(file="pika_clustalw_aligned.fasta.raxml.bestTree")
    plot(tre)
4)Non-parametric bootstrapping 
    raxml-ng --all --msa pika_clustalw_aligned.fasta --model LG+G8+F --bs-trees 10 --prefix pika_clustalw_aligned.fasta.raxml.boostrap
5) Open R and dirty plot
    library(ape)
    tre = read.tree(file="pika_clustalw_aligned.fasta.raxml.boostrap.raxml.support")
    plot(tre)
    nodelabels(tre$node.label)
6)note the tree does not seem to be rooted correctly
    library(ape)
    tre = read.tree(file="pika_clustalw_aligned.fasta.raxml.boostrap.raxml.support")
    plot(tre)
    nodelabels()

    rtre = root(tre, node=33, resolve.root=TRUE)
    plot(rtre)
    nodelabels(rtre$node.label)


- Maximum Likelihood: IQ-Tree
- Description: IQ-TREE is a phylogenetic inference software that uses the maximum likelihood (ML) framework to estimate evolutionary trees from multiple sequence alignments. It features efficient tree-search algorithms and automated model selection (e.g., ModelFinder), and provides measures of branch support such as ultrafast bootstrap and SH-aLRT.
- Assumptions: IQ-TREE assumes that sequences are correctly aligned and homologous, and that sequence evolution follows a specified substitution model (often assuming stationarity, reversibility, and homogeneity across lineages unless otherwise modeled). Sites are generally treated as independent, with rate variation accommodated (e.g., via Gamma or FreeRate models).
- Limitations: Results depend heavily on alignment quality and the appropriateness of the chosen model. Violations of model assumptions (e.g., compositional bias or heterotachy) can bias tree inference. Like other ML methods, the tree search is heuristic and may not find the global optimum. Some support metrics (e.g., ultrafast bootstrap) can be overconfident under certain conditions if model assumptions are violated.
1) Run
    iqtree -s pika_aligned.fasta
2) Open R and run a dirty plot
    library(ape)
    tre = read.tree(file="pika_aligned.fasta.treefile")
    plot(tre)
3) Root the tree again
    plot(tre)
    nodelabels()

    rtre = root(tre, node=31, resolve.root=TRUE)
    plot(rtre)
4) Quantify support for the estimated tree. Add a prefix because it does not let us overright the original files produced.
    iqtree -s pika_aligned.fasta -m HIVb+G4 -b 10 -pre pika_aligned.fasta-iqtree-bootstrap
5) In R, plot the tree again with bootstrap support
    library(ape)
    tre = read.tree(file="pika_aligned.fasta-iqtree-bootstrap.treefile")
    plot(tre)
    nodelabels()

    rtre = root(tre, node=31, resolve.root=TRUE)
    plot(rtre)
    nodelabels(rtre$node.label)










