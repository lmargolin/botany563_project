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
    
## 2026-03-03
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
    
    
## 2026-03-19
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

## 2026-03-25
- Went through pika_raw_sequences and made a new file, pika_11_raw_sequences that only has 1 sequence per species (total of 11 instead of 250)
- Redo Alignment - MAFT
    - Method: MAFFT (--auto)
    -  Algorithm: FFT-based progressive alignment with iterative refinement
    - Assumptions: Sequences are homologous, evolution is mostly vertical
    - Limitations: May misalign highly divergent regions; guide-tree errors propagate
    - command: 
        mafft --maxiterate 1000 --localpair pika_11_raw_sequences.fasta > pika_11_aligned_MAFT.fasta
    - reasoning: The L-INS-i algorithm (--localpair --maxiterate 1000) was chosen because it maximizes alignment accuracy for small datasets by using local pairwise alignment and extensive iterative refinement, which is appropriate for ~11 homologous sequences with potential indels. This approach prioritizes accuracy over speed, aligning with MAFFT’s design for producing reliable alignments in moderately sized, single-gene datasets like this one.
- Redo Alighment: ClustalW
    - Algorithm: Progressive alignment based on a guide tree
    - Assumptions: Sequences are homologous; early alignments are correct
    - Limitations: Errors early in the guide tree propagate; less accurate for divergent sequences
    - command: 
        clustalw -INFILE=pika_11_raw_sequences.fasta -TYPE=DNA -OUTFILE=pika_11_clustalw_aligned.fasta -OUTPUT=FASTA
    - reasoning: ClustalW was used with the DNA sequence setting to ensure appropriate nucleotide substitution scoring for the ALB gene dataset. Its progressive alignment algorithm is suitable for small datasets of closely related sequences, providing a computationally efficient method for generating a reasonable multiple sequence alignment.
- Redo Distance and Parsimony 
    - With MAFT Alignment 
        1) Load Packcages   
            library(ape)
            library(adegenet)
            library(phangorn)
        2) Loading the sample data and concert to pangorn object
            dna <- read.dna("pika_11_aligned_MAFT.fasta", format="fasta")
            dna2 <- as.phyDat(dna)
        3) We need a starting tree for the search on tree space and compute the parsimony score of this 
            tre.ini <- nj(dist.dna(dna, model="raw"))
            parsimony(tre.ini, dna2)
                - result: [1] 8
        4) search for the tree with maximum parsimony 
            tre.pars <- optim.parsimony(tre.ini, dna2)
                - result: n/a
        5) plot tree
            plot(tre.pars, cex=0.6)
                - plot looked good as is
    - With ClustalW Aliightnment 
        1) Load Packcages   
            library(ape)
            library(adegenet)
            library(phangorn)
        2) Loading the sample data and concert to pangorn object
            dna <- read.dna("pika_11_clustalw_aligned.fasta", format="fasta")
            dna2 <- as.phyDat(dna)
        3) We need a starting tree for the search on tree space and compute the parsimony score of this 
            tre.ini <- nj(dist.dna(dna, model="raw"))
            parsimony(tre.ini, dna2)
                - result: [1] 8
        4) search for the tree with maximum parsimony 
            tre.pars <- optim.parsimony(tre.ini, dna2)
                - result: n/a
        5) plot tree
            plot(tre.pars, cex=0.6)
                - plot looked good as is, this one had longer arms and less labels (no species name) compared to the MAFT one
            - note* liked the look of maft better, so using that moving forwards
- Redo Maximum Likelihood
    - Maximum Likelihood (RAxML)
    - Description of RAxML: RAxML-NG (Randomized Axelerated Maximum Likelihood – Next Generation) is a widely used software tool for inferring phylogenetic trees using the maximum likelihood (ML) framework. It estimates the tree topology, branch lengths, and substitution model parameters that maximize the likelihood of observing a given multiple sequence alignment under a specified model of sequence evolution. RAxML-NG is optimized for large datasets and supports parallel computation, making it efficient for analyses with many taxa and sites.
    - Assumptions of RAxML: RAxML assumes that the input sequences are correctly aligned and homologous, meaning each column represents a shared evolutionary position. It relies on an explicit substitution model (e.g., LG+G+F for protein data) that assumes a particular pattern of sequence evolution, including stationarity, reversibility, and homogeneity across the tree (unless partitioned models are used). Sites are typically assumed to evolve independently, with rate heterogeneity modeled (e.g., via a Gamma distribution).
    - Limitations of RAxML: RAxML’s accuracy depends strongly on the quality of the input alignment; misaligned regions or excessive gaps can bias results. Violations of model assumptions (e.g., compositional bias, heterotachy, or non-independence among sites) can lead to incorrect tree inference. The ML search is heuristic, so it does not guarantee finding the global optimum tree, especially for large or complex datasets. Additionally, duplicated or highly similar sequences can influence branch length estimation and support values if not handled appropriately.
    1) check the alignment
        raxml-ng --check --msa pika_11_aligned_MAFT.fasta --model LG+G8+F
    2) find the ML tree
        raxml-ng --msa pika_11_aligned_MAFT.fasta --model LG+G8+F
    3) Open R and dirty plot
        library(ape)
        tre = read.tree(file="pika_11_aligned_MAFT.fasta.raxml.bestTree")
        plot(tre)
    4)Non-parametric bootstrapping (in Terminal, not R)
        raxml-ng --all --msa pika_11_aligned_MAFT.fasta --model LG+G8+F --bs-trees 10 --prefix pika_11_aligned_MAFT.bootstrap
    5) Open R and dirty plot
        library(ape)
        tre = read.tree(file="pika_11_aligned_MAFT.bootstrap.raxml.bestTree")
        plot(tre)
        nodelabels(tre$node.label)
    6)note the tree does not seem to be rooted correctly (R)
        library(ape)
        tre = read.tree(file="pika_11_aligned_MAFT.bootstrap.raxml.bestTree")
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
        iqtree -s pika_11_aligned_MAFT.fasta
    2) Open R and run a dirty plot
        library(ape)
        tre = read.tree(file="pika_11_aligned_MAFT.fasta.treefile")
        plot(tre)
    3) Root the tree again
        plot(tre)
        nodelabels()

        rtre = root(tre, node=31, resolve.root=TRUE)
        plot(rtre)
    4) Quantify support for the estimated tree. Add a prefix because it does not let us overright the original files produced. (terminal, not R)
        iqtree -s pika_11_aligned_MAFT.fasta -m GTR+G4 -b 10 -pre pika_11_aligned_MAFT.fasta-iqtree-bootstrap
    5) In R, plot the tree again with bootstrap support
        library(ape)
        tre = read.tree(file="pika_11_aligned_MAFT.fasta-iqtree-bootstrap.treefile")
        plot(tre)
        nodelabels()

        rtre = root(tre, node=31, resolve.root=TRUE)
        plot(rtre)
        nodelabels(rtre$node.label)


## 2026-04-07
- Working on Bayesian Method via Mr Bayes
- Bayesian Inference: Mr Bayes
- Description: Phylogenetic relationships were inferred using MrBayes, which estimates the posterior distribution of trees and model parameters with a Markov chain Monte Carlo (MCMC) algorithm. Rather than producing a single best tree only, this approach samples many possible trees in proportion to their posterior probability, allowing clade support to be expressed as posterior probabilities. For this dataset, an HKY substitution model with gamma-distributed rate variation among sites was used, with two independent runs of four chains each.
- Assumptions: MrBayes assumes that the input sequences are homologous and correctly aligned, that the chosen substitution model adequately describes sequence evolution, and that sites evolve independently according to the specified model. It also assumes that the MCMC chains converge to the target posterior distribution and mix sufficiently well to provide reliable parameter and tree estimates.
- Limitations: Bayesian phylogenetic inference can be sensitive to model choice, prior settings, and insufficient chain length or poor convergence. Posterior probabilities may appear strongly supported even when the underlying dataset contains limited information. In this analysis, a major limitation is the short alignment length (70 bp), which restricts phylogenetic signal and may reduce confidence in the inferred relationships. Therefore, although convergence diagnostics indicated that the analysis ran successfully, the resulting tree should still be interpreted cautiously because uncertainty may reflect the limited amount of sequence data rather than algorithm performance.
- Prepare input file (input file is pika_11_aligned_MAFT.fasta)
    pip install biopython
    
    python - <<EOF
    from Bio import SeqIO
    
    records = list(SeqIO.parse("pika_11_aligned_MAFT.fasta", "fasta"))
    
    for r in records:
        r.annotations["molecule_type"] = "DNA"

    SeqIO.write(records, "pika_11.nex", "nexus")
    EOF
- Add Mr Bayes code to nexus file (full file looks like this)
    nano pika_11.nex 
    
        #NEXUS
        begin data;
        dimensions ntax=11 nchar=70;
        format datatype=dna missing=? gap=-;
        matrix
        KP292978.1 agaagtgctgtgctgctgctgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcca
        KP292980.1 agaagtgctgtgctgctgatgacaaggaagcctgcttttcagaggaggtactggagctgtgttccytcca
        KP292982.1 agaagtgctgtgctgctgctgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcca
        KP292984.1 agaagtgctgtgctgctgctgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcca
        KP292986.1 agaagtgctgtgctgctgctgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcca
        KP292988.1 agaagtgctgtgctgctgctgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcca
        KP292990.1 agaagtgctgtgctgctgatgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcca
        KP292991.1 agaagtgctgtgctgctgacgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcta
        KP292993.1 agaagtgctgtgctgctgctgacaaggaagcctgctttttagaggaggtactggagctgtgttccctcaa
        KP292995.1 agaagtgctgtgctgccgctgacaaggaagcctgcttttcagaggaggtactggagctgtgctcccttca
        KP292997.1 agaagtgctgtgctgctgctgacaaggaagcctgcttttcagaggaggtactggagctgtgttccctcca
        ;
        end;

        begin mrbayes;
            set autoclose=yes nowarn=yes;
            prset brlenspr=unconstrained:exp(10.0);
            prset shapepr=exp(1.0);
            prset tratiopr=beta(1.0,1.0);
            prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
            lset nst=2 rates=gamma ngammacat=4;
            mcmcp ngen=1000000 samplefreq=100 printfreq=100 diagnfreq=1000 nruns=2 nchains=4 savebrlens=yes;
            mcmc;
            sump burnin=2500;
            sumt burnin=2500;
        end;
- reasoning for the parameters chosen 
    Bayesian phylogenetic inference was performed using MrBayes. A General Time Reversible (GTR) substitution model with gamma-distributed rate variation among sites was employed. Two independent runs of four Markov chains were conducted for 1,000,000 generations, sampling every 100 generations. Convergence was assessed using the average standard deviation of split frequencies (<0.01) and potential scale reduction factor (PSRF ≈ 1). The first 25% of samples were discarded as burn-in. Because the alignment was short (70 bp), results should be interpreted with caution.
- run mr bayes
    mb pika_11.nex
- thoughts on run 
    The Bayesian analysis converged successfully, with an average standard deviation of split frequencies below 0.01 and PSRF values near 1.0. Effective sample sizes (ESS) were high for all parameters, indicating adequate sampling. However, posterior probabilities for most clades were moderate to low, and a large number of trees were present in the credible set, suggesting limited phylogenetic resolution due to the short alignment length (70 bp).
- important output files
    - final tree -> pika_11.nex.con.t
    - convergence and stats -> pika_11.nex.pstat
- Tracer - not done yet? does this need to be completed?

## 2026-04-19 and 2026-04-20
- The Coalescent - ASTRAL
- Description: Phylogenetic relationships were inferred using ASTRAL, a summary method that estimates a species tree from a collection of gene trees under the multispecies coalescent model. Rather than concatenating sequences into a single supermatrix, ASTRAL uses the topologies of individual gene trees to infer the species tree that maximizes agreement among induced quartet relationships (sets of four taxa). This approach explicitly accounts for gene tree discordance, which can arise due to incomplete lineage sorting. In this analysis, one gene tree was inferred per nuclear marker (using maximum likelihood), and these gene trees were combined as input to ASTRAL to produce a species tree with branch support values based on quartet frequencies.
- Assumptions: ASTRAL assumes that the input gene trees are estimated from orthologous loci and are reasonably accurate representations of the true gene histories. It is based on the multispecies coalescent model, which assumes that discordance among gene trees is primarily due to incomplete lineage sorting rather than other processes such as horizontal gene transfer, gene duplication, or hybridization. The method also assumes that taxa are correctly labeled and consistently represented across gene trees, and that each gene tree is unrooted. Additionally, it assumes that loci are independently inherited and that the sampling of loci is sufficient to capture the underlying species tree signal.
- Limitations: Because ASTRAL is a summary method, its accuracy depends heavily on the quality of the input gene trees; errors in gene tree estimation (e.g., due to short alignments, poor model fit, or limited phylogenetic signal) can propagate into the species tree. The method does not use the original sequence data directly, so it cannot correct for biases introduced during gene tree inference. ASTRAL also assumes that incomplete lineage sorting is the primary cause of discordance and does not explicitly model other evolutionary processes such as gene duplication/loss or introgression, which may lead to misleading results if present. In this analysis, a key limitation is the relatively small number of loci (12 nuclear markers), which may limit the ability to fully resolve species relationships and reduce confidence in branches with low quartet support. Therefore, while the inferred species tree provides a coalescent-based estimate of relationships, it should be interpreted alongside results from other methods (e.g., concatenated maximum likelihood or Bayesian inference) to assess consistency and robustness.
- downloaded ASTRAL from zip file (https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip)
- test my installation - it worked
    java -jar /Users/liamargolin/Downloads/Astral/astral.5.7.8.jar
    
    java -jar /Users/liamargolin/Downloads/Astral/astral.5.7.8.jar \
        -i /Users/liamargolin/Downloads/Astral/test_data/song_primates.424.gene.tre \
        -o /Users/liamargolin/Downloads/astral_test_output.tre
        
    ls /Users/liamargolin/Downloads/astral_test_output.tre
    
    cat /Users/liamargolin/Downloads/astral_test_output.tre
- The available 11-taxon FASTA file contained only one marker (ALB), so it was not suitable for a multispecies coalescent analysis in ASTRAL, which requires multiple gene trees from multiple loci.
- Using sample dataset to practice
    in directory -> /Users/liamargolin/Downloads/Astral
    
    java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre
    java -Djava.library.path=./lib/ -jar astralmp.5.7.8.jar -i test_data/song_mammals.424.gene.tre
    
    java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre
    java -Djava.library.path=./lib/ -jar astralmp.5.7.8.jar -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre
    
    *in R - read tree
        library(ape)
        tre = read.tree(file="song_mammals.tre")
        plot(tre)
        nodelabels(text = tre$node.label)
    
    *got an error - ran this instead per chat gpt - got a tree
        tre <- read.tree("/Users/liamargolin/Downloads/Astral/test_data/song_mammals.tre")
        plot(tre)
        nodelabels()

 

