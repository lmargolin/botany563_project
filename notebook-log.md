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


