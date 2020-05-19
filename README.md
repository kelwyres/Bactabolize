# Notes
* Genes can be considered orthologs on the basis of 25% of the protein length
    - not sure if we want to increase requirement for this
    - currently filtering alignments with <25% coverage
    - if a pair of proteins have the highest PID to each other, then they're orthologs
* When capturing unannotated genes with BLASTn, we're only doing a one-way alignment
    - we probably should be taking this hits and alignment them against the reference
        - aligning translated sequence most appropriate here I would think
    - then checking if they qualify as orthologs by definition
* BLASTn hits are only checked for premature stop codon
    - should also check start codon is good
* Implementing as python package rather than nextflow
    - easier to end-users to install
    - writing in a way that this can easily be turned into a nf pipeline
        - idea would be to complement the standalone package


# Current implementation
* Read in model and create set of genes that we will search for in the isolate
* Extract only model gene protein sequences from reference
* Extract all gene protein sequences from isolate
* Bi-directional BLASTp for model proteins and isolate proteins
    - filter on evalue <= 1e-3, coverage >= 25%, pident >= 80%
* Discover orthologs
    - defined as a protein pair that are most similar to each other (by pident)
* Collect model genes that have no ortholog to search at nucleotide level
* Uni-directional BLASTn for unannotated gene detection


# TODO
* Getting 1182 orthologs on K\_quasi\_quasi\_01A030T, tut notes have 1186
    - 1185 after unannotated gene search
    - resolve differences
* Also detecting unannotated genes in K\_quasi\_quasi\_01A030T where tut did not
* Test BLASTn alignment sequence extraction
* Have orthologs.identify create a set of models genes that do NOT have orthologs in isolate
* Determine better scheme for renaming unannotated genes in the new isolate model
    - model genes renamed to match the isolate gene locus\_tag, unnanntated does not have this info


# Planned features
* Validate reference model and genbank
    - primary concern right now are genes present in the model but are missing from gbk
    - some basic code to do this currently but would be better to explicitly do this
* Preferentially take assembly graphs as input
    - to enable collection of deadend count
    - deadends likely important wrt to quality of model
* Reannotate inputs
    - using prodigal model trained on input set
        - computionally intensive for large inputs, optimise conditionally required
* Faster alternative to BLASTp
    - e.g. diamond
    - will need to demonstrate consistency between results
    - must also ensure that reasonable speed up is obtained
* Provide extensive information for draft genomes that fail to optimise
    - determine essential genes from reference and have these output somewhere
        - required only to be generated once (not per strain)


# Planned implementation
* Input as graphs or flat formats (FASTA, genbank, etc)
    - prefering graphs for additional QC
    - (deadends to be demonstrated as important through analysis later)
* QC of assemblies
* Annotation using same prodigal model
* Ortholog detection
* Draft model creatation
* Model optimisation
    - failing models set aside for manual investigation
    - large amount of info provide to help this process
* Growth simulation of specified media
    - should this run if some models failed? I think so


# Planned analysis
* Examine aspects of input assemblies that affect model quality
    - dataset with good quality genomes and biolog data
    - generate assemblies that have varying quality
        - deadends, contigs, completeness, etc
    - create models and test failures
        - quantify on basis of
            - number of genes removed
            - concordance of simulated growth to biolog data
