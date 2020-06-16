# Notes
* Implementing as python package rather than nextflow pipeline
    - easier for end-users to install
    - writing in a way that this can easily be turned into a nf pipeline
        - idea would be to complement the standalone package if people wanted it
* Script is currently designed to process a single isolate per command
    - this can be parallelised itself with comparatively low memory overhead
    - only additional computation in this approach is creating BLAST databases for the reference
        - this is negligible though


# Queries
* Genes can be considered orthologous from only 25% of the protein sequence
    - currently filtering alignments with <25% coverage
    - if a pair of proteins have the highest PID to each other, then they're orthologs
    - not sure if we want to increase coverage requirement for this
* When capturing unannotated genes with BLASTn, we're only doing a one-way alignment
    - we probably should be taking these hits and align them against the reference
        - aligning translated sequence is most appropriate here I would think
    - then checking if they qualify as orthologs by the standard definition
* BLASTn hits are only checked for premature stop codon
    - I think we should also check start codon is good
    - the whole approach of identifying genes that couldn't be annotated is dicey


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
    - filter on evalue <= 1e-3, coverage >= 80%, pident >= 80%
    - require translated sequence to not have a truncating mutation
* Any BLASTn hit that passes filtering is automatically considered orthologous
    - I don't think this is the best approach
* Remove models genes that do not have an ortholog in the isolate
    - Artificial genes are excepted here
* Rename identified orthologs to match locus\_tags in isolate
* Write model to disk


# TODO
* Check required programs are in PATH and have correct versions
* Getting 1182 orthologs on K\_quasi\_quasi\_01A030T, tut notes have 1186
    - 1185 after unannotated gene search (tut detected no annotated genes)
    - resolve differences
* Test BLASTn alignment sequence extraction
    - tested this to some extent and I think it's correct but tut notes have different offsets...
* Determine better scheme for renaming unannotated genes in the new isolate model
    - model genes renamed to match the isolate gene locus\_tag, unnannotated does not have this info
    - currently naming as '{old\_name}\_unannotated'
* Validate input model format (json)
    - check all genes are present in reference?
        - this may not be desired - certainly not for pan-models
* Fix help message
    - all arguments appear in optional section
* Test quality of models using flux balance analysis
    - using a objective function and a matrix of metabolite stoichiometric coeffs.
        - stoichiometric coefficients should be experimentally derived
        - object function is just a vector of linear objective coefficients
            - default is 1
        - combined create system to linear eqs. to solve flux for different components
    - calculates flow of metabolites through metabolic network
        - specifically optimises for reaction rates (fluxes)
    - enable prediction of growth rate of organism or production rate of specific metabolite
    - for our purposes objective function will return growth rate
    - presumably growth media simulation takes the same process but constrains available precursors
        - we should then use a simple FBA to quickly perform shallow validation of critical pathways
            - i.e. validate on minimal media before proceeding to FBA for specific media types
* Troubleshooting
    - gather all essential metabolic genes from model
        - this should only be done once
        - must be done in a way that is decoupled from the current single-isolate command
        - in the future, the current single-isolate command will be wrapped by other code
            - which manages execution of draft model creation
        - check which are missing
    - collect essential reactions
        - see above notes
        - discover missing
    - the two above points appear to be done by cobra gapfill
        - might be reaction-specific and not indicate which genes are missing
        - solver is heuristic - run several times to check for convergence/ consistency
    - retain all BLAST results
    - check for more distance homologs
    - run Bandage graph blast
        - help diagnose bad assembly
    - check in BiGG database for different genes have the same functionality
    - report known problem genes such as those in capsule loci


# Planned features
* Validate reference model and genbank
    - primary concern right now are genes present in the model but are missing from gbk
    - some basic code to do this currently but would be better to explicitly do this
* Preferentially take assembly graphs as input
    - to enable collection of deadend count
    - deadends likely important wrt to quality of model
* Reannotate inputs
    - using prodigal model trained on input set
        - computationally intensive for large inputs, optimise conditionally required
* Faster alternative to BLASTp
    - e.g. diamond
    - will need to demonstrate consistency between results
    - must also ensure that reasonable speed up is obtained
* Provide extensive information for draft genomes that fail to optimise
    - determine essential genes from reference and have these output somewhere
        - required only to be generated once (not per strain)


# Planned implementation
* Input as graphs or flat formats (FASTA, genbank, etc)
    - preferring graphs for additional QC
    - (deadends to be demonstrated as important through analysis later)
* QC of assemblies
* Annotation using same prodigal model
* Ortholog detection
* Draft model creation
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
