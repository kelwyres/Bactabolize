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


# Current implementation
* Input model validation
* Assembly QC
* Annotation
    - match with existing ORFs from input reference
        - require 80% overlap of ORF range for a match
    - note differences between annotation bounds in qualifiers
    - transfer qualifiers to new annotation (e.g. locus tag, gene name, gene product)
    - add existing annotations that were unmatched
* Draft model creation
    - Bi-directional BLASTp for model proteins and isolate proteins
        - filter on evalue <= 1e-3, coverage >= 25%, pident >= 80%
    - Discover orthologs
        - defined as a protein pair that are most similar to each other (by pident)
    - Collect model genes that have no ortholog to search at nucleotide level
    - Uni-directional BLASTn for unannotated gene detection
        - filter on evalue <= 1e-3, coverage >= 80%, pident >= 80%
        - require translated sequence to not have a truncating mutation
    - Any BLASTn hit that passes filtering is automatically considered orthologous
        - I don't think this is the best approach
    - Remove models genes that do not have an ortholog in the isolate
        - Artificial genes are excepted here
    - Rename identified orthologs to match locus\_tags in isolate
    - Write model to disk
* Draft model assessment
    - FBA on M9
    - on failure:
        - gapfilling to identify missing genes
        - collate information for debugging
        - exit
* FBA for carbon sources
    - iterate all metabolites that contain carbon
        - need to convince myself on this one


# TODO
* Collate compiled troubleshooting info into a single file
* Output results of draft model assessment FBA to file
* For single-stage execution
    - annotation
        - convert genbank to FASTA, at least in nf
        - add code for annotation matching and transfer
    - model\_fba: create, consider if this is optional or not
        - builtin test sufficient? idk
        - using dummy in nf for now
* Gracefully handle draft model failure in nf
    - either rework draft model code, change exit code
    - or allow specific exit code failure to be ignored in nf
    - alternatively, split channel into pass and into failed
        - add some custom resume function (intended to run after fixing draft issues)
            - e.g. search output directory, minimal FBA check, then proceed to full
* Check required programs are in PATH and have correct versions
* Getting 1182 orthologs on K\_quasi\_quasi\_01A030T, tut notes have 1186
    - 1185 after unannotated gene search (tut detected no annotated genes)
    - resolve differences
* Test BLASTn alignment sequence extraction
    - tested this to some extent and I think it's correct but tut notes have different offsets...


# Planned features
* Faster alternative to BLASTp
    - e.g. diamond
    - will need to demonstrate consistency between results
    - must also ensure that reasonable speed up is obtained
* Provide extensive information for draft genomes that fail to optimise
    - determine essential genes from reference and have these output somewhere
        - required only to be generated once (not per strain)
