# Queries
* Genes can be considered orthologous from only 25% of the protein sequence
    - currently filtering alignments with <25% coverage
    - if a pair of proteins have the highest PID to each other, then they're orthologs
    - not sure if we want to increase coverage requirement for this
* When capturing unannotated genes with BLASTn, we're only doing a one-way alignment
    - we probably should be taking these hits and align them against the reference
        - aligning translated sequence is most appropriate here I would think
    - then checking if they qualify as orthologs by the standard definition
* For full FBA, selection of extracellular metabolites are done by iterating reactions that exchange mass with the exterior
  compartment (generally extracellular)
    - this is a method of the cobra.model but does not seem to include reactions occuring in both extracellular and periplasm
        - this is mostly diffusion reactions
        - from the reaction annotation these do not appear to be exchange reactions
    - methods in Shigella paper with John suggests only exchange reactions are considered, not these diffusion reactions


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
    - iterate all metabolites that contain carbon, sulfur, nitrogen, etc
    - FBA on specific media set
    - FBA on user-provided spec (JSON format)

# Planned features
* Faster alternative to BLASTp
    - e.g. diamond
    - will need to demonstrate consistency between results
    - must also ensure that reasonable speed up is obtained
