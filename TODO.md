# TODO
* Output results FBA to file
    - draft model assessment
    - media FBA
    - individual rx
* For single-stage execution
    - annotation
        - convert genbank to FASTA, at least in nf
        - add code for annotation matching and transfer
    - model\_fba: create, consider if this is optional or not
        - builtin test sufficient? idk
        - using dummy in nf for now
* Gracefully handle draft model failure in nf
    - allow specific exit code failure to be ignored in nf
    - split channel into pass and into failed
        - add some custom resume function (intended to run after fixing draft issues)
            - e.g. search output directory, minimal FBA check, then proceed to full
* In draft model creation, return all blast results for troubleshooting
* More troubleshooting information
    - for matches to ORFs that did not pass filtering:
        - is it complete
        - has it been reannotated
* Handle spontaneous genes better
* Check required programs are in PATH and have correct versions
* Getting 1182 orthologs on K\_quasi\_quasi\_01A030T, tut notes have 1186
    - 1185 after unannotated gene search (tut detected no annotated genes)
    - resolve differences
* Test BLASTn alignment sequence extraction
    - tested this to some extent and I think it's correct but tut notes have different offsets...
