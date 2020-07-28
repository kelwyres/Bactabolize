# TODO
* Investigate lack of translations in CDS featuers for reannotated assemblies
* Clean up file i/o in annotation
    - specifically matching ORFs
* Gracefully handle draft model failure in nf
    - allow specific exit code failure to be ignored in nf
    - split channel into pass and into failed
        - add some custom resume function (intended to run after fixing draft issues)
            - e.g. search output directory, minimal FBA check, then proceed to full
* Add automated resume/repair command
    - best approach would seem to allow user to specify what changes to make via a file
        - this would be done after user has completed manual investigation and found a solution
    - enable reproducible workflow - save file for rerunning pipeline, etc
    - easy for user to set what needs to be fixed and easy to apply in code
* Handle spontaneous genes better
* Check required programs are in PATH and have correct versions
* Getting 1182 orthologs on K\_quasi\_quasi\_01A030T, tut notes have 1186
    - 1185 after unannotated gene search (tut detected no annotated genes)
    - resolve differences
* Test BLASTn alignment sequence extraction
    - tested this to some extent and I think it's correct but tut notes have different offsets...
