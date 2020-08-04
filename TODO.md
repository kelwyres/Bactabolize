# TODO
* Warn once when exchange is missing from model
* Draft output
    - optimisiation value
    - command/code to collate all optimisation codes into single table
        - pass/fail, optimisation value, media
* Patch function
    - flag that a model has failed somewhere
    - for a failed model, require a patch file
    - apply patch
    - run minimal media assessment
    - for now, do not continue to FBA
* For pipeline
    - two stages: drafting and then FBA
        - for now run them separately
        - in future, we could continue to second stage of FBA if none fail
* Investigate lack of translations in CDS featuers for reannotated assemblies
* Clean up file i/o in annotation
    - specifically matching ORFs
* Clean up code in general
* Test fba\_spec\_fp and fba\_types in nf
* Add automated resume/repair command
    - resume with a patch file to apply to model so that it produces biomass on minimal media
* Handle spontaneous genes better
* Check required programs are in PATH and have correct versions
* Getting 1182 orthologs on K\_quasi\_quasi\_01A030T, tut notes have 1186
    - 1185 after unannotated gene search (tut detected no annotated genes)
    - resolve differences
* Test BLASTn alignment sequence extraction
    - tested this to some extent and I think it's correct but tut notes have different offsets...
