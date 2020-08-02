# TODO
* For FBA, allow specification of minimal media
    - behaviour may need to change for supplying default element sources
        - e.g. if a sulfur source is in the media, use that?
        - might be best to allow explicity specification of default element sources
* Check for netmeta metanetx identifiers
    - requried for exchange source assignment in FBA (processes compund forumla)
* Investigate lack of translations in CDS featuers for reannotated assemblies
* Clean up file i/o in annotation
    - specifically matching ORFs
* Clean up code in general
* Report failed draft models in nf somewhere
* Test fba\_spec\_fp and fba\_types in nf
* Add automated resume/repair command
    - best approach would seem to allow user to specify what changes to make via a file
        - this would be done after user has completed manual investigation and found a solution
    - enable reproducible workflow - save file for rerunning pipeline, etc
    - easy for user to set what needs to be fixed and easy to apply in code
    - would need to integrate into NF somehow
        - config option to search a path for incomplete models?
* Handle spontaneous genes better
* Check required programs are in PATH and have correct versions
* Getting 1182 orthologs on K\_quasi\_quasi\_01A030T, tut notes have 1186
    - 1185 after unannotated gene search (tut detected no annotated genes)
    - resolve differences
* Test BLASTn alignment sequence extraction
    - tested this to some extent and I think it's correct but tut notes have different offsets...
