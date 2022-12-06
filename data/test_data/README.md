# Test data

This directory contains reproducible test data for users. The [working_model](https://github.com/kelwyres/Bactabolize/tree/main/data/test_data/working_model) contains the annotated K. pneumoniae KPPR1 genome which will successfully produce a model when the [KpSC_pan model](https://github.com/kelwyres/KpSC-pan-metabolic-model) is used as a reference, along with expected outputs.
The [model_patching_example](https://github.com/kelwyres/Bactabolize/tree/main/data/test_data/model_patching_example) contains the annotated K. africana KPN2285 genome  will not produce a working model when the [KpSC_pan model](https://github.com/kelwyres/KpSC-pan-metabolic-model) is used as a reference. This can be used as an example for the [patch_model](https://github.com/kelwyres/Bactabolize/wiki/6.-Trouble-shooting-models) command, along with expected outputs.


## draft_model

In-depth instructions found [here](https://github.com/kelwyres/Bactabolize/wiki/3.-Build-a-model)
- KpSC reference model found here: [KpSC_pan-metabolic_model_v1.json](https://github.com/kelwyres/KpSC-pan-metabolic-model/blob/main/KpSC_pan-metabolic_model_v1.json)
- KpSC reference gene data found here: [KpSC_pan-metabolic_model_v1_nucl.fna](https://github.com/kelwyres/KpSC-pan-metabolic-model/blob/main/KpSC_pan-metabolic_model_v1_nucl.fna)
- KpSC reference protein data found here: [KpSC_pan-metabolic_model_v1_prots.faa](https://github.com/kelwyres/KpSC-pan-metabolic-model/blob/main/KpSC_pan-metabolic_model_v1_prots.faa)
- KPN2285 patch file found here: [KPN2285_patch.json](https://github.com/kelwyres/Bactabolize/blob/main/data/test_data/model_patching_example/KPN2285_patch.json)

Successful model building test data example:
```
bactabolize draft_model --assembly_fp data/working_model/K.pneumoniae_KPPR1.gbk --ref_genes_fp KpSC_pan-metabolic_model_v1_nucl.fna --ref_proteins_fp KpSC_pan-metabolic_model_v1_prots.faa --ref_model_fp KpSC_pan-metabolic_model_v1.json --biomass_reaction_id BIOMASS_Core_Oct2019 --output_fp K.pneumoniae_KPPR1
```

Non-successful model building test data, requiring patching example:
```
bactabolize draft_model --assembly_fp data/model_patching_example/KPN2285.gbk --ref_genes_fp KpSC_pan-metabolic_model_v1_nucl.fna --ref_proteins_fp KpSC_pan-metabolic_model_v1_prots.faa --ref_model_fp KpSC_pan-metabolic_model_v1.json --biomass_reaction_id BIOMASS_Core_Oct2019 --output_fp KPN2285
```


## patch_model for non-functioning draft models

In-depth instructions found [here](https://github.com/kelwyres/Bactabolize/wiki/6.-Trouble-shooting-models)

Successful model patching example:

1. Create a `KPN2285_patch.json` file using the missing reactions from the output `*troubleshoot_summary.txt`. Use [KPN2285_patch.json](https://github.com/kelwyres/Bactabolize/blob/main/data/test_data/model_patching_example/KPN2285_patch.json) as template
2. Run patch_model

```
bactabolize patch_model --draft_model_fp data/model_patching_example/KPN2285_model.json --ref_model_fp KpSC_pan-metabolic_model_v1.json --patch_fp KPN2285_patch.json --biomass_reaction_id BIOMASS_Core_Oct2019 --output_fp KPN2285_patched.json
```
