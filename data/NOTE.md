# Formula data
Very few metabolites in the BiGG Kp model have a chemical formula. To correctly categorise metabolites by source type (e.g.
carbon, sulfur, phosphate, etc) I have downloaded chemical formulas from metanetx, which BiGG metabolites externally link to.

```bash
# Formulas (~500MB). To reduce this file size I have selected only the identifier and formula column, then compressed the
result
wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv

# Deprecated and new identifiers; many MNX identifiers in the BiGG models are deprecated and need to be resolved
wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_depr.tsv
```
