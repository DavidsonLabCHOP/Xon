## Intropolis README

### Purpose:

These intropolis database search tools are being released as a part of Monteys et al., Nature (2021). They are used to query splice junctions of interest in the Intropolis database of exon-exon junctions to identify the datasets in which a junction has been previously reported.

### Contents:

- intropolis_search.R
- ExtractIntropolisLine.py

### Use:

The first script `intropolis_search.R` is used to extract the intropolis entry for a splice junction of interest. This process is started by obtaining a copy of the intropolis database containing only the first 4 columns (rowNumber, chromosome, startPosition, endPosition).  Call this file "smallIntropDB.tsv".  The `intropolis_search.R` script is then used to identify the line of the intropolis database represinting a splice event of interest.  Once this line number is obtained we then use a second python-based script "ExtractIntropolisLine.py" to extract (by rowNumber) the number and identity of datasets that contained that splice event from the full sized intropolis database.

### Intropolis:

Intropolis contains a list of exon-exon junctions found across 21,504 RNA-Seq Datasets. It can be accessed here https://github.com/nellore/intropolis
