## Clone this repository
`mkdir ~/Experiments
cd ~/Experiments/
gh repo clone NolanBentley/Tspicata`

## Create a directory called ./data/data_ignored/ 
This is where you will store large datasets that can't be uploaded to github

## Add the following datasets to this directory:
T_spicatus_V2.2.fasta
Tspicata_both.d8.merged.vcf
TS_all_filtered.recode_wide1.tsv

## To generate fasta file for MEGA run:
`cd ~/Experiments/Tspicata #Or a modified path`
`bash ./scripts/tsvToFasta.R`

## To generate indel subset of vcf run:
`cd ~/Experiments/Tspicata #Or a modified path`
`rscript ./scripts/subsettingToIndels.R`

## To subset indels to gel-amendable markers run:
`cd ~/Experiments/Tspicata #Or a modified path`
`rscript ./scripts/findingIndelMarkers.R`

