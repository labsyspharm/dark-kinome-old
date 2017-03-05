# DarkKinome

The associated files are available on Synapse at: https://www.synapse.org/#!Synapse:syn8377973/files/

## Prerequisites
The scripts utilize the following R libraries: dplyr, synapseClient.

## Wrangling the summary table
The `summary.tsv` table can be produced by running the following on the command line: `Rscript wrangle.R`

The Overall Score is computed according to the following criteria: `s1 + 2*s2 + s3 + s4 + 5*s5`, where
- `s1 = 1/4 ( I(Recombinant Kinase) + I(3D Structure) + I(Antibody) + I(ThermoFisher)`
- `s2 = 1/2 ( I(Disease) + Normalized Genomic Aberration Score)`
- `s3 = 1/2 ( I(Small Molecule Ligand) + Normalized Specificity Score)`
- `s4 = Normalize Biological Activity Score`
- `s5 = I( Expressed in LSP datasets )`

