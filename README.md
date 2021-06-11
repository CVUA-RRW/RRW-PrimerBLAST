https://zenodo.org/badge/DOI/10.5281/zenodo.4926145.svg

# RRW-PrimerBLAST

This pipeline emulates a Primer-BLAST to find primer matching sequences in a nucleotide 
database and recover the sequences flanked by these matching sites.

With this pipeline you can:
* Filter the query database using an ancestor taxid node
* Set the minimal coverage and identity of the primer sequences for the BLAST step
* Scan the database with an arbitrary number of primers 
* Recover a fasta file with flanked sequences
* Produce a BLAST-formated database of flanked sequences

## Getting started 

### Prerequisites

RRW-PrimerBLAST runs in a UNIX environment with BASH (tested on Debian GNU/Linux 10 (buster)) and requires conda and an internet 
connection (at least for the first run).

### Installing

Start by getting a copy of this repository on your system, either by downloading and unpacking the archive, 
or using 'git clone':

```bash
cd path/to/repo/
git clone --recurse-submodules https://github.com/CVUA-RRW/RRW-PrimerBLAST.git
```

Set up a conda environment containing snakemake, python and the pandas library and activate it:

```bash
conda create --name snakemake -c bioconda -c anaconda snakemake pandas biopython
conda activate snakemake
```

### Getting the databases

RRW-PrimerBLAST requires several databases to run, all are available from the NCBI ftp servers:

* [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)
* [taxdb](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz)
* Any nucleotide collection you want to use, this needs to be a searchable BLAST database with taxonomy information. For this
 you can build a local database from a subset of sequences, for exemple from the BOLD database. Check the
 [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279688/) to know how to do this.

### Running RRW-PrimerBLAST

RRW-PrimerBLAST should be run using the snakemake command-line application.
For this you will need to manually fill the config.yaml file with the paths to the required files.
You can also modify the parameters already present in the file.

Then run the pipeline with:

```bash 
snakemake -s /path/to/FooDMe/Snakefile --configfile path/to/config.yaml --use-conda --conda-prefix path/to/your/conda/envs
```

Consult [snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more details.

### Configuration file

The configuration file contains the following parameters:

```yaml
# Fill in the path belows with your own specifications:
workdir:                    # Path to output directory
blast_db:                   # Path to BLAST-formated database
taxdb:                      # Path to the folder containing the taxdb files
rankedlineage_dmp:          # Path to rankedlineage.dmp
nodes_dmp:                  # Path to nodes.dmp
primers:                    # Path to the primer sequence files (in FASTA), any number of primers is acceptable

# Modify the parameters below:
parent_node: 32524          # Ancestor node to pre-filter BLAST database, 1 to ignore
primerBlast_coverage: 100   # Minimal query coverage for primer BLAST (0-100)
primerBlast_identity: 100   # Minimal sequence identity value for primer BLAST (0-100)
```

## Definition of flanked sequences

Flanked sequences can be seen as in silico PCR amplicons. However note that only
the first amplicon for each database sequence will be returned!
A flanked sequence is simply defined as the first sequences between a primer match 
on the plus strand followed by a match on the minus strand. See examples below:

```
-> indicates a match on the '+' strand 
<- indicates a match on the '-' strand

primer matches:      <- ->  ->    <-  ->
database sequence: =========================
flanked sequence:           ========

primer matches:      -> <-  ->    <-  ->
database sequence: =========================
flanked sequence:    =====       

primer matches:      <- <-  ->    ->  ->
database sequence: =========================
flanked sequence:  
```

## Credits

RRW-PrimerBLAST is built with [Snakemake](https://snakemake.readthedocs.io/en/stable/) and uses the [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 

## Contributing

For new features or to report bugs please submit issues directly on the online repository.

## License

This project is licensed under a BSD 3-Clauses License, see the LICENSE file for details.

## Author

For questions about the pipeline, problems, suggestions or requests, feel free to contact:

Grégoire Denay, Chemisches- und Veterinär-Untersuchungsamt Rhein-Ruhr-Wupper 

<gregoire.denay@cvua-rrw.de>
