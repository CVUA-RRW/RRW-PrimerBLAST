# Generates a BLAST formatted database of predicted amplicons from an existing
# BLAST database and a fasta file containing primers.

import os

shell.executable("bash")


# Settings ------------------------------------------------------------------------------------------------------------------

workdir: config["workdir"]

# Functions -----------------------------------------------------------------------------------------------------------------

def generate_db_name(wildcards=None):
    path, dbname = os.path.split(config["blast_db"])
    path, primername = os.path.split(config["primers"])
    return dbname + '_' + primername.split('.')[0]

# Input rule ----------------------------------------------------------------------------------------------------------------
 
rule all:
    input: 
        # "primer_blast/missing_barcodes.txt",
        "primer_blast/no_barcodes.txt",
        "blast_db/{name}.fasta".format(name = generate_db_name()),
        expand("blast_db/{name}.{ext}", name = generate_db_name(), ext= ["nto", "ntf", "nsq", "not", "nos", "nog", "nin", "nhr", "ndb"]),
        "report.txt"
        
# Workflow ------------------------------------------------------------------------------------------------------------------

rule get_taxid_from_db:
    output:
        temp("db_filtering/taxid_list.txt")
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    message: "Extracting taxid list from database"
    conda: "./envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastdbcmd -db {params.blast_DB} -tax_info -outfmt '%T' > {output}
        """

rule filter_taxid:
    input:
        "db_filtering/taxid_list.txt"
    output:
        mask = "db_filtering/taxid_mask.txt",
        failed = "db_filtering/taxid_missing.txt"
    params:
        taxid = config["parent_node"],
        lineage = config["rankedlineage_dmp"],
        nodes = config["nodes_dmp"]
    message: "Extracting taxids under parent node"
    script:
        "./scripts/make_blast_mask.py"

rule get_seqidlist:
    input:
        "db_filtering/taxid_mask.txt"
    output:
        seqids = temp("db_filtering/seqids.txt"),
        binary = temp("db_filtering/seqids.acc"),
        id_table = temp("db_filtering/table.tsv")
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    message: "Retrieving SeqIDs to search"
    conda: "./envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastdbcmd -db {params.blast_DB} -taxidlist {input} -outfmt '%a\t%T\t%S' > {output.id_table}
        
        cat {output.id_table} | cut -d$'\t' -f1 > {output.seqids}
        
        blastdb_aliastool -seqid_file_in {output.seqids} -seqid_file_out {output.binary}
        """

rule find_primer_matches:
    input:
        seqids = "db_filtering/seqids.txt",
        binary = "db_filtering/seqids.acc",
        primers = config["primers"]
    output:
        "primer_blast/primer_blast.tsv"
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"],
        cov = config["primerBlast_coverage"],
        identity = config["primerBlast_identity"]
    threads: workflow.cores
    message:
        "Blasting primers"
    conda: "./envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastn -db {params.blast_DB} \
            -query {input.primers} \
            -task blastn-short \
            -seqidlist {input.binary} \
            -outfmt '6 saccver qseqid staxid sstart send length sstrand mismatch' \
            -ungapped -qcov_hsp_perc {params.cov} -perc_identity {params.identity} \
            -subject_besthit \
            -max_target_seqs  1000000000 \
            -num_threads {threads} \
                | sort -k1 | sed '1 i\seqid\tquery\ttaxid\tstart\tend\tlength\tstrand\tmismatch' > {output}
        """

rule extract_barcodes_pos:
    input:
        "primer_blast/primer_blast.tsv"
    output:
        "primer_blast/barcode_pos.tsv"
    message: "Extracting barcodes sequences"
    script:
        "./scripts/extract_barcodes.py"

rule extract_barcodes_seq:
    input:
        "primer_blast/barcode_pos.tsv"
    output:
        "blast_db/{name}.fasta"
    message: "Extracting barcode sequences"
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    conda: "./envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        while IFS=$'\t' read -r acc tax start stop length; do
            blastdbcmd -entry $acc \
                -db {params.blast_DB} \
                -range $start-$stop \
                -outfmt %f | sed -e 's/:[[:digit:]-]*//' >> {output}
        done < {input}
        """

rule missing_barcodes: 
    input:
        seqids = "db_filtering/seqids.txt",
        barcodes = "primer_blast/barcode_pos.tsv",
        table = "db_filtering/table.tsv" 
    output:
        acc = temp("primer_blast/no_barcodes.txt"),
        full = "primer_blast/missing_barcodes.txt"
    message: "Identifying missing barcodes"
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    conda: "./envs/blast.yaml"
    shell:
        """
        
        comm -13 <(cat {input.barcodes} | cut -d$'\t' -f1 | sort -k1) \
                <(cat {input.seqids} | sort -k1) \
                | tr -d "\t" \
                > {output.acc}
        
         join --nocheck-order <(sort -b -k1d {input.table}) <(sort -b -k1d {output.acc}) > {output.full}
        """

rule make_barcode_db:
    input: 
        fasta = "blast_db/{name}.fasta".format(name=generate_db_name()),
        id_table = "db_filtering/table.tsv"
    output:
        taxid_mapper = temp("blast_db/taxmap.tsv"),
        DB = expand("blast_db/{name}.{ext}", name = generate_db_name(), 
                    ext= ["nto", "ntf", "nsq", "not", "nos", "nog", "nin", "nhr", "ndb"])
    params: 
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"],
        dbname = generate_db_name
    message: "Formatting barcodes to BLAST database"
    conda: "./envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        cat {input.id_table} | cut -d$'\t' -f1,2 > {output.taxid_mapper}
        
        makeblastdb -in {input.fasta} -dbtype nucl -parse_seqids -blastdb_version 5 -taxid_map {output.taxid_mapper} -out blast_db/{params.dbname}
        """

rule write_report:
    input:
        db_txd = "db_filtering/taxid_list.txt",
        txd_mask = "db_filtering/taxid_mask.txt",
        txd_failed = "db_filtering/taxid_missing.txt",
        primers = config["primers"],
        fasta = "blast_db/{name}.fasta".format(name=generate_db_name()),
        missing = "primer_blast/no_barcodes.txt"
    output:
        "report.txt"
    params:
        blast_DB = config["blast_db"],
        parent = config["parent_node"]
    message: "Logging session info"
    shell:
        """
        date > {output}
        echo "Filtering BLAST database {params.blast_DB}" >> {output}
        echo "    Extracting $(grep -c . {input.db_txd}) taxid" >> {output}
        echo "    $(grep -c . {input.txd_mask}) are descendent of taxid {params.parent}" >> {output}
        echo "    $(grep -c . {input.txd_failed}) were missing from taxdump and ignored" >> {output}
        echo "Blasting primers {input.primers}" >> {output}
        echo "    $(grep -c '^>' {input.fasta}) barcode sequences were found" >> {output}
        echo "    $(grep -c . {input.missing}) sequences did not yield a barcode" >> {output}
        """