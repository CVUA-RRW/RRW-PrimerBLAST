from Bio import Seq, SeqIO
from itertools import product

# https://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence
# Jivan's answer
def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   d = Seq.IUPAC.IUPACData.ambiguous_dna_values
   return list(map("".join, product(*map(d.get, seq))))

def primers_to_fasta(name, seq_list):
    """return fasta string of primers with tracing newline"""
    fas = ""
    for i in range(len(seq_list)):
        fas += f">{name}[{i}]\n{seq_list[i]}\n"
    return fas

def main(fastain, fastaout):
    with open(fastain, "r") as fin, open(fastaout, "w") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            explicit = extend_ambiguous_dna(record.seq)
            fasta = primers_to_fasta(record.id, explicit)
            fout.write(fasta)

if __name__ == '__main__':
    main(snakemake.input[0], snakemake.output[0])
    