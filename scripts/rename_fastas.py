import os
from Bio import SeqIO

new_recs = []
with open(snakemake.input.fasta, "r") as f:
    species = os.path.splitext(os.path.basename(snakemake.input.fasta))[0]
    recs = list(SeqIO.parse(f, "fasta"))
    for i,rec in enumerate(recs):
        some_id_orf, *desc = rec.id.split(" ")
        some_id, orf = some_id_orf.split("_")
        rec.id = "{}_{}".format(species, orf)
        rec.description = ""
        new_recs.append(rec)

with open(snakemake.output.fasta, "w") as f:
    SeqIO.write(new_recs, f, "fasta")