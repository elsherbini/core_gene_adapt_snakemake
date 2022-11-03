library(Biostrings)
library(tidyverse)

forward_primer <- DNAString("ACCCDATGGATCTKAAGCGYGG")
reverse_primer <- DNAString("CCRCGRTCGAACTGCATVCCTT")

bloomers <- tribble(~genotype, ~amplicon_sequence,
"cycl_bloomer", "TATCGACAAAGCGGTTATCGCGGCTGTTGAAGAACTGAAGAACCTTTCTGTTCCATGTGCAGATACGAAAGCTATCGCGCAAGTAGGTACTATCTCTGCAAACTCTGACGCAACAGTGGGCAACATTATTGCTGAAGCGATGGAAAAAGTAGGCCGTGATGGCGTTATCACTGTTGAAGAAGGTCAGGCTCTACAAGACGAGCTAGACGTAGTTG",
"splen_bloomer", "CATCGACAAAGCGGTTATCGCGGCTGTTGAAGAGCTGAAGAACCTTTCTGTTCCTTGTTCAGACACGAAAGCTATCGCGCAAGTAGGTACTATCTCTGCGAACTCTGATTCAACAGTAGGTAACATCATTGCTGAAGCGATGGAAAAAGTAGGTCGTGATGGTGTAATCACGGTTGAAGAAGGTCAGGCTCTGCAAGACGAGCTAGACGTAGTTG")

orfs <- read_csv("output/orf_table.csv")
hmmer_results <- read_csv("output/hmmer_cpn60/cpn60_hmmsearch.tbl.csv")

get_amplicon_sequence <- function(nt_seq) {
  forward_match <- (pairwiseAlignment(pattern=forward_primer, subject=DNAString(nt_seq), type="overlap"))@subject@range
  reverse_match <- (pairwiseAlignment(pattern=reverse_primer, subject=reverseComplement(DNAString(nt_seq)), type="overlap"))@subject@range
  forward_pos <- forward_match@start + forward_match@width
  reverse_pos <- -1*(reverse_match@start + reverse_match@width) + str_length(nt_seq) + 1
  str_sub(nt_seq, forward_pos, reverse_pos)
}

hmmer_results %>%
  filter(full_E < 0.1) %>%
  separate(target, into=c("strain"), sep="_", extra="drop", remove=FALSE) %>%
  group_by(strain) %>%
  arrange(full_E) %>%
  mutate(hsp60_match=n()) %>%
  ungroup() %>%
  left_join(select(orfs, target=protein_accession, orf_seq, fasta_id, contig, orf)) %>%
  mutate(amplicon_sequence = map_chr(orf_seq, get_amplicon_sequence)) %>%
  mutate(amplicon_length = str_length(amplicon_sequence)) %>%
  select(strain, hsp60_match, full_E, amplicon_length, amplicon_sequence, target, contig, orf, everything()) %>%
  left_join(bloomers) %>%
  write_csv(snakemake@output[["amplicon_sequences"]])
