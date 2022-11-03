import pandas as pd
from collections import namedtuple
import os

INPUT_TABLE = snakemake.input.table
OUTPUT_CSV = snakemake.output.csv
POPULATION = os.path.splitext(os.path.basename(snakemake.input.table))[0]

def parse_info(info):
    info = info.split()
    return HitInfo(*info)._asdict()

hit_info_fields = (
    'target, target_accession, query, query_accession,'
    'full_E, full_score, full_bias,'
    'best_domain_E, best_domain_score, best_domain_bias,'
    'exp, num_discrete_region, num_multidomain,'
    'num_overlap, num_envelope, num_defined_domain,'
    'num_domain_reporting_threshold, num_domain_inclusion_threshold, description')
HitInfo = namedtuple("HitInfo", hit_info_fields)


def parse_description(description):
    orf_start, orf_end, strand = description
    orf_length = abs(int(orf_start) - int(orf_end))
    return Description(orf_start, orf_end, orf_length, strand)._asdict()
description_fields = 'orf_start, orf_end, orf_length, strand'
Description = namedtuple("Description", description_fields)


def parse_prodigal_info(prodigal_info):
    (prod_id, prod_partial,
        prod_start_type, prod_rbs_motif, prod_rbs_spacer,
        prod_gc_cont) = [field.split("=")[1] for field in prodigal_info.split(";")]
    prod_partial_left, prod_partial_right = prod_partial
    return ProdInfo(prod_id, prod_partial_left, prod_partial_right,
                    prod_start_type, prod_rbs_motif, prod_rbs_spacer,
                    prod_gc_cont.strip())._asdict()
prodigal_info_fields = (
    'prod_id, prod_partial_left, prod_partial_right,'
    'prod_start_type, prod_rbs_motif, prod_rbs_spacer, prod_gc_cont')
ProdInfo = namedtuple("ProdInfo", prodigal_info_fields)


def parse_tblout(tblout):
    with open(tblout, "r") as f:
        for line in f:
            if line[0] == "#":
                next
            else:
                yield(parse_hit(line))


def parse_hit(line):
    #info, *description, prodigal_info = line.split("#")
    info = line
    parsed_info = parse_info(info)
    #parsed_description = parse_description(description)
    #parsed_prodigal_info = parse_prodigal_info(prodigal_info)
    #return {'population': POPULATION, **parsed_info, **parsed_description, **parsed_prodigal_info}
    return {'population': POPULATION, **parsed_info}


df = pd.DataFrame(parse_tblout(INPUT_TABLE))
if not df.empty:
    df.to_csv(OUTPUT_CSV, index=False)
else:
    with open(OUTPUT_CSV, "w") as f:
        f.write("target,hmm_name\n")