#! /usr/bin/env python3

"""
Script to produce from a SAM/BAM file, a text with a string of all references
that fall under each taxonomic level specified by the user.

For example, find all reference IDs that fall under each 'species' level.
Strains or sub-species reference IDS are then put in the Species level output
file.

Requires ete3 python taxonomy database in ~/.etetoolkit/taxa.sqlite is already
installed and updated!
"""

## Libraries
import os
import pysam
import ete3
import entrezpy
import argparse
import textwrap
from ete3 import NCBITaxa


from pysam import AlignmentFile

## Custom functions

## function to get unique values
## https://www.geeksforgeeks.org/python-get-unique-values-list/#:~:text=Method%202%20%3A%20Using%20Set,a%20list%20to%20print%20it.

def unique(list1):

    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = list(list_set)
    return unique_list

## Custom Data
mode = {"sam": "r", "bam": "rb", "cram": "rc"}
ranks_basic = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
ranks_order = [
    "superkingdom",
    "kingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "superorder",
    "order",
    "suborder",
    "superfamily",
    "family",
    "subfamily",
    "genus",
    "subgenus",
    "section",
    "species",
    "subspecies",
    "strain",
    "morph",
    "forma specialis",
    "varietas",
]
taxa_to_retain = {}
taxa_to_ignore = {}
collapsed_taxa = {}

## Arguments
parser = argparse.ArgumentParser(
    prog="collapse_sam_taxonomy.py",
    usage="python3 %(prog)s [-h] [--version] -i /<path>/<to>/<input>.{bam/sam} -t <taxonomic_level> -o <output_dir>/",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=textwrap.dedent(
        """\
  Description: %(prog)s produces text files containing BAM or SAM reference IDs
    that group under a user-specific taxonomic level. BAM file should have NCBI
    accession IDs as reference.

  """
    ),
    epilog="Author: James A. Fellows Yates (jfy133[at]gmail.com)",
)

parser.add_argument("-v", "--version", action="version", version="v0.0.1")
parser.add_argument(
    "-i",
    "--input",
    nargs=1,
    required=True,
    help="Sorted SAM/BAM file with NCBI accessions as reference sequence name.",
)
parser.add_argument(
    "-t",
    "--taxonomic_level",
    default="species",
    choices=ranks_basic,
    help="Taxonomic level to group reference sequence names by.",
)
parser.add_argument(
    "-o",
    "--output",
    nargs=1,
    default="results/",
    help="Output directory to place output files.",
)
parser.add_argument(
    "-m",
    "--verbose",
    default=False,
    help="Give verbose warning messages.",
    action="store_true",
)
parser.add_argument(
    "-e",
    "--etetaxdb",
    help="Alternative location of ete3 toolkit taxonomy database if other than default of ~/.etetoolkit/taxa.sqlite"
)

args = parser.parse_args()

## Load NCBI Taxonomy, if nothing supplied assume default location (~/.etetoolkit/taxa.sqlite)
if args.etetaxdb is not None: 
    ncbi = NCBITaxa(dbfile=str(args.etetaxdb))
else:
    ncbi = NCBITaxa()

## Load BAM/SAM
filetype = args.input[0].split(".")[-1]
samfile = pysam.AlignmentFile(args.input[0])

## Parse Limits
min_rank = ranks_order.index(args.taxonomic_level)
refs = unique(list(samfile.references))

if isinstance(args.output, list):
    outdir = args.output[0].rstrip("/")
else:
    outdir = args.output.rstrip("/")

## Make output directory
if not os.path.exists(outdir):
    os.mkdir(outdir)

## Search for ranks
for ref in refs:
    ## Find rank of reference sequence
    if len(ref.split("|")) < 2:
        if args.verbose:
            print("[sam_collapse_taxonomy.py] warn:", ref, "has no tax ID. Ignoring")
        taxa_to_ignore[ref] = ["Reason: no tax ID", {"Reference": ref}]
        continue
    else:
        found_id = int(ref.split("|")[2])
        rank = ncbi.get_rank([found_id])

        if len(rank) == 0:
            if args.verbose:
                print(
                    "[sam_collapse_taxonomy.py] warn: NCBI taxonomy rank for",
                    ref,
                    "not found in database. Ignoring.",
                )
            taxa_to_ignore[ref] = [
                "Reason: no NCBI taxonomy rank",
                {"Reference": ref, "Tax ID": found_id},
            ]
            continue

    ## Check if the rank is a 'valid' one (matching the internal lists) and if 
    ## not exclude
    if rank[found_id] in ranks_order:
        rank_depth = ranks_order.index(rank[found_id])
    else:
        if args.verbose:
            print(
                "[sam_collapse_taxonomy.py] warn: taxid for",
                ref,
                "has ambiguous taxonomic rank. Ignoring.",
            )
        taxa_to_ignore[ref] = [
            "Reason: unrecognised rank",
            {"Reference": ref, "Tax ID": found_id, "Rank": rank},
        ]
        continue

    ## Discard ranks higher in hierarchy than users' requested
    if rank_depth >= min_rank:
        taxa_to_retain[ref] = [found_id, rank[found_id]]
    else:
        if args.verbose:
            print(
                "[sam_collapse_taxonomy.py] warn: taxid for",
                ref,
                "rank higher in hierarchy than max specified! Ignoring!",
            )
        taxa_to_ignore[ref] = [
            "Reason: rank higher in hierarchy than specified",
            {"Reference": ref, "Tax ID": found_id, "Rank": rank},
        ]

## Collapse rank to user specificed level
for ref in taxa_to_retain.keys():
    ## for given taxon ID, find level rank in hierarchy
    rank = taxa_to_retain[ref][1]
    id = taxa_to_retain[ref][0]
    tax_level = ranks_order.index(rank)

    ## if equal or higher retain find equivalent taxa level requested
    ## by user otherwise discard as uninformative
    if tax_level >= min_rank:
        lineage = ncbi.get_rank(ncbi.get_lineage(id))

        if args.taxonomic_level not in lineage.values():
            if args.verbose:
                print(
                    "[sam_collapse_taxonomy.py] warn: taxid for",
                    ref,
                    "doesn't have the specified rank in it's NCBI taxonomy lineage (e.g. unspecifiedSphingomonadales)! Ignoring",
                )
            taxa_to_ignore[ref] = [
                "Reason: specified rank does not exist for this taxid",
                {"Reference": ref, "Tax ID": id, "Rank": rank},
            ]
            continue
        else:
            tax_limit_pos = list(lineage.values()).index(args.taxonomic_level)

        result_taxid = list(lineage)[tax_limit_pos]
        result_name = str(list(ncbi.get_taxid_translator([result_taxid]).values())[0])
        if str(result_name) not in collapsed_taxa:
            collapsed_taxa[str(result_name)] = [ref]
        else:
            collapsed_taxa[str(result_name)].append(ref)
    else:
        if args.verbose:
            print(
                "[sam_collapse_taxonomy.py] warn: unable to get specified rank for taxid ",
                ref,
                "! Ignoring!",
            )
        taxa_to_ignore[ref] = [
            "Reason: rank lower than specified",
            {"Reference": ref, "Tax ID": id, "Rank": rank},
        ]
        taxa_to_retain.pop(ref)

## Write failed collapsing
with open(outdir + "/unresolved_taxonomic_ids.txt", "w") as f:
    for key, value in taxa_to_ignore.items():
        f.write("%s:%s\n" % (key, value))

## Write collapsable taxonomic IDs, one per feile
for key in collapsed_taxa:
    with open(
        outdir + "/%s.txt" % str(key.replace(" ", "_").replace("/", "_")), "w"
    ) as f:
        for name in collapsed_taxa[key]:
            f.write("%s " % name)

