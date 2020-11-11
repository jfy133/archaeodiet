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

## TODO
## Add input BAM exists verification

#### Libraries #################################################################

import os
import pysam
import ete3
import entrezpy
import argparse
import textwrap
import time
import datetime

## import specific functions/objects
from tqdm import tqdm, trange
from ete3 import NCBITaxa
from Bio import Entrez
from pysam import AlignmentFile

#### Custom functions ##########################################################

## function to get unique values: https://www.geeksforgeeks.org/python-get-unique-values-list/#:~:text=Method%202%20%3A%20Using%20Set,a%20list%20to%20print%20it.

def unique(list1):

    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = list(list_set)
    return unique_list

## Split list into X chunks: By oremj from https://stackoverflow.com/a/1751478
def chunk(l, n):
    n = max(1, n)
    return (l[i:i+n] for i in range(0, len(l), n))

## Small function for using join with map
def comma_join(x):
    return ','.join(x)

## Pull taxid from efetch results
def extract_taxids_from_efetch(x):
    res = {}
    acc = x['Caption']
    taxid = x['TaxId']
    res[acc] = taxid
    return res

#### Custom Data ###############################################################

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
    "serogroup",
    "serotype",
    "morph",
    "forma specialis",
    "varietas",
    "isolate",
    "biotype",
    "no rank"
]

ranks_exclude = [
    "serogroup",
    "serotype",
    "isolate",
    "no rank",
    "biotype"
]

#### User Arguments ############################################################

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
    "--outputdir",
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
    "-d",
    "--etetaxdb",
    help="Alternative location of ete3 toolkit taxonomy database if other than default of ~/.etetoolkit/taxa.sqlite"
)
parser.add_argument(
    "-a",
    "--entrez_apikey",
    nargs=1,
    required=True,
    help="Your Entrez API key to allow fast querying. See: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/"
)
parser.add_argument(
    "-e",
    "--entrez_email",
    nargs=1,
    required=True,
    help="Email address associated with your API key."
)

args = parser.parse_args()

#### Variable Preparation ######################################################

## Load NCBI Taxonomy, if nothing supplied assume default location (~/.etetoolkit/taxa.sqlite)
if args.etetaxdb is not None: 
    ncbi = NCBITaxa(dbfile=str(args.etetaxdb))
else:
    ncbi = NCBITaxa()

## Prepare Entrez variables
Entrez.api_key = args.entrez_apikey
Entrez.email = args.entrez_email
Entrez.max_tries = 5
Entrez.sleep_between_tries = 30

## Parse Limits
min_rank = ranks_order.index(args.taxonomic_level)

## Make output directory
if isinstance(args.outputdir, list):
    outdir = args.outputdir[0].rstrip("/")
else:
    outdir = args.outputdir.rstrip("/")

if not os.path.exists(outdir):
    os.mkdir(outdir)

## Required internal pre-prepped objects
query_dict = {}
search_results = list()
queries_flattened = {}
queries_valid = {}
result_list = {}

#### START #####################################################################

## Load BAM/SAM
filetype = args.input[0].split(".")[-1]
samfile = pysam.AlignmentFile(args.input[0])

## Extract reference sequence IDs
refs = unique(list(samfile.references))

for ref in refs:
    query_dict[ref.split('|')[0].split('.')[0]]=ref

if args.verbose:
    print("[collapse_taxonomy.py] info: Detected", len(refs) , "reference sequences")

## Clean up refs and combine to single list
refs = [i.split('|')[0] for i in refs] 
','.join(refs)

chunked_res = list(chunk(refs, 1000))
search_queries = map(comma_join, chunked_res)

#### Entrez search #############################################################

if args.verbose:
    print("[collapse_taxonomy.py] info: Generated", len(chunked_res), "entrez search queries")

n = 1

for i in tqdm(search_queries, desc="[collapse_taxonomy.py] info: submitting entrez queries", unit=" queries", ascii=True, total=len(chunked_res)):
    handle = Entrez.efetch(db = 'nucleotide', id=i, rettype="docsum")
    res = Entrez.read(handle)
    search_results.append(res)
    handle.close()
    n += 1


## Flatten and clear
search_results_flat = [y for x in search_results for y in x]

if args.verbose:
    print("[collapse_taxonomy.py] info: Received", len(list(search_results_flat)), "entrez search query results")

queries_flattened = list(map(extract_taxids_from_efetch, search_results_flat))

for i in queries_flattened:
    queries_valid.update(i)

## Add results back with ref dictionary
for k in query_dict:
    query_dict[k] = {'ref': query_dict[k], 'taxid': queries_valid.get(k, {})}

#### ete3 search ###############################################################

## Find the tax level of tax_id and the corresponding rank number
for query in query_dict:
    if query_dict[query]['taxid'] == {}:
        query_dict[query] = {'ref': query_dict[query], 'taxid': 'no_tax_id', 'rank': "none", 'rank_level': -1}
        continue
    
    rank = ncbi.get_rank([query_dict[query]['taxid']])
    
    if len(rank) == 0:
        if args.verbose:
            print("[collapse_taxonomy.py] warn: no taxonomic rank found for:", query)
        query_dict[query] = {'ref': query_dict[query], 'taxid': queries_valid.get(query, {}), 'rank': "none", 'rank_level': -1}
    else:
        tax_level = list(rank.values())[0]
        rank_level = ranks_order.index(tax_level)
        if args.verbose:
            print("[collapse_taxonomy.py] info: taxonomic rank for", query, "is", tax_level)
        query_dict[query] = {'ref': query_dict[query]['ref'], 'taxid': queries_valid.get(query, {}), 'rank': rank.values(), 'rank_level': rank_level}

## Evaluate whether query can be collapsed up to the request taxonomic level
## Translate user specification so we can find the number in our rank order
min_rank = ranks_order.index(args.taxonomic_level)

for query in query_dict:
    taxid = query_dict[query]['taxid']
    rank_level = query_dict[query]['rank_level']

    if rank_level <= min_rank:
        if args.verbose:
            print("[collapse_taxonomy.py] warn: taxonomic level too high or unknown for:", query)
        query_dict[query]['collapsed_rank'] = "too_high_or_unknown"
        continue

    else:

        lineage = ncbi.get_rank(ncbi.get_lineage(taxid))
        
        if args.taxonomic_level not in lineage.values():
            if args.verbose:
                print("[collapse_taxonomy.py] warn: requested rank does not exist for:", query)
            query_dict[query]['collapsed_rank'] = "assigned_rank_noname"
            continue
        else:
            ## Look up the requested taxonomic level in lineage, and find corresponding name from taxid
            if args.verbose:
                print("[collapse_taxonomy.py] info: requested rank has been found for:", query)
            tax_limit_pos = list(lineage.values()).index(args.taxonomic_level)
            result_taxid = list(lineage)[tax_limit_pos]
            result_name = str(list(ncbi.get_taxid_translator([result_taxid]).values())[0])

            ## Save
            query_dict[query]['collapsed_rank'] = result_name

#### Output preparation ########################################################

## Group all reference accession IDs together per found taxonomic level
for i in query_dict:
    
    taxon_name = query_dict[i]['collapsed_rank']
    ref_name = query_dict[i]['ref']
    
    print
    
    if taxon_name in list(result_list.keys()):
        result_list[taxon_name].append(ref_name)
    else:
        result_list[taxon_name] = list()
        result_list[taxon_name].append(ref_name)

if args.verbose:
    print("[collapse_taxonomy.py] info: saving output files for", len(result_list), "taxonomic levels")


for key in result_list:
    with open(
        outdir + "/%s.txt" % str(key.replace(" ", "_").replace("/", "_")), "w"
    ) as f:
        for name in result_list[key]:
            f.write("%s " % name)
