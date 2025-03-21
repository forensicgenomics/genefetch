# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
global_defaults.py

Global variables used during the fetching process. Imported by mutliple modules.
Adjust these in this file to complement your process.

The vars:
   MAX_NUM, BATCH_SIZE, FETCH_PARALLEL, NUM_WORKERS, SOFT_RESTART, SEARCH_TERM
are also settable as a command line argument when running fetch.py

Author: Noah Hurmer as part of the mitoTree Project.
"""

from datetime import datetime
import os


####
# output dir
if not os.path.exists("data"):
    os.makedirs("data")
####
# output file locations
IDS_FILE = "data/ids_list.txt"
METADATA_FILE = "data/metadata.csv"
REMOVED_IDS_FILE = 'data/removed_ids.csv'
####
# fasta file locations
SEQS_DIR = "data/seqs"
if not os.path.exists(SEQS_DIR):
    os.makedirs(SEQS_DIR)
####
# current time for writing log and processed files
RUN_TIME = datetime.now()
TIMESTAMP = RUN_TIME.strftime("%d_%m_%Y_%H_%M_%S")
####
#################################################################################
####
# vars relevant to the fetching process
#
# term to eSearch the ncbi database with
# see the ncbi nuccore docu on how to adjust this
SEARCH_TERM = "mitochondrion[filter] AND Homo Sapiens[Organism] AND 16500:99999999[SLEN] NOT Homo sapiens neanderthalensis"
# whether to use threading to fetch in parallel and with how many threads
# this speeds up the fetching by a lot as the Entrez response can be rather slow
FETCH_PARALLEL = False
NUM_WORKERS = 16
# how many profiles to fetch
MAX_NUM = 500
# batch size to save them in
BATCH_SIZE = 100
# LIMIT_NUM is set to return all search results
#   (used to constrain the max value for soft restarting
#   and for the post-process check to check if any profiles were not downloaded correctly)
LIMIT_NUM = 100000
# fasta check flag
# if set, the fetcher not also checks for the existenz of a profile in the 'ids_list.txt'
# but also if the fasta file for that profile is actually downloaded.
CHECK_FOR_FASTA = True
####
#################################################################################
####
# soft restart
#
# soft restart flag
# If set, skips profiles in the latest 'processed_ids' file
# useful when fetching a lot of data manually in smaller batches without needing to raise MAX_NUM
SOFT_RESTART = False
# processed ids file location
PROCESSED_IDS_DIR = "data/processed_ids"
if not os.path.exists(PROCESSED_IDS_DIR):
    os.makedirs(PROCESSED_IDS_DIR)
CURRENT_PROCESSED_IDS_FILE = os.path.join(PROCESSED_IDS_DIR, f"processed_ids_{TIMESTAMP}.txt")
####
# flag to remove profile sequence, metadata etc, if not within current call fetch.
CLEAN_DIR = False
#################################################################################
####
# filters which profiles to remove while fetching
# current filters to remove:
#   - length <16500 ; (In theory unnecessary, as this can be dealt with in the search string)
#   - more than 30 missing bases in the sequence
#   - more than 30 ambiguous bases in the sequence
#   - human filter ; (necessary as the search string also returns 'Homo sapiens' subspecies)
#   - manual exclusions are added dynamically
#
# adjust / delete / append these filter list entries as needed
FILTERS = [
    {
        'description': '<16519',
        'fun': lambda record: len(record.seq) < 16519
    },
    {
        'description': '>30 N/n',
        'fun': lambda record: record.seq.count('N') + record.seq.count('n') > 30
    },
    {
        'description': '>30 amb',
        'fun': lambda record: sum(1 for base in record.seq if base not in 'ACTGactg') > 30
    },
    {
        'description': 'species',
        'fun': lambda record: record.annotations['organism'] != "Homo sapiens"
    }
]
# manual exclusions are added dynamically with all files in the following dir
EXCLUSIONS_DIR = "exclusions"
####
#################################################################################
####
# secrets locations
API_DEST = "secrets/ncbi_api_key.txt"
EMAIL_DEST = "secrets/ncbi_email.txt"
####
#################################################################################
####
# logging setup
# produces logfile of two different levels datestamped in the log dir
LOG_DIR = "data/logs"
if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)
LOG_FILE = f"{LOG_DIR}/fetchlog_{TIMESTAMP}.log"
####
# debug directory
DEBUG_DIR = os.path.join(LOG_DIR, "debug_info")
if not os.path.exists(DEBUG_DIR):
    os.makedirs(DEBUG_DIR)
DEBUG_LOG_FILE = f"{DEBUG_DIR}/fetchlog_{TIMESTAMP}_DEBUG.log"
####
#################################################################################
####
# file to save date of the last run to
# if present the date is used to fetch metadata for changed entries between then and now
LAST_RUN_PATH = "data/last_run_date.txt"
####
#################################################################################
####
# provides the template dictionary for the metadata to fetch.
# caution: adjusting this does not automatically include or exclude the values!
METADATA_TEMPLATE = {
    "accession": None,
    "pub_title": None,
    "first_aut": None,
    "pubmed_id": None,
    "pub_date": None,
    "geo_origin": None,
    "asm_method": None,
    "seq_tech": None
}
####
