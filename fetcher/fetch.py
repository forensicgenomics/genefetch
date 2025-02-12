# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
fetch.py

Main script for fetching mitochondrial genome data from NCBI, processing it,
and saving results. Supports parallel and sequential fetching, soft restarts,
and batch processing.

Usage:
    Run this script directly with the location as working directory, or
    execute this script via command line using:
    `
        python fetch.py --search-term "mitochondrion complete genome AND Homo Sapiens[Organism]"
            --max-num 1000 --batch-size 100 --fetch-parallel --num-workers 20 --soft-restart
    `

    Adjust other variables within the `global_defaults.py` script if needed.

Dependencies: argparse, random, HTTPError, BioPython, pandas, concurrent.futures, ratelimit,
 custom modules for fetching, file I/O, and metadata processing.

Author: Noah Hurmer as part of the mitoTree Project.
"""

import argparse
import random
import warnings
from io import StringIO
from urllib.error import HTTPError
from Bio import SeqIO, Entrez
import concurrent.futures
import threading
import time
from ratelimit import limits, sleep_and_retry

from .fetch_tools import (set_entrez_rate,
                         set_entrez_globals,
                         filter_changed_profiles,
                         fetch_profile_accs,
                         readd_recently_modified_profiles)
from .file_io import (save_batch_info,
                      filter_unprocessed_ids,
                      load_processed_ids,
                      load_local_versions,
                      load_removed_versions,
                      cleanup_old_files,
                      get_last_run_date,
                      write_last_run_date,
                      write_seq_as_fasta,
                      post_process_metadata,
                      clean_profiles_from_data)
from .metadata_tools import (get_pubmed_info,
                            get_assembly_info,
                            get_geo_info)
from .global_defaults import (LIMIT_NUM,
                              SEARCH_TERM,
                              PROCESSED_IDS_DIR,
                              LOG_DIR,
                              RUN_TIME,
                              NUM_WORKERS,
                              MAX_NUM,
                              BATCH_SIZE,
                              SOFT_RESTART,
                              FETCH_PARALLEL,
                              METADATA_TEMPLATE,
                              CLEAN_DIR,
                              DEBUG_DIR)
from .filter_tools import load_filters
from .logger_setup import get_logger
from .post_process_check import main as post_process_check


####
# global vars
VERBOSE = False

# set up logging
logger = get_logger()
logger.info(f"For a debug level log file, view the logfile within {DEBUG_DIR}.\n")

# file lock writing for parallel fetching
write_lock = threading.Lock()

# dynamically load all filters
FILTERS = load_filters()

### all of this needs to stay here otherwise it causes some parallel problems
set_entrez_globals()
MAX_CALLS_PER_SECOND = set_entrez_rate()
rate_limit_semaphore = threading.Semaphore(MAX_CALLS_PER_SECOND)

def wait_helper(attempt_num):
    sleep_time = 2 ** attempt_num + random.uniform(0, 1)
    time.sleep(sleep_time)

@sleep_and_retry
@limits(calls=MAX_CALLS_PER_SECOND, period=1)
def rate_limited_call(api_call, *args, **kwargs):
    """
    Rate-limited call wrapper with retries

    Args:
        api_call (callable): Entrez API function to call.
        *args, **kwargs: Arguments to pass to the API call.

    Returns:
        object: The result of the API call.
    """

    with rate_limit_semaphore:
        time.sleep(1 / MAX_CALLS_PER_SECOND)
        for attempt in range(3):
            try:
                return api_call(*args, **kwargs)
            except HTTPError as e:
                # TODO maybe we should just blanket retry all HTTP errors and retry regardeless of type?
                if e.code == 429:  # Too Many Requests
                    if logger:
                        logger.warning(f"HTTP error: {e} ; Retrying. Attempt {attempt + 1}/3.")
                    wait_helper(attempt)
                    continue
                elif e.code == 400:  # Bad Request
                    if logger:
                        logger.warning(f"HTTP error: {e} ; Retrying. Attempt {attempt + 1}/3.")
                    wait_helper(attempt)
                    continue
                elif e.code == 500:  # Internal Server Error
                    if logger:
                        logger.warning(f"HTTP error: {e} ; Retrying. Attempt {attempt + 1}/3.")
                    wait_helper(attempt)
                    continue
                else:
                    if logger:
                        logger.error(f"HTTP error: {e}")
                    raise
            except Exception as e:
                if "Remote end closed connection without response" in str(e): # this error comes up sometimes and a retry fixes it
                    if logger:
                        logger.warning(f"Error: {e.__class__.__name__}. Retrying. Attempt {attempt + 1}/3.")
                    wait_helper(attempt)
                    continue
                if "urlopen error" in str(e):
                    if logger:
                        logger.warning(f"Error: {e.__class__.__name__}: {e}. Attempt {attempt + 1}/3.")
                    wait_helper(attempt)
                    continue
                else:
                    if logger:
                        logger.error(f"Unexpected error: {e}")
                    raise
        if logger:
            logger.error("Max retries exceeded.")
        raise Exception("Max retries exceeded.")


def rate_limited_fetch(*args, **kwargs):
    """
    Wrapper for Entrez.efetch with rate-limiting.

    Args:
        *args, **kwargs: Arguments for the Entrez.efetch function.

    Returns:
        object: Result of the Entrez.efetch API call.
    """
    return rate_limited_call(Entrez.efetch, logger=logger, *args, **kwargs)


# moved here to keep rate limit fetch in one spot
def pubmed_api_fetch(pubmed_id):
    """
    Fetch metadata from PubMed for a given PubMed ID.

    Args:
        pubmed_id (str): The PubMed ID.

    Returns:
        object: Result of the Entrez.efetch API call for PubMed metadata.
    """
    pubmed_handle = rate_limited_call(Entrez.esummary, logger=logger, db="pubmed", id=pubmed_id)
    pub_med_record = Entrez.read(pubmed_handle)
    pubmed_handle.close()

    return pub_med_record


def fasta_fetch(record_id):
    """
    Fetch Sequence via fasta from nucleotide db.

    Args:
        record_id (str): The accession Number

    Returns:
        object: Result of the Entrez.efetch API call for Fasta Sequence.
    """
    if logger:
        logger.info(f"Downloading FASTA from db for ID: {record_id}")
    fasta_handle = rate_limited_call(Entrez.efetch, logger=logger, db="nucleotide", id=record_id,
                               rettype="fasta", retmode="text")
    fasta_sequence = fasta_handle.read()
    fasta_handle.close()

    return fasta_sequence


#######


def get_metadata(record):
    """
    Extract metadata from a GenBank record.

    Args:
        record (SeqRecord): The GenBank record to extract metadata from.

    Returns:
        dict: Metadata extracted from the record.
    """
    metadata = METADATA_TEMPLATE.copy()
    pubmed_record = False
    try:
        ref = record.annotations['references'][0]
        if ref.pubmed_id:
            pubmed_record = pubmed_api_fetch(ref.pubmed_id)
    except Exception as e:
        if logger:
            logger.error(f"Error extracting Reference info for ID {record.id}: {e}")


    metadata.update({
        "accession": record.id,
        **get_pubmed_info(record, pubmed_record, logger),
        **get_assembly_info(record, logger),
        **get_geo_info(record, logger)
    })

    return metadata


def process_entry(acc_id):
    """
    Main processing function of a single accession ID.
    Fetches, filters, gathers metadata and saves fasta for a profile.

    Args:
        acc_id (str): The accession number.

    Returns:
        tuple: (kept, accession_version, metadata or removal_reason)
    """
    try:
        logger.debug(f"Fetching information for ID: {acc_id}")
        if VERBOSE:
            print(f"Fetching information for ID: {acc_id}")
        start_time_fetch = time.time()
        handle = rate_limited_fetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        if VERBOSE:
            print(f"Fetching {acc_id} took {time.time() - start_time_fetch} seconds.")

        # start_time_proc = time.time()

        # we need to catch weird only 'N' Sequences
        # and replace them with fasta fetched sequences before (!) filtering
        record = ensure_correct_seq(record)

        # qc of seq
        rem_reason = apply_filters(record)
        if rem_reason:
            return False, acc_id, rem_reason

        metadata = get_metadata(record)

        write_seq_as_fasta(record)

        # print(f"Processing {acc_id} took {time.time() - start_time_proc} seconds.")

        return True, acc_id, metadata

    except Exception as e:
        logger.error(f"Error processing entry for ID {acc_id}: {e} ; skipping.")
        # necessary to keep `future.result()` iterable
        return []


def process_entries_sequential(id_list: list, batch_size: int):
    """
    Fetching wrapper function.
    Processes a list of accession numbers sequentially.

    Args:
        id_list (list): List of accession numbers to process.
        batch_size (int): Number of entries to process per batch.
    """
    filtered_entries = []
    removed = []
    metas = []

    if not id_list:
        return filtered_entries

    for index, accession_version in enumerate(id_list):
        result = process_entry(accession_version)

        if not result:
            continue

        kept, accession_version, result = result
        if kept:
            filtered_entries.append(accession_version)
            metas.append(result)
        else:
            removed.append({
                "accession": accession_version,
                "filter": result})

        # save metadata in batches
        if (index + 1) % batch_size == 0 or index == len(id_list) - 1:
            save_batch_info(index, filtered_entries, removed, metas, logger=logger)
            metas = []
            filtered_entries = []
            removed = []

    return 0


def process_entries_parallel(id_list: list, batch_size: int, NUM_WORKERS):
    """
    Fetching wrapper function for parallel fetching.
    Processes a list of accession numbers in parallel.

    Args:
        id_list (list): List of accession numbers to process.
        batch_size (int): Number of entries to process per batch.
        NUM_WORKERS (int): Number of worker threads for parallel processing.
    """
    start_time_batch = time.time()
    filtered_entries = []
    removed = []
    metas = []

    if not id_list:
        return

    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
        try:
            futures = {executor.submit(process_entry, acc_id): acc_id for acc_id in id_list}

            for index, future in enumerate(concurrent.futures.as_completed(futures)):
                try:
                    result = future.result(timeout=60)
                    if not result:
                        continue

                    kept, accession_version, data = result

                    if kept:
                        filtered_entries.append(accession_version)
                        metas.append(data)
                    else:
                        removed.append({
                            "accession": accession_version,
                            "filter": data})

                except concurrent.futures.TimeoutError:
                    logger.warning(
                        f"Timeout occurred for accession {futures[future]}. Task skipped."
                    )
                    continue

                # save batch info and reset lists
                if (index + 1) % batch_size == 0 or index == len(id_list) - 1:
                    try:
                        logger.info(f"Saving progress after {batch_size} Profiles, {index + 1} out of {len(id_list)}.")
                        start_time_writing = time.time()
                        save_batch_info(index, filtered_entries, removed, metas, logger=logger)
                        if VERBOSE:
                            print(f"Writing took {time.time() - start_time_writing} seconds.")
                    except Exception as e:
                        logger.critical(f"Fatal error during batch writing: {e}")
                        raise

                    metas = []
                    filtered_entries = []
                    removed = []
                    print(f"Processed {batch_size} Profiles in {time.time() - start_time_batch} seconds.")
                    start_time_batch = time.time()

        # ensure workers quit if there exists a problematic error, f.e. during batch writing
        except Exception as e:
            logger.critical(f"Shutting down due to fatal error: {e}")
            executor.shutdown(wait=False, cancel_futures=True)
            raise

    return 0


def apply_filters(record):
    """
    Apply quality control filters to a GenBank record.

    Args:
        record (SeqRecord): The GenBank record to filter.

    Returns:
        str or bool: Removal reason if the record is filtered out, or False if it passes all filters.
    """
    remove_reason = None
    for criteria in FILTERS:
        if criteria['fun'](record):
            remove_reason = criteria['description']
            break
    if remove_reason:
        logger.warning(f"Sequence ID {record.id} removed: {remove_reason}.")
        return remove_reason
    return False


def ensure_correct_seq(record):
    """
    Checks if `record.seq` is malformed.
    If so, refetches sequence via fasta from nucleotide db and replaces record's sequence.

    Args:
        record (SeqRecord): The GenBank record.

    Returns:
        record (SeqRecord): The GenBank record with possibly replaced sequence.
    """
    if record.seq.__class__.__name__ == "UnknownSeq" or record.seq.count("N") == len(record.seq):
        # recall up to n times if sequence is not valid
        try:
            fetched_record = find_clean_seq_fetcher(record)

            if fetched_record is None:
                raise ValueError("No clean sequence found.")

        except Exception as e:
            if logger:
                logger.error(f"Could not find clean FASTA for {record.id}.")
            raise ValueError(f"Error refetching Sequence for {record.id} : {e}.")

        # replace old seq with fasta fetched one
        record.seq = fetched_record.seq

    return record


# This was added on 25/01/16 because apparently on some pcs, running StringIO read operations
# while saving the batch info causes an incomplete read here.
# this is probably not handled well. For one because of the length of code necessary and secondly,
# because this essentially nests with the retry funion of `rate_limited_call`. So essentially,
# a 9x wait helper run through would be possible, which is a lot of time!
def find_clean_seq_fetcher(record):
    """
    Wrapper to refetch and read the sequence of a record as a fasta.
    Excepts and retries, if:
        - it fetches a seq of 'unknownSeq' class
        - reading the sequence via StringIO throws an error

    Args:
        record (SeqRecord): The GenBank record.

    Returns:
        record (Seq | String): The refetched Sequence if it finds a clean one within retries.
    """
    for attempt in range(3):
        try:
            fasta_str = fasta_fetch(record.id)
            if fasta_str.__class__.__name__ == "UnknownSeq":
                if logger:
                    logger.warning(f"Refetched Sequence of type `UnknownSeq`, retrying."
                                   f" Attempt {attempt + 1}/3.")
                wait_helper(attempt)
                continue
            record_io = StringIO(fasta_str)
            fetched_record = SeqIO.read(record_io, "fasta")
            record_io.close()

            return fetched_record

        except Exception as e:
            if logger:
                logger.warning(f"Error refetching Sequence: {e} ; retrying. Attempt {attempt + 1}/3.")
            wait_helper(attempt)
            continue

    raise ValueError(f"Max retries.")


def process_profiles(id_list:list, batch_size:int, parallel:bool, NUM_WORKERS):
    """
    Parent fetching wrapper function.
    Processes all profiles of a list of accession numbers, either sequentially or in parallel.

    Args:
        id_list (list): List of accession numbers to process.
        batch_size (int): Number of entries to process per batch.
        parallel (bool): Whether to process entries in parallel.
        NUM_WORKERS (int): Number of worker threads for parallel processing.
    """
    print(f"Fetching {len(id_list)} Profiles.")
    print(f"Fetching in Batch sizes of {batch_size}.")
    if parallel:
        process_entries_parallel(id_list, batch_size, NUM_WORKERS)
    else:
        process_entries_sequential(id_list, batch_size)


def soft_restart(id_list, MAX_NUM, SEARCH_TERM):
    """
    Filter out previously processed IDs and fetch new ones if necessary.

    Args:
        id_list (list): List of accession Numbers.
        MAX_NUM (int): Maximum number of profiles to fetch.
        SEARCH_TERM (str): NCBI search term for fetching IDs.

    Returns:
        list: Updated list of accession Numbers to process.
    """
    prev_processed = load_processed_ids(logger)
    if len(prev_processed):
        print(f"Restarting softly with previously processed Profiles.")

        id_list = filter_unprocessed_ids(id_list, prev_processed, logger)

        if len(id_list) < MAX_NUM:
            # fill up to MAX_NUM number to process
            # This was edited on 24/01/16 to not return 0 new profile sto fetch, when you make consecutive
            # max nums smaller. However, this is noted here, as i had to 'fix' this before to the previous version.
            # So perhaps this version now will again bring back the previous problem.
            max_num = len(prev_processed) + (MAX_NUM - len(id_list))
            logger.info(f"Some or all fetches already in the processed file. Raising up to {MAX_NUM} new ones and refetching.")
            if max_num > LIMIT_NUM:
                # edited on 25/01/16 to no longer throw an error, but to use LIMIT_NUM instead
                limit_msg = ("Achieved maximum size for automatic size increase with soft restarts.\n"
                             "This may indicate unwanted use of --soft-restart ?")
                logger.warning(limit_msg)
                warnings.warn(limit_msg)
                logger.warning(f"Changing to fetch up to {LIMIT_NUM} profiles only.")

                max_num = LIMIT_NUM

            id_list = filter_unprocessed_ids(fetch_profile_accs(SEARCH_TERM, max_num=max_num, logger=logger),
                                             prev_processed, logger)

    return id_list


def parse_args():
    parser = argparse.ArgumentParser(description="Fetch mitochondrial data.")
    parser.add_argument("--max-num", type=int, default=MAX_NUM, help="Maximum number of profiles to fetch.")
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE, help="Batch size for processing.")
    parser.add_argument("--fetch-parallel", action="store_true", default=FETCH_PARALLEL, help="Enable parallel fetching.")
    parser.add_argument("--num-workers", type=int, default=NUM_WORKERS, help="Number of workers for parallel fetching.")
    parser.add_argument("--soft-restart", action="store_true", default=SOFT_RESTART, help="Restart softly with previously processed profiles.")
    parser.add_argument("--search-term", type=str, default=SEARCH_TERM, help="NCBI search term.")
    parser.add_argument("--clean-dir", action="store_true", default=CLEAN_DIR, help="Set this flag to clean the output directory,"
                                                                                    " removing profiles not within the current calls fetch.")
    parser.add_argument("--y", action="store_true", default=False,
                        help="Answer all user inputs with yes.")

    return parser.parse_args()


def check_valid_inputs(arg_dict):
    # check we dont exceed limit number of fetches
    if arg_dict["max_num"] > LIMIT_NUM or arg_dict["max_num"] < 1:
        raise ValueError(f"Supplied `max-num` is not within the bounds of [1, {LIMIT_NUM}].\n"
                         f"Choose a number of Profiles to fetch with in those bounds.")

    # TODO possibly add more checking other values here


def main():
    # CLI arguments
    args = parse_args()

    MAX_NUM = args.max_num
    BATCH_SIZE = args.batch_size
    FETCH_PARALLEL = args.fetch_parallel
    NUM_WORKERS = args.num_workers
    SOFT_RESTART = args.soft_restart
    SEARCH_TERM = args.search_term
    CLEAN_DIR = args.clean_dir
    FORCE = args.y

    check_valid_inputs({"max_num": MAX_NUM, "batch_size": BATCH_SIZE, "num_workers":NUM_WORKERS,
                        "search_term": SEARCH_TERM, "clean_dir": CLEAN_DIR})

    start_time = time.time()
    try:
        id_list = fetch_profile_accs(SEARCH_TERM, max_num=MAX_NUM, logger=logger)
        num_all = len(fetch_profile_accs(SEARCH_TERM, max_num=LIMIT_NUM, logger=logger))

        print(f"\nProvided Search Term returns {num_all} total profiles.\n"
              f"{len(id_list)} will be attempted to be fetched.\n")

        if not FORCE:
            clean_msg = ("\033[93m\nExecuting the fetcher with the `clean-dict` flag set will remove any profiles from all data files "
                     "if they are not within the current search Query! This will potentially delete a lot of profile data!\n"
                     f"Only those within the {num_all} total returned matches will be kept after this action.\n"
                     "If you are unsure if this is the correct action to take, "
                     "make sure to make a copy of the data directory.\033[0m\n") if CLEAN_DIR else ""

            if input(f"If this seems incorrect, adjust the `max_num` or other parameters.\n"
                     f"{clean_msg}\n"
                     f"Do you wish to continue? (y/n)") != "y":
                print("\nAborting.")
                exit()

        print("\nExecuting Fetch!\n")

    except Exception as e:
        exit_msg = (f"\033[91mFatal Error when initially fetching IDs from Databank: {e}. Aborting.\n"
                     f"This is likely caused by the genebank api, please retry again later or "
                     f"contact the developers if this error persists.\033[0m")
        logger.error(exit_msg)
        print(exit_msg)

        raise SystemExit(1)

    # soft restart by loading processed IDs and filtering the list
    if SOFT_RESTART:
        id_list = soft_restart(id_list, MAX_NUM, SEARCH_TERM)

    # filter out profiles that do not need to be fetched, as their current version is up to date
    local_versions = load_local_versions(logger)
    removed_local = load_removed_versions(logger)
    id_list = filter_changed_profiles(id_list, local_versions, removed_local, logger=logger)

    # readd profiles, whose metadata has been changed since the last run
    last_run_date = get_last_run_date(logger=logger)
    if last_run_date:
        id_list = readd_recently_modified_profiles(SEARCH_TERM, id_list, last_run_date)

    # execute fetching and writing process
    process_profiles(id_list, BATCH_SIZE, FETCH_PARALLEL, NUM_WORKERS)
    print("Finished fetching all Profiles.")
    print(f"Processed {len(id_list)} Profiles in {time.time() - start_time} seconds.")

    # remove duplicates from metadata file
    post_process_metadata(logger)

    # cleanup
    write_last_run_date(run_date=RUN_TIME, logger=logger)

    try:
        full_id_list = fetch_profile_accs(SEARCH_TERM, max_num=LIMIT_NUM, logger=logger)
    except Exception as e:
        logger.error(f"Problem fetching full accession list for post-process check: {e}. Continuing without.")
        full_id_list = []

    if CLEAN_DIR:
        print("Cleaning data directory based on current query.")
        logger.info(f"Cleaning data directory by removing any profile data not present in the current query with the "
                   f"search-term: {SEARCH_TERM} and a max-num of {LIMIT_NUM}.")
        if full_id_list:
            clean_profiles_from_data(full_id_list, logger=logger)
        else:
            logger.error("Fetch query missing or empty, Cannot clean directory without successful ncbi query. Skipping.")

    post_process_check(logger=logger, ids_list=full_id_list)
    cleanup_old_files(PROCESSED_IDS_DIR, 3, logger)
    cleanup_old_files(LOG_DIR, 5, logger)
    cleanup_old_files(DEBUG_DIR, 8, logger)

    return 0


if __name__ == "__main__":
    main()
