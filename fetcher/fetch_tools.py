# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
fetch_tools.py

Utility Tools related to Entrez interaction for mitoFetch Pipeline

This script provides helper functions for interacting with the NCBI Entrez API and managing
mitochondrial profile fetching. It includes utilities for setting API credentials, filtering
profiles based on version updates or removal, and fetching metadata for modified or newly added
profiles.

Dependencies: Biopython.

Usage:
- Import the functions into other scripts to streamline profile fetching and filtering.
- Ensure NCBI API credentials (API key or email) are configured in the appropriate secret files.

Author: Noah Hurmer as part of the mitoTree Project.
"""

# helpers for pipeline
import os.path
from datetime import date
import time
from Bio import Entrez
from .global_defaults import (API_DEST,
                              EMAIL_DEST,
                              SEQS_DIR,
                              CHECK_FOR_FASTA)


def set_entrez_globals(api:bool=True, email:bool=True):
    """
    Configure Entrez global settings for API access.

    Args:
        api (bool): Whether to use the NCBI API key for authentication. Defaults to True.
        email (bool): Whether to use the email address for Entrez access. Defaults to True.

    Raises:
        Exception: If no API key or email is provided or if the specified files are missing.
    """
    if (not api and not email) or (not os.path.isfile(API_DEST) and not os.path.isfile(EMAIL_DEST)):
        raise Exception("Supply either API key or email address to access NCBI api via Entrez.\n")

    if os.path.isfile(API_DEST) and api:
        Entrez.api_key = read_ncbi_api_key()
    if os.path.isfile(EMAIL_DEST) and email:
        Entrez.email = read_ncbi_email()


def parse_secret(dest:str, sec_type:str):
    """
    Read and parse secret information from a file.

    Args:
        dest (str): Path to the secret file.
        sec_type (str): Type of the secret, e.g., "email" or "api_key".

    Returns:
        str: The secret value.
    """
    try:
        f = open(dest, "r")
        auth_key = f.read().strip()
    except FileNotFoundError:
        raise Exception(f"file '{dest}' not found. Save your {sec_type} key to (that) location and retry.\n")
    else:
        f.close()

    return auth_key


def read_ncbi_email(dest:str=EMAIL_DEST) -> str:
    """
    Read and return the NCBI email address from the secret file.

    Args:
        dest (str): Path to the email secret file. Defaults to EMAIL_DEST global variable.

    Returns:
        str: The email address.
    """
    entrez_email = parse_secret(dest, "email")

    assert isinstance(entrez_email, str), "Key must be string"

    return entrez_email


def read_ncbi_api_key(dest:str=API_DEST) -> str:
    """
    Read and return the NCBI API key from the secret file.

    Args:
        dest (str): Path to the API key secret file. Defaults to API_DEST.

    Returns:
        str: The API key.
    """
    auth_key = parse_secret(dest, "api_key")

    # check for length
    assert isinstance(auth_key, str), "Key must be string"
    assert len(auth_key) == 36, "Auth key of incorrect length"

    return auth_key


def set_entrez_rate():
    """
    Determine the number of API calls allowed per second based on authentication.

    Returns:
        int: The rate of API calls per second. Defaults to 1 without credentials,
             increases to 3 with email, and 10 with an API key.
    """
    calls_per_second = 1
    if Entrez.api_key:
        calls_per_second = 10
    elif Entrez.email:
        calls_per_second = 3

    return calls_per_second


def filter_changed_profiles(accession_list, local_versions, local_removed_versions=None, logger=None):
    """
    Filter out profiles that do not need to be fetched based on their version that is already present locally.

    Args:
        accession_list (list): List of accession numbers to process.
        local_versions (dict): Mapping of accession numbers to their latest local versions.
        local_removed_versions (dict, optional): Mapping of removed accession IDs to their versions. Defaults to None.
        logger (logging.Logger, optional): Logger for logging messages. Defaults to None.

    Returns:
        list: Filtered list of accession IDs to process.
    """
    filtered_list = []
    filtered_out = []

    for accession in accession_list:
        accession_num, version = accession.split('.')
        version = int(version)
        local_version = local_versions.get(accession_num, 0)

        # skip if the version is not newer or fasta is not present if desired
        fasta_check = os.path.exists(f"{SEQS_DIR}/{accession.split('.')[0]}.fasta") if CHECK_FOR_FASTA else True
        if version <= local_version and fasta_check:
            if logger:
                logger.debug(f"Accession ID {accession} already downloaded and is up to date.")
            filtered_out.append(accession)
            continue
        # also remove if that version is in the removed list
        removed_version = local_removed_versions.get(accession_num, 0) if local_removed_versions else 0
        if version <= removed_version:
            if logger:
                logger.debug(f"Accession ID {accession} of current version previously removed.")
            filtered_out.append(accession)
            continue
        if local_version != 0 and version > local_version:
            if logger:
                logger.debug(f"Accession ID {accession} is an update of a previous version.")

        filtered_list.append(accession)

    if logger:
        logger.info(f"{len(accession_list) - len(filtered_list)} IDs dont need to be processed and are filtered out.")

    return filtered_list, filtered_out


def _accn_or_query(ids):
    ids = [i.strip() for i in ids if i and i.strip()]
    if not ids:
        return ""
    return "(" + " OR ".join(f"{i}[ACCN]" for i in ids) + ")"


def recently_modified_profiles(currently_filtered_out:list, last_change_date:date, max_num, logger=None,
                               chunk_size = 100):
    """
    Returns ids of currently_filtered_out list that were modified since the last run.

    Args:
        currently_filtered_out (list): List of current filtered out accession numbers.
        last_change_date (datetime.date): Date of the last metadata modification.
        logger (logging.Logger, optional): Logger for logging messages. Defaults to None.

    Returns:
        list: modified profiles.
    """
    if not currently_filtered_out:
        return []

    num_days_passed = abs((date.today() - last_change_date).days) + 1
    if logger:
        logger.info(f"Fetching modified Profiles from the last {num_days_passed} days to refetch metadata.")

    modified = set()
    # fetch by id in chunks using helper fun
    for i in range(0, len(currently_filtered_out), chunk_size):
        batch = currently_filtered_out[i:i+chunk_size]
        term = _accn_or_query(batch)
        hits = fetch_profile_accs(term, n_days=num_days_passed, logger=logger) or []
        modified.update(hits)

    # intersect with modified set
    updates = [acc for acc in currently_filtered_out if acc in modified]

    return updates


def readd_recently_modified_profiles(ids_list:list, currently_filtered_out:list, last_change_date:date, max_num, logger=None):
    """
        Re-add profiles that are filtered out from fetching as they are downloaded and have the current accession version,
        yet have been modified since the last run to the processing list. Used to fetch altered metadata.

        Args:
            ids_list (list): List of current accession numbers.
            currently_filtered_out (list): List of currently filtered out accession numbers from existing data)
            last_change_date (datetime.date): Date of the last metadata modification.
            logger (logging.Logger, optional): Logger for logging messages. Defaults to None.

        Returns:
            list: Updated list of accession numbers including modified profiles.
    """
    updates = recently_modified_profiles(currently_filtered_out, last_change_date, max_num, logger=logger)

    if updates:
        msg = f"Adding {len(updates)} modified Profiles to be updated in the process."
        print(msg + "\n")
        if logger:
            logger.info(msg)

        ids_list = list(set(updates).union(set(ids_list)))

    return ids_list


# TODO hard coded max num
def fetch_profile_accs(search_term, max_num=100000, n_days=None, logger=None):
    """
    Fetch accession numbers using the NCBI Entrez API.

    Args:
        search_term (str): Search term to use for fetching accessions.
        max_num (int, optional): Maximum number of IDs to fetch. Defaults to 10k since Entrez defaults this to 20 otherwise!
        n_days (int, optional): Fetch profiles modified in the last `n_days`. Defaults to None.
        logger (logging.Logger, optional): Logger for logging messages. Defaults to None.

    Returns:
        list: List of accession numbers fetched from Entrez.
    """
    attempt = 0
    retries = 3
    while attempt < retries:
        try:
            if logger:
                logger.info(f"Fetching IDs with search term: {search_term}")

            if n_days:
                handle = Entrez.esearch(db="nucleotide", idtype="acc", term=search_term, retmax=max_num, datetype="mdat",
                                        reldate=n_days)
            else:
                handle = Entrez.esearch(db="nucleotide", idtype="acc", term=search_term, retmax=max_num)

            record = Entrez.read(handle)
            handle.close()
            if logger:
                logger.info(f"Found {record.get('Count')} matches.")
            return record['IdList']

        except Exception as e:
            attempt += 1
            if attempt == retries:
                if logger:
                    logger.error(f"Error fetching IDs: {e}")
                raise e
            else:
                if logger:
                    logger.warning(f"Fetching IDs Attempt {attempt} failed: {e} ; retrying.")
                time.sleep(1)