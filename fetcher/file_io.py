# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
file_io.py

File Input/Output Utility Functions for mitoFetch

This script contains functions for managing file operations related to processed IDs,
metadata, removed IDs, and sequence files. It also includes utilities for cleaning up
old files and post-processing metadata.

Key Features:
- Save and load processed IDs, removed ids and metadata.
- Write sequences as cleaned FASTA files.
- Post-process metadata entries, including duplicate removal and version control.
- Clean up outdated files.

Dependencies: Biopython, Pandas.

Usage:
- Import functions as needed into other scripts for processing and file management.
- Ensure appropriate paths and global constants are configured in `global_defaults`.

Author: Noah Hurmer (aka. minimops) as part of the mitoTree Project.
"""


import os
import pandas as pd
from datetime import date, datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .global_defaults import (CURRENT_PROCESSED_IDS_FILE,
                             METADATA_FILE,
                             LAST_RUN_PATH,
                             SEQS_DIR,
                             REMOVED_IDS_FILE,
                             IDS_FILE,
                             PROCESSED_IDS_DIR)


def save_processed_ids(processed_ids, logger=None):
    """
    Save processed sequence IDs to a new processed IDs file,
    carrying over the processed IDs from the old file.

    Args:
        processed_ids (list): List of sequence IDs to save.
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.
    """
    file_path = CURRENT_PROCESSED_IDS_FILE
    old_processed_ids = load_processed_ids()

    # combine old and new IDs and write them to file
    # remove duplicates
    all_processed_ids = list(dict.fromkeys(processed_ids + old_processed_ids))

    with open(file_path, 'w') as f:
        for seq_id in all_processed_ids:
            f.write(f"{seq_id}\n")
    if logger:
        logger.info(f"Saved {len(all_processed_ids)} processed seq_ids to {file_path}. (Includes {len(old_processed_ids)} previously processed IDs.)")


def find_latest_processed_ids_file(logger=None):
    """
    Find the most recent processed IDs file in the processed IDs directory.

    Args:
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.

    Returns:
        str or None: Path to the most recent processed IDs file, or None if none exist.
    """
    try:
        files = [
            os.path.join(PROCESSED_IDS_DIR, file)
            for file in os.listdir(PROCESSED_IDS_DIR)
            if file.endswith(".txt")
        ]
        if not files:
            if logger:
                logger.warning("No processed IDs file found for soft restart.")
            return None
        latest_file = max(files, key=os.path.getctime)
        if logger:
            logger.info(f"Using the most recent processed IDs file: {latest_file}")
        return latest_file
    except Exception as e:
        if logger:
            logger.error(f"Error finding latest processed IDs file: {e}")
        return None


def load_processed_ids(logger=None):
    """
    Load processed sequence IDs from the most recent processed IDs file.

    Args:
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.

    Returns:
        list: List of processed sequence IDs.
    """
    file_path = find_latest_processed_ids_file(logger)
    if file_path:
        with open(file_path, 'r') as f:
            processed_ids = [line.strip() for line in f]
        if logger:
            logger.info(f"Loaded {len(processed_ids)} processed seq_ids from {file_path}.")
        return processed_ids
    return []


def filter_unprocessed_ids(id_list, processed_ids, logger=None):
    """
    Filter out already processed sequence IDs from the supplied list of accession numbers.

    Args:
        id_list (list): List of sequence IDs to process.
        processed_ids (list): List of already processed sequence IDs.
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.


    Returns:
        list: List of unprocessed sequence IDs.
    """
    unprocessed_ids = [seq_id for seq_id in id_list if seq_id not in processed_ids]
    if logger:
        logger.info(f"Filtered out {len(id_list) - len(unprocessed_ids)} already-processed IDs."
                    f" Remaining: {len(unprocessed_ids)}.")
    return unprocessed_ids


def save_batch_info(index, filtered_entries, removed, metas, logger=None):
    """
    Save batch information of fetched profiles.

    Args:
        index (int): Index of the current batch.
        filtered_entries (list): List of sequence IDs that passed filters.
        removed (list): List of sequence IDs that failed filters with reasons.
        metas (list): Metadata for the processed sequences.
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.
    """
    print("Saving batch info.")

    update_local_versions(filtered_entries, logger)
    save_metadata(metas, logger)
    save_removed_versions(removed, logger)

    all_processed_ids = [rem.get("accession") for rem in removed] + filtered_entries
    save_processed_ids(all_processed_ids)

    print(f"Progress saved after processing {index + 1} entries.")



def split_accession(accession):
    """
    Helper to split an accession string into ID and version.

    Args:
        accession (str): Accession string in the format 'ABC123456.1'.

    Returns:
        tuple: (accession_id, version) where version is an integer.
    """
    try:
        accession_id, version = accession.split('.')
        return accession_id, int(version)
    except ValueError:
        raise ValueError(f"Invalid accession format: {accession}")


def save_metadata(new_metas, logger=None):
    """
    Save metadata to the metadata file.

    Args:
        new_metas (list or DataFrame): New metadata entries to save.
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.
    """
    meta_df = pd.DataFrame(new_metas)
    meta_df.to_csv(METADATA_FILE, mode='a', header=not os.path.exists(METADATA_FILE), index=False)
    if logger:
        logger.info(f"{len(new_metas)} metadata entries written to {METADATA_FILE}.")


def load_local_versions(logger=None):
    """
    Load local versions of sequence IDs from the local versions file.

    Args:
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.

    Returns:
        dict: Dictionary mapping accession IDs to their versions.
    """
    if not os.path.exists(IDS_FILE):
        return {}
    local_versions = {}
    with open(IDS_FILE, 'r') as f:
        for line in f:
            accession, version = line.strip().split('.')
            local_versions[accession] = int(version)
    if logger:
        logger.info(f"Loaded {len(local_versions)} local versions.")
    return local_versions


def duplicate_removal(entries, logger=None):
    """
    Remove duplicate entries based on accession and version.

    Args:
        entries (DataFrame): DataFrame of entries with accession and version columns.
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.

    Returns:
        DataFrame: Deduplicated DataFrame.
    """
    # remove identical accession+version, keeping last
    prev_len = len(entries)
    entries.drop_duplicates(subset="accession", keep="last", inplace=True)
    if logger:
        logger.warning(f"{len(entries) - prev_len} duplicate entries dropped from entries.")
    entries[["id", "version"]] = entries['accession'].apply(lambda x: pd.Series(split_accession(x)))
    # remove older version of same accession
    prev_len = len(entries)
    entries = entries.loc[entries.groupby('id')['version'].idxmax()].drop(columns=['id', 'version'])
    if logger:
        logger.warning(f"{len(entries) - prev_len} older versions dropped, where newer entries exist.")

    # TODO catch if we download an older version?

    return entries


def update_local_versions(entries, logger=None):
    """
    Update the local versions file with new entries, removing duplicates and older versions.

    Args:
        entries (list): List of new entries to add to the local versions file.
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.
    """
    with open(IDS_FILE, 'a') as f:
        for entry in entries:
            f.write(f"{entry}\n")

    # remove duplicates here
    with open(IDS_FILE, 'r') as f:
        ex_ids = f.read().strip().split("\n")

    if not len(ex_ids):
        return

    ids_df = duplicate_removal(pd.DataFrame(ex_ids, columns=["accession"]))
    ids_df.to_csv(IDS_FILE, index=False, header=False)
    if logger:
        logger.info(f"Updated local versions in {IDS_FILE}.")


def load_removed_versions(logger=None):
    """
    Load removed sequences from the removed ids file.

    Args:
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.

    Returns:
        dict: Dictionary mapping accession IDs to their versions.
    """
    if not os.path.exists(REMOVED_IDS_FILE):
        return []
    local_removed = {}
    with open(REMOVED_IDS_FILE, 'r') as f:
        next(f)
        for line in f:
            accession_num, _ = line.strip().split(',')
            accession, version = accession_num.strip().split('.')
            local_removed[accession] = int(version)
    if logger:
        logger.info(f"Loaded {len(local_removed)} local removed versions.")
    return local_removed


def save_removed_versions(removed_entries, logger=None):
    """
    Save removed sequences to the removed ids file.

    Args:
        removed_entries (list): List of removed entries with accession and filter reason.
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.
    """
    if not len(removed_entries):
        return

    removed_entries = pd.DataFrame(removed_entries, columns=['accession', 'filter'])
    if os.path.exists(REMOVED_IDS_FILE):
        prev_removed = pd.read_csv(REMOVED_IDS_FILE)
    else:
        prev_removed = pd.DataFrame(columns=removed_entries.columns)

    removed = pd.concat([prev_removed, removed_entries], ignore_index=True)
    removed = duplicate_removal(removed)

    # save back to the file
    removed.to_csv(REMOVED_IDS_FILE, index=False)
    if logger:
        logger.info(f"Saved {len(removed_entries)} entries to {REMOVED_IDS_FILE}.")


def cleanup_old_files(directory, keep_last=3, logger=None):
    """
    Cleans up old files in the specified directory, keeping only the last `keep_last` files.

    Args:
        directory (str): Path to the directory containing the files.
        keep_last (int): Number of most recent files to keep. Defaults to 3.
        logger (logging.Logger, optional): Logger for logging messages. Defaults to None.
    """
    try:
        # List all files in the directory with the specified prefix
        print(f"Cleaning up files in {directory}.")

        files = [
            os.path.join(directory, file)
            for file in os.listdir(directory)
        ]
        if len(files) <= keep_last:
            if logger:
                logger.info(f"No cleanup needed for {directory}. Found {len(files)} files.")
            return

        # Sort files by their modification time, keeping the last `keep_last`
        files.sort(key=os.path.getctime, reverse=True)
        files_to_remove = files[keep_last:]

        # Remove old files
        for file in files_to_remove:
            os.remove(file)
            if logger:
                logger.info(f"Removed old file: {file}")

        if logger:
            logger.info(f"Cleanup complete for {directory}. Kept {keep_last} most recent files.")

    except Exception as e:
        if logger:
            logger.error(f"Error during cleanup of {directory}: {e}")


def write_seq_as_fasta(record, logger=None):
    """
    Cleans the sequence by replacing all 'D' with '-' and writes it to a FASTA file.

    Args:
        record (SeqRecord): A Biopython SeqRecord containing the sequence data.
        logger (logging.Logger, optional): Logger for logging errors. Defaults to None.

    Returns:
        bool: True if the operation was successful, False otherwise.
    """
    try:
        # this is a bit annoying,but in order to actually manipulate the sequence,
        # the record needs to be remade.
        # alternative would be to read the file lines as strings, not sure if thats better
        cleaned_sequence = clean_sequence(record.seq)

        cleaned_record = SeqRecord(
            Seq(cleaned_sequence),
            id=record.id,
            name=record.name,
            description=record.description
        )

        SeqIO.write(cleaned_record, f"{SEQS_DIR}/{record.id.split('.')[0]}.fasta", "fasta")

        if logger:
            logger.debug(f"FASTA for {record.id} written.")
        return True

    except Exception as e:
        if logger:
            logger.error(f"Error writing FASTA for ID {record.id}: {e}")
        else:
            print(f"Error writing FASTA for ID {record.id}: {e}")
        return False


def clean_sequence(sequence):
    """
    Clean a sequence string by defined rules.
    Currently just replacing all occurrences of 'D' with '-'.

    Args:
        sequence (str): The sequence string to clean.

    Returns:
        str: Cleaned sequence string.
    """
    return sequence.replace('D', '-')


def get_last_run_date(file_path=LAST_RUN_PATH, logger=None):
    """
    Reads the last run date from the specified file.

    Args:
        file_path (str): Path to the file containing the last run date.
        logger (logging.Logger, optional): Logger for writing messages. Defaults to None.

    Returns:
        datetime.date: The last run date, or None if an error occurred or the file does not exist.
    """
    try:
        # Check if the file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File '{file_path}' not found.")

        # Read and process the date
        with open(file_path, 'r') as f:
            date_str = f.read().strip()

        # Convert to datetime.date
        last_run_date = datetime.strptime(date_str, "%Y-%m-%d").date()

        if logger:
            logger.info(f"Last run date fetched: {last_run_date}.")

        return last_run_date

    except Exception as e:
        message = f"Failed to read date from {file_path}: {e}"
        if logger:
            logger.error(message)
        else:
            print(message)

    return None


def write_last_run_date(file_path=LAST_RUN_PATH, run_date=None, logger=None):
    """
    Writes the given date or the current date to the specified file.

    Args:
        file_path (str): Path to the file where the date will be saved.
        run_date (datetime.date, optional): The date to save. If None, the current date is used. Defaults to None.
        logger (logging.Logger, optional): Logger for writing messages. Defaults to None.

    Returns:
        bool: True if the date was written successfully, False otherwise.
    """
    try:
        if run_date is None:
            run_date = date.today()

        run_date = run_date.strftime("%Y-%m-%d")
        with open(file_path, 'w') as f:
            f.write(run_date)

        if logger:
            logger.info(f"Date {run_date} written to {file_path}.")

        return True

    except Exception as e:
        message = f"Failed to write date to {file_path}: {e}"
        if logger:
            logger.error(message)
        else:
            print(message)

        return False


def post_process_metadata(logger=None):
    """
    Perform post-processing on the metadata file:
    - Remove duplicate accessions.
    - Keep the highest version per ID.
    - Log warnings for irregularities.

    Args:
        logger (logging.Logger, optional): Logger for logging progress. Defaults to None.
    """
    if not os.path.exists(METADATA_FILE):
        if logger:
            logger.error(f"Metadata file '{METADATA_FILE}' not found.")
        return

    meta_df = pd.read_csv(METADATA_FILE)
    if logger:
        logger.info(f"Loaded {len(meta_df) - 1} rows from metadata file.")

    print("Performing post-processing of metadata file.")
    if logger:
        logger.info(f"Performing Post-processing of metadata file.")

    # add helper columns
    meta_df['index'] = meta_df.index # order of fetched
    meta_df[['id', 'version']] = meta_df['accession'].apply(lambda x: pd.Series(split_accession(x)))

    meta_df = meta_df.sort_values(by=['id', 'version', 'index'], ascending=[True, False, True])

    # keep the highest version per group
    rows_to_keep = []
    for _, group in meta_df.groupby('id'):
        # highest version row
        best_row = group.iloc[0]
        # check if highest version is the most recent
        if best_row['index'] != group['index'].max():
            if logger:
                logger.warning(f"Highest version for ID {best_row['id']} ({best_row['accession']}) "
                               f"is not the most recently added row."
                               f"Keeping highest version anyway.")

        # check highest version row has the most filled fields
        for _, row in group.iloc[1:].iterrows():
            if row.notna().sum() > best_row.notna().sum():
                if logger:
                    logger.warning(f"Row {row['index']} with ID {row['id']} and version {row['version']} "
                                   f"has more filled fields than the highest version"
                                   f" ({best_row['accession']} @ {best_row['index']}).")

        # split up any tied highest versions (same ID, same version) by most recently added
        max_version_rows = group[group['version'] == best_row['version']]
        best_rows = max_version_rows.loc[max_version_rows['index'] == max_version_rows['index'].max()]
        rows_to_keep.append(best_rows)

    final_df = pd.concat(rows_to_keep).drop(columns=['id', 'version', 'index']).reset_index(drop=True)
    if logger:
        logger.info(f"Removed {len(meta_df) - len(final_df)} duplicate rows.")

    # save cleaned metadata
    final_df.to_csv(METADATA_FILE, index=False)
    if logger:
        logger.info(f"Post-processing complete. Saved {len(final_df) - 1} rows to '{METADATA_FILE}'.")
