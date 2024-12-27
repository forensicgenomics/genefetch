# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
post_process_check.py

This script performs post-processing checks on files and directories used in the data pipeline.
It verifies consistency, detects duplicates, and ensures the expected relationships between
ids_list, metadata, and sequence files are maintained.

Dependencies: pandas.

Author: Noah Hurmer as part of the mitoTree Project.
"""

import os
import pandas as pd

from .global_defaults import IDS_FILE, METADATA_FILE, REMOVED_IDS_FILE, SEQS_DIR


def check_duplicates_in_file(file_path, columns=None, logger=None):
    """
    Check for duplicate rows in a specified file based on given columns.

    Args:
        file_path (str): Path to the file to check.
        columns (list or str): Columns to check for duplicates. None checks entire rows.
        logger (logging.Logger, optional): Logger instance for logging errors and warnings.
    """
    try:
        if not os.path.exists(file_path):
            if logger:
                logger.warning(f"File '{file_path}' not found. Skipping duplicate check.")
            return

        df = pd.read_csv(file_path)

        if columns:
            duplicate_rows = df[df.duplicated(subset=columns, keep=False)]
        else:
            duplicate_rows = df[df.duplicated(keep=False)]

        if not duplicate_rows.empty:
            if logger:
                logger.warning(f"Duplicate rows found in '{file_path}' (columns: {columns}):\n{duplicate_rows}")

    except Exception as e:
        if logger:
            logger.error(f"Error checking duplicates in file '{file_path}': {e}")


def check_removed_and_ids_list(ids_list, logger=None):
    """
    Verify that the provided ids_list (eSearch of all matches) matches the combination of accessions in ids_list
    and removed_ids in those two files.

    Args:
        ids_list (list): List of accessions to verify.
        logger (logging.Logger, optional): Logger instance for logging errors and warnings.
    """
    try:
        # load ids_list
        ids_df = pd.read_csv(IDS_FILE, header=None, names=["accession"])
        ids_accessions = set(ids_df['accession'])

        # load removed_ids
        removed_df = pd.read_csv(REMOVED_IDS_FILE)
        removed_accessions = set(removed_df['accession'])

        # combine
        combined_accessions = set(ids_accessions) | removed_accessions

        ids_set = set(ids_list)
        missing_in_combined = ids_set - combined_accessions
        extra_in_combined = combined_accessions - ids_set

        # how many to show
        n = 5
        if missing_in_combined:
            logger.warning(
                f"{len(missing_in_combined)} accessions missing in combined list. Showing first {n}: {list(missing_in_combined)[:n]}")
        if extra_in_combined:
            logger.warning(
                f"{len(extra_in_combined)} extra accessions in combined list. Showing first {n}: {list(extra_in_combined)[:n]}")

    except Exception as e:
        logger.error(f"Error checking ids_list and removed_ids: {e}")


def check_metadata_and_ids_list(logger=None):
    """
    Check that all accessions in ids_list file are present in the metadata file
    and vice versa.

    Args:
        logger (logging.Logger, optional): Logger instance for logging errors and warnings.
    """
    try:
        # Load ids_list
        ids_df = pd.read_csv(IDS_FILE, header=None, names=["accession"])
        ids_accessions = set(ids_df['accession'])

        # Load metadata
        metadata_df = pd.read_csv(METADATA_FILE)
        metadata_accessions = set(metadata_df['accession'])

        # Compare
        if ids_accessions == metadata_accessions:
            msg = "ids_list accessions match exactly with metadata accessions."
            logger.info(msg) if logger else print(msg)
        else:
            missing = ids_accessions - metadata_accessions
            extra = metadata_accessions - ids_accessions
            msg = f"Mismatch between ids_list and metadata accessions.\nMissing: {missing}\nExtra: {extra}"
            logger.warning(msg) if logger else print(msg)

    except Exception as e:
        emsg = f"Error checking ids_list and metadata: {e}"
        logger.error(emsg) if logger else print(emsg)


def check_seqs_directory(logger=None):
    """
    Verify that the sequence directory contains exactly the expected files
    matching the accessions in ids_list.

    Args:
        logger (logging.Logger, optional): Logger instance for logging errors and warnings.
    """
    try:
        # load ids_list
        ids_df = pd.read_csv(IDS_FILE, header=None, names=["accession"])
        ids_accessions = {acc.split(".")[0] for acc in ids_df['accession']}

        # list of files in seqs directory
        seq_files = {f.replace('.fasta', '') for f in os.listdir(SEQS_DIR) if f.endswith('.fasta')}

        if seq_files == ids_accessions:
            msg = "Sequence directory has exactly the expected files."
            logger.info(msg) if logger else print(msg)

        else:
            missing = ids_accessions - seq_files
            extra = seq_files - ids_accessions
            msg = f"Mismatch in sequence directory.\nMissing files: {missing}\nExtra files: {extra}"
            logger.warning(msg) if logger else print(msg)

    except Exception as e:
        emsg = f"Error checking sequence directory: {e}"
        logger.error(emsg) if logger else print(emsg)


def main(logger=None, ids_list=None):
    print("Performing postprocess check.")
    if logger:
        logger.info("Performing post-process check.")

    check_duplicates_in_file(IDS_FILE, logger=logger)
    check_duplicates_in_file(REMOVED_IDS_FILE, columns="accession", logger=logger)
    if ids_list:
        check_removed_and_ids_list(ids_list, logger=logger)
    check_metadata_and_ids_list(logger=logger)
    check_seqs_directory(logger=logger)

    print("postprocess check complete.")
    if logger:
        logger.info("post-process check complete.")


if __name__ == "__main__":
    import logging

    # setup logger
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger("PostProcessLogger")

    main(logger)
