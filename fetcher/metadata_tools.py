# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
metadata_tools.py

This script provides utility functions to extract and process metadata information
from records. Functions extract metadata such as publication
details, geographic origin, assembly information, and author details from provided
record objects.
Additional helper functions are provided for parsing dates and author names.

Author: Noah Hurmer as part of the mitoTree Project.
"""

import re
from datetime import datetime


def get_pubmed_info(record, pubmed_record, logger=None):
    """
    Extract PubMed-related metadata from a record and PubMed record object.

    Args:
        record (Bio.SeqRecord): The biological record.
        pubmed_record (dict): The PubMed record object from Entrez.
        logger (logging.Logger, optional): Logger instance for logging errors and warnings.

    Returns:
        dict: Metadata including publication title, first author, and publication date.
    """
    pubmed_info = {}
    try:
        ref = record.annotations['references'][0]

        pubmed_info["pubmed_id"] = ref.pubmed_id

        if pubmed_record:
            pubmed_info |= pubmed_info_fill(ref.pubmed_id, pubmed_record)
        else:
            # if logger:
            #     logger.info(f"No pubmed id for {record.id}, fetching author from ncbi.")

            pubmed_info["pub_title"] = ref.title
            if ref.journal:
                pubmed_info["pub_date"] = extract_date(ref.journal)
            if ref.authors:
                first_aut = extract_first_author(ref.authors)
                # TODO for testing
                if first_aut[:3] != ref.authors[:3]:
                    if logger:
                        logger.warning(f"Author extraction seems to have failed for ID: {record.id}")
                pubmed_info["first_aut"] = extract_first_author(ref.authors)

    except Exception as e:
        if logger:
            logger.error(f"Error filtering entry for ID {record.id}: {e}")

    return pubmed_info


def pubmed_info_fill(pubmed_id, pubmed_record, logger=None):
    """
    Fill in PubMed metadata fields using a PubMed record object.

    Args:
        pubmed_id (str): PubMed ID of the record.
        pubmed_record (dict): The PubMed record object from Entrez.
        logger (logging.Logger, optional): Logger instance for logging errors.

    Returns:
        dict: Metadata including title, date, and first author.
    """
    pubmed_api_info = {}
    try:
        pubmed_api_info["pub_title"] = pubmed_record[0]["Title"]
        pubmed_api_info["pub_date"] = pubmed_record[0]["PubDate"]
        pubmed_api_info["first_aut"] = pubmed_record[0]["AuthorList"][0]

    except Exception as e:
        if logger:
            logger.error(f"Error gathering pubmed info for ID {pubmed_id}: {e}")

    return pubmed_api_info


def get_geo_info(record, logger):
    """
    Extract geographic origin metadata from a record.

    Args:
        record (Bio.SeqRecord): The biological record.
        logger (logging.Logger): Logger instance for logging errors and warnings.

    Returns:
        dict: Geographic origin metadata, including the geo_loc_name field.
    """
    geo_info = {}

    try:
        if 'geo_loc_name' in record.features[0].qualifiers.keys():
            if len(record.features[0].qualifiers['geo_loc_name']) > 1:
                logger.warning(f"geo_loc_name list longer than 1 element for ID {record.id}")
            geo_info["geo_origin"] = record.features[0].qualifiers['geo_loc_name'][0]
        elif 'origin' in record.features[0].qualifiers.keys():
            # TODO testing
            logger.info(f"Field origin populated for {record.id}")
            if len(record.features[0].qualifiers['origin']) > 1:
                logger.warning(f"origin list longer than 1 element for ID {record.id}")
            geo_info["geo_origin"] = record.features[0].qualifiers['origin'][0]
        elif 'note' in record.features[0].qualifiers.keys():
            if len(record.features[0].qualifiers['note']) > 1:
                logger.warning(f"origin list longer than 1 element for ID {record.id}")
            note = record.features[0].qualifiers['note']
            if isinstance(note, list):
                note = " ".join(note)
            match = re.search(r'origin_locality:([^;]+)', note)
            if match:
                 geo_info["geo_origin"] = match.group(1).strip()

    except Exception as e:
        logger.error(f"Error fetching geo info for ID {record.id}: {e}")

    return geo_info


def get_assembly_info(record, logger=None):#
    """
    Extract assembly-related metadata from a record.

    Args:
        record (Bio.SeqRecord): The biological record.
        logger (logging.Logger, optional): Logger instance for logging errors.

    Returns:
        dict: Metadata including assembly method and sequencing technology.
    """
    assembly_info = {}

    try:
        struc_comment = record.annotations['structured_comment']
        if "Assembly-Data" in struc_comment.keys():
            asm_dict = struc_comment["Assembly-Data"]
        elif "Genome-Assembly-Data" in struc_comment.keys():
            asm_dict = struc_comment["Genome-Assembly-Data"]
        else:
            # TODO are there yet more differing naming schemes?
            raise KeyError

        try:
            assembly_info["asm_method"] = asm_dict['Assembly Method']
        except KeyError:
            pass

        try:
            assembly_info["seq_tech"] = asm_dict['Sequencing Technology']
        except KeyError:
            pass

    except Exception as e:
        if logger:
            logger.error(f"Error fetching assembly info for ID {record.id}: {e}")

    return assembly_info


def extract_date(journal_field):
    """
    Extract a date from a journal field string, if present.

    Args:
        journal_field (str): String containing potential date information.

    Returns:
        str: Extracted date in the format 'mm-yyyy' or 'yyyy', or None if no valid date is found.
    """
    # match various date formats that were found manually
    patterns = [
        r'\b(\d{2}-\d{2}-\d{4})(\b|$)',  # 'dd-mm-yyyy'
        r'\b(\d{2},\d{4})(\b|$)',  # 'mm,yyyy'
        r'\b([A-Za-z]+,\d{4})(\b|$)',  # 'month,yyyy'
        r'\b([A-Za-z]+ \d{4})(\b|$)',  # 'month yyyy'
        r'\b\((\d{4})\)(\b|$)',  # '(yyyy)'
        r'\b(\d{2}-[A-Za-z]{3}-\d{4})(\b|$)'  # 'dd-MMM-yyyy' (e.g., '04-OCT-2024')
    ]

    def is_valid_year(year):
        """Helper function to check if the year is valid
         (between 1900 and the current year)."""
        current_year = datetime.now().year
        return 1900 <= year <= current_year

    for pattern in patterns:
        match = re.search(pattern, journal_field)
        if match:
            date_str = match.group(0)
            try:
                # check matches are valid dates and return as mm-yyyy or yyyy
                if '-' in date_str:
                    if re.match(r'\d{2}-[A-Za-z]{3}-\d{4}', date_str):
                        date_obj = datetime.strptime(date_str, '%d-%b-%Y')
                    elif re.match(r'\d{2}-\d{2}-\d{4}', date_str):
                        date_obj = datetime.strptime(date_str, '%d-%m-%Y')
                    else:
                        continue
                    if is_valid_year(date_obj.year):
                        return date_obj.strftime('%m-%Y')
                elif ',' in date_str:
                    if date_str[0].isdigit():
                        date_obj = datetime.strptime(date_str, '%m,%Y')
                    else:
                        date_obj = datetime.strptime(date_str, '%B,%Y')
                    if is_valid_year(date_obj.year):
                        return date_obj.strftime('%m-%Y')
                elif date_str.isdigit() and len(date_str) == 4:
                    year = int(date_str)
                    if is_valid_year(year):
                        return str(year)
                else:
                    date_obj = datetime.strptime(date_str, '%B %Y')
                    if is_valid_year(date_obj.year):
                        return date_obj.strftime('%m-%Y')

            except ValueError:
                continue

    return None


# trying to catch different author versions
def extract_first_author(authors):
    """
    Extract the first author's name from an authors string and standardize its format.

    Args:
        authors (str): String of authors separated by commas or 'and'.

    Returns:
        str: First author's name in the format 'Lastname, Initials.', or None if no author is found.
    """
    # split by ", " or " and "
    author_list = re.split(r'\s+and\s+|,\s(?=\S)', authors)

    if author_list:
        first_author = author_list[0].strip()
        # parse name to unify format
        if ',' in first_author:
            # assume format: "lastname, initials"
            lastname, rest = first_author.split(',', 1)
            lastname = lastname.strip()
            rest = rest.strip()
            return f"{lastname}, {rest}"
        else:
            # assume "firstname lastname" format
            name_parts = first_author.split()
            if len(name_parts) > 1:
                lastname = name_parts[-1]
                initials = " ".join([name[0] + '.' for name in name_parts[:-1]])
                return f"{lastname}, {initials}"
            else:
                # assume 'lastname'
                return first_author.strip()

    return None
