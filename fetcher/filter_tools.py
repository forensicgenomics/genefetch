# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
exclusion_filters.py

This script provides functions to dynamically generate filters for record exclusion
based on user-provided exclusion lists in global EXCLUSIONS_DIR,
as well as combine them with predefined static filters from global_defaults.py.

Functions:
- build_exclusion_filters: Dynamically creates filters from exclusion files.
- load_filters: Combines static filters with dynamically created filters.

Author: Noah Hurmer as part of the mitoTree Project.
"""

import os

from .global_defaults import EXCLUSIONS_DIR, FILTERS


def build_exclusion_filters(exclusions_dir):
    """
    For each .txt file in EXCLUSIONS_DIR, create a dict with:
        {
            "description": filename_without_ext,
            "fun": lambda record: record.id in <loaded_id_list>
        }

    Args:
        exclusions_dir (str): Path to the directory containing exclusion files.

    Returns:
        list[dict]: A list of filter dictionaries with `description` and `fun` keys.
    """
    exclusion_filters = []

    if not os.path.isdir(exclusions_dir):
        return exclusion_filters

    for file_name in os.listdir(exclusions_dir):
        # skip non-txt files
        if not file_name.endswith(".txt"):
            continue

        full_path = os.path.join(exclusions_dir, file_name)
        # e.g. file_name = "manual_exclusion.txt" -> description: "manual_exclusion"
        description = os.path.splitext(file_name)[0]

        with open(full_path, "r") as f:
            exclusion_ids = [line.strip() for line in f if line.strip()]

        # create filter dictionary
        exclusion_filters.append({
            "description": description,
            "fun": lambda record, ids=exclusion_ids: record.id in ids
        })

    return exclusion_filters


def load_filters():
    """
    Load and combine static and dynamically created filters.

    Combines:
        - Predefined static filters (`FILTERS`) from `global_defaults.py`.
        - Dynamically generated filters from exclusion files in `EXCLUSIONS_DIR`.

    Returns:
        list[dict]: A combined list of static and dynamically generated filters.
    """
    dynamic_filters = build_exclusion_filters(exclusions_dir=EXCLUSIONS_DIR)
    all_filters = FILTERS + dynamic_filters
    return all_filters
