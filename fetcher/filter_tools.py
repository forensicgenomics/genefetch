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


def read_exclusion_files():
    """
    Read all of the profiles to be excluded from all files in EXCLUSIONS_DIR.

    Returns:
        {"basename_without_ext": set_of_ids, ...} ignores blank lines and lines starting with '#'
    """
    out = {}
    if not os.path.isdir(EXCLUSIONS_DIR):
        return out
    for fname in os.listdir(EXCLUSIONS_DIR):
        if not fname.endswith(".txt"):
            continue
        reason = os.path.splitext(fname)[0]
        path = os.path.join(EXCLUSIONS_DIR, fname)
        ids = set()
        with open(path, "r", encoding="utf-8") as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith("#"):
                    continue
                ids.add(s)
        if ids:
            out[reason] = ids
    return out


def build_exclusion_filters(exclusions_dir=None, id_map=None):
    """
    if id_map is provided, use it; else read files from exclusions_dir.
    returns: [{"description": <filebase>, "fun": lambda rec: rec.id in ids}, ...]
    """
    """
    If id_map is given, use it, otherwise read each .txt file in EXCLUSIONS_DIR.
    Creates a dict with:
        {
            "description": filename_without_ext,
            "fun": lambda record: record.id in <loaded_id_list>
        }

    Args:
        exclusions_dir (str): Path to the directory containing exclusion files.

    Returns:
        list[dict]: A list of filter dictionaries with `description` and `fun` keys.
    """
    if id_map is None:
        if exclusions_dir is None:
            raise ValueError("either id_map or exclusions_dir must be provided")
        id_map = read_exclusion_files(exclusions_dir)

    filters = []
    for desc, ids in id_map.items():
        filters.append({
            "description": desc,
            "fun": (lambda record, ids=ids: record.id in ids)
        })
    return filters


def load_filters(exclusion_id_map=None):
    """
    Load and combine static and dynamically created filters.

    Combines:
        - Predefined static filters (`FILTERS`) from `global_defaults.py`.
        - Dynamically generated filters from exclusion files in `EXCLUSIONS_DIR`.

    Returns:
        list[dict]: A combined list of static and dynamically generated filters.
    """
    dynamic_filters = build_exclusion_filters(exclusions_dir=EXCLUSIONS_DIR, id_map=exclusion_id_map)
    all_filters = FILTERS + dynamic_filters
    return all_filters
