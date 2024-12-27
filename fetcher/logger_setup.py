# This file is part of the mitoTree project and authored by Noah Hurmer.
#
# Copyright 2024, Noah Hurmer & mitoTree.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
logging_setup.py

This script configures a logger for the mitoFetch pipeline. The logger supports two levels of output:
1. DEBUG logs, saved to a dedicated debug file for detailed tracing and debugging.
2. INFO and higher-level logs, saved to a general log file for operational insights.

The logger is used throughout the project to log events, warnings, and errors.

Author: Noah Hurmer (aka. minimops) as part of the mitoTree Project.
"""

import logging
from .global_defaults import LOG_FILE, DEBUG_LOG_FILE

def get_logger():
    """
    Set up and return a logger instance for the mitoFetch pipeline.

    The logger writes:
    - DEBUG-level and above logs to a debug log file (DEBUG_LOG_FILE).
    - INFO-level and above logs to a general log file (LOG_FILE).

    Returns:
        logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger("mitoFetchLogger")
    logger.setLevel(logging.DEBUG)

    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # handler 1: write all logs (DEBUG+) to debug.log
    debug_handler = logging.FileHandler(DEBUG_LOG_FILE, mode='a')
    debug_handler.setLevel(logging.DEBUG)
    debug_handler.setFormatter(formatter)
    logger.addHandler(debug_handler)

    # handler 2: write only INFO+ logs to LOG_FILE (as specified in global_defaults)
    info_handler = logging.FileHandler(LOG_FILE, mode='a')
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(formatter)
    logger.addHandler(info_handler)

    logger.propagate = False

    return logger