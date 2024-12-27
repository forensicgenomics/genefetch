
import os

from .global_defaults import EXCLUSIONS_DIR, FILTERS


def build_exclusion_filters(exclusions_dir):
    """
    For each .txt file in EXCLUSIONS_DIR, create a dict with:
        {
            "description": filename_without_ext,
            "fun": lambda record: record.id in <loaded_id_list>
        }
    Returns a list of these dicts.
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


# combine static + dynamic filters:
def load_filters():
    dynamic_filters = build_exclusion_filters(exclusions_dir=EXCLUSIONS_DIR)
    all_filters = FILTERS + dynamic_filters
    return all_filters
