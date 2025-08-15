
# genefetch

**genefetch** is a Python-based tool for fetching, processing, and locally saving Profiles from the NCBI database.
Based on a given search term, it fetches profiles, filters them based on given rules, extracts metadata and fasta files
and stores everything into a convenient structure.

---

## **Features**

- Fetch sequences using customizable NCBI search queries.
- Support for parallel and sequential processing.
- Metadata extraction with publication information, assembly and sequencing technology, and geographic information.
- Automated handling of sequence updates and filtering based on custom criteria.
- Post-processing checks for data consistency and integrity.

---

## **Installation**

### **Prerequisites**
- Python 3.8 or higher
- NCBI API key (optional for increased rate limits)

### **Setup**

1. Clone the repository:
   ```bash
   git clone https://github.com/forensicgenomics/genefetch.git
   cd genefetch
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Add your email (and optionally NCBI API key):
   - Save your email in `secrets/ncbi_email.txt`.
   - Save your NCBI API key in `secrets/ncbi_api_key.txt`.

4. Adjust defaults
   - Change global vars in `fetcher/global_defaults.py` if desired. 
   - Add txt files of accession numbers to be excluded into `exclusions/`
---

## **Usage**

You can use genefetch as a command-line tool or integrate it into a larger pipeline.

### **Command-Line Usage**
Run the fetcher with customizable options:

```bash
python run_fetcher.py \
  --search-term "mitochondrion complete genome AND Homo Sapiens[Organism]" \
  --max-num 1000 \
  --batch-size 100 \
  --fetch-parallel \
  --num-workers 20 \
  --soft-restart
```

or as a module:   

```bash
python -m fetcher.main \
  --search-term "mitochondrion complete genome AND Homo Sapiens[Organism]" \
  --max-num 1000 \
  ...
```

### **Options**
| Argument              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--search-term`       | Term to eSearch the ncbi database with for fetching profiles.<br/>See the ncbi nuccore docu on how to adjust this.                                                                                                                                                                                                                                                                                                                                                                  |
| `--max-num`           | Maximum number of profiles to fetch (Maximum of 10k as per Entrez).                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `--batch-size`        | (optional, defaults to 100) Number of profiles to process before saving.                                                                                                                                                                                                                                                                                                                                                                                                            |
| `--num-workers`       | (optional, defaults to 16) Number of workers for parallel fetching.                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `--fetch-parallel`    | (flag) Enable parallel fetching.                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `--soft-restart`      | (flag) Restart with previously processed profiles.<br/>Use this if you are building upon the last run with the same search term and dont want to fetch all profiles at once.<br/>If set, all profiles of the previous fetched are not checked for updates!                                                                                                                                                                                                                          |
| `--clean-dir`         | (flag) Removes profile data from all files (ids_list / removed_ids / metadata / fasta_files) if they are not returned by the databank query with the current search parameters.<br/>As the `max-num` parameter, the `LIMIT_NUM` is used.<br/>Set this flag if you made a mistake with the last fetch search term for example and dont want to start over or manually remove them.<br/>Be careful though, as this can remove a lot of data if your current search term is incorrect! |
| `--update-exclusions` | (flag) Call with this flag if you think that the exclusion dir has changed since the last fetch. It removes any profiles from (ids_list / metadata / fasta_files) if they are new in the exclusons and adds them to removed_ids. Any profiles in removed_ids that are now no longer found in the exclusions dir are removed from the removed file. Be careful though, this currently also removes dynamically filtered out profiles from filters.                                   |
| `--y`                 | (flag) Automatic 'yes' to all user questions. Use when integrating into other scripts / workflows.                                                                                                                                                                                                                                                                                                                                                                                  |


You cannot use `soft-restart` in combination with `clean-dir` or `update-exclusions`.

Default values are stored in `genefetch/global_defaults.py`.

---

## **Integrated Workflow**

1. **Fetch Profiles**: Fetch profiles using the specified search term.
2. **Process Profiles**:
   - Filter out unneeded profiles.
   - Fetch metadata for each sequence.
   - Save processed sequences and metadata.
3. **Post-Processing**:
   - Remove duplicates.
   - Validate metadata and sequence files.
   - Clean up old logs and intermediate files.

---

## **Configuration**

MitoFetch uses a `global_defaults.py` file to define global settings like directories, rate limits, and search terms.
Update this file to customize your fetcher.

If you have any profiles you want to exclude during fetching, you can add any number of `txt` files to `exclusions/`.
Filters for these will be created dynamically.

---

## **Development**

### **Project Structure**
```
genefetch/
│
├── fetcher/
│   ├── __init__.py
│   ├── fetch.py                # Main fetcher script
│   ├── fetch_tools.py          # Helper functions for Entrez queries
│   ├── file_io.py              # File handling utilities
│   ├── metadata_tools.py       # Metadata extraction functions
│   ├── logger_setup.py         # Logger setup
│   ├── global_defaults.py      # Global configuration
│   ├── filter_tools.py         # Dynamic filtering rules
│   └── post_process_check.py   # Post-processing validation
│  
├── exclusions/                 # Manual exclusion lists
├── secrets/                    # NCBI API keys and email 
│ 
│── data/
│   ├── logs/                   # Logs for debugging and auditing
│       └── debug_info          # dditional files that may be useful for debugging
│   ├── processed_ids/          # Processed ID lists for soft restarting
│   ├── seqs/                   # Sequence FASTA files
│   ├── ids_list.txt            # Fetched profile accession numbers
│   ├── removed_ids.csv         # Filtered out profiles
│   ├── metadata.csv            # Metadata df of fetched profiles
│   └── last_run_date.txt       # Date stored when last run occured
│
├── tests/                      # Unit tests
├── requirements.txt            # Python dependencies
├── README.md                   
└── LICENSE                     
```


## **License**

This project is licensed under the MPL 2.0 License. See the [LICENSE](./LICENSE) file for details.

---

## **Acknowledgments**

genefetch leverages the NCBI Entrez API and Biopython for efficient sequence and metadata management.

genefetch was created as part of the [mitoTree Project](https://genomics.gmi.tirol/projects/mitoTree/).
