
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
python run_fetcher.py
  --search-term "mitochondrion complete genome AND Homo Sapiens[Organism]"
  --max-num 1000
  --batch-size 100
  --fetch-parallel
  --num-workers 20
  --soft-restart
```

or as a module:   

```bash
python -m fetcher.main
  --search-term "mitochondrion complete genome AND Homo Sapiens[Organism]"
  --max-num 1000
  ...
```

### **Options**
| Argument         | Description                                                |
|------------------|------------------------------------------------------------|
| `--max-num`      | Maximum number of profiles to fetch.                      |
| `--batch-size`   | Number of profiles to process in each batch.              |
| `--fetch-parallel` | Enable parallel fetching.                               |
| `--num-workers`  | Number of workers for parallel fetching.                  |
| `--soft-restart` | Restart with previously processed profiles.               |
| `--search-term`  | NCBI search term for fetching mitochondrial sequences.    |

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
│   ├── processed_ids/          # Processed ID lists
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

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

---

## **Acknowledgments**

genefetch leverages the NCBI Entrez API and Biopython for efficient sequence and metadata management.

genefetch was created as part of the [mitoTree Project](https://genomics.gmi.tirol/projects/mitoTree/).
