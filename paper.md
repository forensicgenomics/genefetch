---
title: "genefetch: automating retrieval and inspection of mitochondrial DNA data from GenBank"
authors:
  - name: Noah Hurmer
    affiliation: 1
    orcid: 0009-0006-8574-9011
  - name: Nicole Huber
    affiliation: 1
    orcid: 0000-0001-9267-6910
affiliations:
  - name: Institute of Legal Medicine, Medical University of Innsbruck, Innsbruck, Austria
    index: 1
date: 2025-10-23
bibliography: paper.bib
---

# Summary

GenBank hosts a continuously expanding collection of mitochondrial DNA (mtDNA) sequences that form a cornerstone for phylogenetic, population, and forensic genetics research. Regularly incorporating new records is essential to maintain up-to-date datasets, yet repeated manual downloads or ad-hoc scripts can make it difficult to reproduce earlier results. Small changes in search terms or filters can alter the retrieved data and obscure the provenance of analyses.

genefetch is a lightweight, open-source Python tool that automates the retrieval of mtDNA sequences and associated metadata—such as geographic origin and sequencing technology—from GenBank. It allows users to define search parameters once, re-execute identical queries over time, and manage persistent exclusion lists to omit undesired samples. By combining automation with transparent configuration files and structured output, genefetch provides a reproducible and easily extendable workflow for assembling well-documented mtDNA datasets suitable for phylogenetic and forensic genomics applications.

# Statement of need

Phylogenetic analyses often rely on periodic downloads of mtDNA data from GenBank using custom scripts or manual queries. When these procedures are repeated over time, even small changes in query parameters or filtering can lead to inconsistent datasets and unclear provenance.

genefetch offers a reproducible framework for these tasks by storing search parameters, exclusion lists, and metadata extraction rules in a single configuration. This approach makes it straightforward to update datasets when new sequences appear and ensures that both sequence and metadata retrieval are performed in a consistent, documented manner.

# Functionality

genefetch is written in **Python 3.8+**, distributed under the **MPL-2.0** license, and hosted at [https://github.com/forensicgenomics/genefetch](https://github.com/forensicgenomics/genefetch).

**Main features**
- Command-line and API interfaces;  
- Parallel or sequential downloading using the NCBI Entrez API;  
- Metadata parsing and export to CSV/TSV;  
- Configurable exclusion and update lists;  
- Logging and progress tracking;  
- Example configuration and test datasets included.

**Example usage**
```bash
python run_fetcher.py \
  --search-term "mitochondrion complete genome[Title] AND Homo sapiens[Organism]" \
  --fetch-parallel \
  --max-num 70000 \
  --batch-size 250 \
  --clean-dir

