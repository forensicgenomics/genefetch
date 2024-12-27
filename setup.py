from setuptools import setup, find_packages

setup(
    name="fetcher",
    version="1.0.0",
    description="A module for fetching genetic profiles, their metadata and sequences from the NCBI db.",
    author="Noah Hurmer",
    author_email="minimops.projects@gmail.com",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "biopython",
        "ratelimit",
    ],
    entry_points={
        "console_scripts": [
            "fetch=fetch.fetch:main",
        ]
    },
)
