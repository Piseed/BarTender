from setuptools import setup, find_packages

setup(
    name="BarTender",
    version="1.0.0",
    description="CLI tool for soft FASTQ demultiplexing",
    author="",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "edlib",
        "pandas",
        "openpyxl",
        "matplotlib",
        "seaborn"
    ],
    include_package_data=True,
    package_data={
        "bartender": ["barcode_db.sh", "demux.py"],
    },

    entry_points={
        "console_scripts": [
            "bartender=bartender.cli_manager:main",
        ],
    },
)
