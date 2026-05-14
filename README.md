# BarTender
## Installation
1. Download and unarchive, enter the directory:  
   ```
    cd BarTender
   ```
2. Install:  
   ```
    pip install .
   ```
## Usage
```
    barcode_db        Create a barcode database necessary for demultiplexing
    demux             Generates files with demultiplexed and trimmed reads,
                      each corresponding to a sample ID.
```
### barcode_db
```
usage: bartender barcode_db [-h] -b BARCODES -o OUTDIR

options:
  -h, --help            show this help message and exit
  -b BARCODES, --barcodes BARCODES
                        Path to barcode FASTA file. File must contain
                        sequences of barcodes WITH adapters, with barcodes
                        facing left. Headers should begin with an F for
                        forward barcodes, with an R for reverse barcodes
  -o OUTDIR, --outdir OUTDIR
                        Output directory.
```
### demux
```
usage: bartender demux [-h] -i INPUT_FASTQ -m MANIFEST -d DBDIR -o OUTDIR
                       [-f FORWARD] [-r REVERSE]
                       [--max_errors_adapter MAX_ERRORS_ADAPTER]
                       [--max_errors_barcode MAX_ERRORS_BARCODE]
                       [--indel_buffer INDEL_BUFFER]

options:
  -h, --help            show this help message and exit
  -i INPUT_FASTQ, --input_fastq INPUT_FASTQ
                        Path to FASTQ file with reads.
  -m MANIFEST, --manifest MANIFEST
                        Path to the manifest - a 3-column (Barcode1, Barcode2,
                        SampleID) CSV where sets of barcodes are matched with
                        corresponding sample IDs
  -d DBDIR, --dbdir DBDIR
                        Barcode database directory (database can be generated
                        with barcode_db).
  -o OUTDIR, --outdir OUTDIR
                        Output directory.
  -f FORWARD, --forward FORWARD
                        Forward barcoding adapter sequence.
  -r REVERSE, --reverse REVERSE
                        Reverse barcoding adapter sequence.
  --max_errors_adapter MAX_ERRORS_ADAPTER
                        Max errors in adapter, default: 5
  --max_errors_barcode MAX_ERRORS_BARCODE
                        Max errors in barcode, default: 2
  --indel_buffer INDEL_BUFFER
                        Max distance of barcode from adapter, default: 5
```
