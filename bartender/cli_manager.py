#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys


def main():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    BASH_SCRIPT_PATH = os.path.join(package_dir, "barcode_db.sh")
    PYTHON_PIPELINE_PATH = os.path.join(package_dir, "demux.py")

    parser = argparse.ArgumentParser(
        description="Bartender — a sensitive demultiplexing instrument for long reads"
    )
    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Available commands"
    )

    # Subcommand 1: barcode_db 
    parser_db = subparsers.add_parser(
        "barcode_db", help="Create a barcode database necessary for demultiplexing"
    )
    parser_db.add_argument(
        "-b",
        "--barcodes",
        required=True,
        help="Path to barcode FASTA file. File must contain sequences of barcodes WITH adapters, with barcodes facing left. Headers should begin with an F for forward barcodes, with an R for reverse barcodes",
    )
    parser_db.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory.",
    )

    # Subcommand 2: demux
    parser_demux = subparsers.add_parser(
        "demux", help="Generates files with demultiplexed and trimmed reads, each corresponding to a sample ID."
    )
    parser_demux.add_argument(
        "-i", "--input_fastq", required=True, help="Path to FASTQ file with reads."
    )
    parser_demux.add_argument(
        "-m", "--manifest", required=True, help="Path to the manifest - a 3-column (Barcode1, Barcode2, SampleID) CSV where sets of barcodes are matched with corresponding sample IDs"
    )
    parser_demux.add_argument(
        "-d",
        "--dbdir",
        required=True,
        help="Barcode database directory (database can be generated with barcode_db).",
    )
    parser_demux.add_argument(
        "-o", "--outdir", required=True, help="Output directory."
    )
    parser_demux.add_argument(
        "--max_errors_adapter",
        type=int,
        default=5,
        help="Max errors in adapter, default: 5",
    )
    parser_demux.add_argument(
        "--max_errors_barcode",
        type=int,
        default=2,
        help="Max errors in barcode, default: 2",
    )
    parser_demux.add_argument(
        "--indel_buffer",
        type=int,
        default=5,
        help="Max distance of barcode from adapter, default: 5",
    )

    args = parser.parse_args()

    # LOGIC

    if args.command == "barcode_db":
        outdir = os.path.abspath(args.outdir)
        os.makedirs(outdir, exist_ok=True)

        print("=== [Bartender] Barcode database generation ===")
        shutil.copy(
            os.path.abspath(args.barcodes),
            os.path.join(outdir, "input_barcodes.fasta"),
        )

        try:
            subprocess.run(
                ["bash", BASH_SCRIPT_PATH, "input_barcodes.fasta"],
                cwd=outdir,
                check=True,
            )
            if os.path.exists(os.path.join(outdir, "input_barcodes.fasta")):
                os.remove(os.path.join(outdir, "input_barcodes.fasta"))
            print(f"-> Database created successfully in {outdir}")
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    elif args.command == "demux":
        print("=== [Bartender] Demultiplexing ===")
        outdir = os.path.abspath(args.outdir)
        os.makedirs(outdir, exist_ok=True)


        py_cmd = [
            "python3",
            PYTHON_PIPELINE_PATH,
            "-i",
            os.path.abspath(args.input_fastq),
            "-m",
            os.path.abspath(args.manifest),
            "-d",
            os.path.abspath(args.dbdir),
            "-o",
            outdir,
            "--max_errors_adapter",
            str(args.max_errors_adapter),
            "--max_errors_barcode",
            str(args.max_errors_barcode),
            "--indel_buffer",
            str(args.indel_buffer),
        ]

        try:
            subprocess.run(py_cmd, check=True)
            print(f"-> Finished demultiplexing. Results in {outdir}")
            print("\n=== [Bartender] Генерация графиков статистики ===")
            from demux_package.plots_generator import generate_demux_plots
            
            report_path = os.path.join(outdir, "demux_report.tsv")
            generate_demux_plots(report_path, outdir)
        except subprocess.CalledProcessError as e:
            print(f"Pipeline error: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
