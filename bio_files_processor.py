import argparse
import os
from typing import Iterator, Tuple, Optional


def safe_derive_output_path(input_path: str,
                            suffix: str = ".oneline.fasta") -> str:
    """
    Derive a safe output filename next to the input.
    Examples:
      "reads.fasta" -> "reads.oneline.fasta"
      "reads.fna"   -> "reads.oneline.fasta"
      "weirdfile"   -> "weirdfile.oneline.fasta"
    """
    base_name, extension = os.path.splitext(input_path)
    fasta_like_extensions = {".fa", ".fasta", ".fna", ".faa"}
    if extension.lower() in fasta_like_extensions:
        return f"{base_name}{suffix}"
    return f"{input_path}{suffix}"


def open_for_write_safely(path_to_write: str, allow_overwrite: bool = False):
    """
    Open a file for writing in a safe way:
    - refuses to overwrite an existing file unless explicitly allowed;
    - creates missing parent directories when needed.
    """
    if os.path.exists(path_to_write) and not allow_overwrite:
        raise FileExistsError(
            f"Output file already exists: {path_to_write}. "
            f"Use --force to allow overwriting."
        )

    parent_dir = os.path.dirname(path_to_write)
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)

    return open(path_to_write, "w", encoding="utf-8", newline="\n")


def iter_fasta_records(fasta_path: str) -> Iterator[Tuple[str, str]]:
    """
    Stream FASTA records as (header, sequence) pairs.

    - Lines starting with '>' are treated as headers and kept verbatim.
    - Sequence lines are concatenated until the next header.
    - All whitespace characters inside sequence lines are removed.
    - Empty lines are ignored.
    """
    with open(fasta_path, "r", encoding="utf-8") as file_handle:
        current_header: Optional[str] = None
        sequence_chunks: list[str] = []

        for raw_line in file_handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(sequence_chunks)
                current_header = line
                sequence_chunks = []
            else:
                cleaned_line = "".join(ch for ch in line if not ch.isspace())
                sequence_chunks.append(cleaned_line)

        if current_header is not None:
            yield current_header, "".join(sequence_chunks)


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: Optional[str] = None,
    *,
    allow_overwrite: bool = False,
    keep_case: bool = False
) -> str:
    """
    Convert a multi-line FASTA into a one-line-per-sequence FASTA.
    Each record in the output has exactly two lines:
      1) header line (verbatim)
      2) single sequence line (concatenated)

    Parameters
    ----------
    input_fasta : str
        Path to the input FASTA file.
    output_fasta : Optional[str]
        Path to the output FASTA file. If not provided, a name will be derived
        next to the input: e.g., "reads.fasta" -> "reads.oneline.fasta".
    allow_overwrite : bool
        If True, allow overwriting an existing output file.
    keep_case : bool
        If True, preserve case of sequence letters. Otherwise, convert to upper case.

    Returns
    -------
    str
        Path to the created output file.
    """
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input file not found: {input_fasta}")

    target_path = output_fasta or safe_derive_output_path(input_fasta)

    with open_for_write_safely(target_path, allow_overwrite=allow_overwrite) as out_handle:
        for header, sequence in iter_fasta_records(input_fasta):
            sequence_out = sequence if keep_case else sequence.upper()
            out_handle.write(f"{header}\n{sequence_out}\n")

    return target_path


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bio_files_processor.py",
        description="Utilities for processing bioinformatics files.",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        metavar="command",
        help="Available subcommands. Use 'convert-fasta' to convert FASTA to one-line sequences."
    )

    convert_parser = subparsers.add_parser(
        "convert-fasta",
        help="Convert multi-line FASTA records to one-line sequences.",
        description="Concatenates sequence lines of each FASTA record into a single sequence line."
    )
    convert_parser.add_argument(
        "--input", "-i",
        dest="input_fasta",
        required=True,
        help="Path to the input FASTA file."
    )
    convert_parser.add_argument(
        "--output", "-o",
        dest="output_fasta",
        default=None,
        help="Path to the output FASTA file. If omitted, a name will be derived automatically."
    )
    convert_parser.add_argument(
        "--force",
        dest="allow_overwrite",
        action="store_true",
        help="Allow overwriting the output file if it already exists."
    )
    convert_parser.add_argument(
        "--keep-case",
        dest="keep_case",
        action="store_true",
        help="Preserve original case of sequence letters (default: convert to upper case)."
    )

    return parser


def main() -> None:
    parser = build_cli_parser()
    args = parser.parse_args()

    if args.command == "convert-fasta":
        try:
            output_path = convert_multiline_fasta_to_oneline(
                input_fasta=args.input_fasta,
                output_fasta=args.output_fasta,
                allow_overwrite=args.allow_overwrite,
                keep_case=args.keep_case,
            )
            print(f"[OK] Saved to: {output_path}")
        except FileNotFoundError as error:
            print(f"[Error] {error}")
            raise SystemExit(2)
        except FileExistsError as error:
            print(f"[Error] {error}")
            raise SystemExit(3)


if __name__ == "__main__":
    main()
