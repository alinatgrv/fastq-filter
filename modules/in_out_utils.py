from pathlib import Path
from typing import Generator, Tuple, Union
import gzip


FastqRecord = Tuple[str, str, str, str]


def open_text_maybe_gzip(file_path: Union[str, Path], mode: str):
    """
    Open a plain text file or a .gz file in text mode.

    Parameters
    ----------
    file_path : Union[str, Path]
        Path to the file to open.
    mode : str
        Either "rt" for reading text or "wt" for writing text.

    Returns
    -------
    IO
        A file-like object opened in the requested mode.
    """
    path_obj = Path(file_path)
    if path_obj.suffix == ".gz":
        return gzip.open(path_obj, mode, encoding="utf-8", newline="")
    return open(path_obj, mode, encoding="utf-8", newline="")


def iterate_fastq_records(input_path: Union[str, Path]) -> Generator[FastqRecord, None, None]:
    """
    Stream FASTQ (or FASTQ.GZ) records one by one, validating structure.

    Yields tuples without trailing newlines: (header, sequence, plus_line, qualities).

    Parameters
    ----------
    input_path : Union[str, Path]
        Path to an input FASTQ or FASTQ.GZ file.

    Yields
    ------
    FastqRecord
        A tuple of four strings representing one FASTQ record.

    Raises
    ------
    ValueError
        If the file is truncated or does not follow FASTQ structure.
    """
    path_obj = Path(input_path)
    line_counter = 0

    with open_text_maybe_gzip(path_obj, "rt") as input_handle:
        while True:
            header_line = input_handle.readline()
            if not header_line:
                break

            sequence_line = input_handle.readline()
            plus_line = input_handle.readline()
            qualities_line = input_handle.readline()
            line_counter += 4

            if not (sequence_line and plus_line and qualities_line):
                raise ValueError(
                    f"Truncated FASTQ: '{path_obj}' ended mid-record "
                    f"(last read ended at line {line_counter})."
                )

            header_line = header_line.rstrip("\n")
            sequence_line = sequence_line.rstrip("\n")
            plus_line = plus_line.rstrip("\n")
            qualities_line = qualities_line.rstrip("\n")

            if not header_line.startswith("@"):
                raise ValueError(
                    f"Invalid FASTQ: header line does not start with '@' "
                    f"(around line {line_counter - 3} in '{path_obj}')."
                )
            if not plus_line.startswith("+"):
                raise ValueError(
                    f"Invalid FASTQ: third line of a record does not start with '+' "
                    f"(around line {line_counter - 1} in '{path_obj}')."
                )
            if len(sequence_line) != len(qualities_line):
                raise ValueError(
                    "Invalid FASTQ: sequence and quality lengths differ "
                    f"(lines {line_counter - 2} and {line_counter} in '{path_obj}')."
                )

            yield header_line, sequence_line, plus_line, qualities_line


def ensure_output_directory(directory_path: Union[str, Path] = "filtered") -> Path:
    """
    Create an output directory if it does not exist.

    Parameters
    ----------
    directory_path : Union[str, Path], optional
        Directory to create or ensure, by default "filtered".

    Returns
    -------
    Path
        The ensured directory path.
    """
    output_dir = Path(directory_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def build_safe_output_path(
    file_name: str,
    directory: Union[str, Path] = "filtered",
    force: bool = False,
) -> Path:
    """
    Build a safe output path that avoids clobbering existing files.

    If a file already exists and force=False, a numeric suffix like _1, _2, ...
    is appended to the stem until a free name is found.

    Parameters
    ----------
    file_name : str
        Desired output file name (e.g., "out.fastq" or "out.fastq.gz").
    directory : Union[str, Path], optional
        Output directory, by default "filtered".
    force : bool, optional
        If True, allow overwriting existing files, by default False.

    Returns
    -------
    Path
        A safe, non-existing output path (unless force=True).
    """
    output_dir = ensure_output_directory(directory)
    candidate_path = output_dir / file_name

    if force:
        return candidate_path
    if not candidate_path.exists():
        return candidate_path

    suffix_index = 1
    base_name = candidate_path.stem
    extension = candidate_path.suffix

    while True:
        next_candidate = output_dir / f"{base_name}_{suffix_index}{extension}"
        if not next_candidate.exists():
            return next_candidate
        suffix_index += 1


def write_fastq_record(output_handle, record: FastqRecord) -> None:
    """
    Write a single FASTQ record (4 lines) to an open file handle.

    Parameters
    ----------
    output_handle : IO
        Open file-like object in text write mode.
    record : FastqRecord
        Tuple (header, sequence, plus_line, qualities).
    """
    header_line, sequence_line, plus_line, qualities_line = record
    output_handle.write(f"{header_line}\n")
    output_handle.write(f"{sequence_line}\n")
    output_handle.write(f"{plus_line}\n")
    output_handle.write(f"{qualities_line}\n")


def open_output_fastq(output_path: Union[str, Path]):
    """
    Open an output FASTQ (or .gz) in text mode for writing.

    The caller is responsible for closing the returned handle.

    Parameters
    ----------
    output_path : Union[str, Path]
        Path to the output file.

    Returns
    -------
    IO
        A file-like object opened in text write mode.
    """
    return open_text_maybe_gzip(output_path, "wt")
