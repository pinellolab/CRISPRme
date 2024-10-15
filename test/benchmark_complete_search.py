"""
"""

from typing import List, Tuple, Dict, Set
from time import time

import sys
import os


def parse_commandline(args: List[str]) -> Tuple[str, str]:
    # TODO: fix benchmark to point zenodo's NGG 1000G+HGDP dataset
    if len(args) != 2:
        raise ValueError("Too many/few input arguments")
    benchmark_fname, test_fname = args
    if not os.path.isfile(benchmark_fname):
        raise FileNotFoundError(f"Unable to find {benchmark_fname}")
    if not os.path.isfile(test_fname):
        raise FileNotFoundError(f"Unable to find {test_fname}")
    return benchmark_fname, test_fname


def load_dataset(dataset_fname: str) -> List[List[str]]:
    try:  # load the input dataset and split each line in individual fields
        with open(dataset_fname, mode="r") as infile:
            header = [infile.readline().strip().split()]  # read header
            data = [line.strip().split() for line in infile]
            return header + data  # append data content to header line
    except FileNotFoundError as e:  # catch file existence errors
        raise FileNotFoundError(f"Unable to find {dataset_fname}") from e
    except PermissionError as e:  # catch permission errors
        raise PermissionError(f"Permission denied: cannot open {dataset_fname}") from e
    except OSError as e:  # catch os-related errors
        raise OSError(f"An error occurred while accessing {dataset_fname}") from e
    except Exception as e:  # catch other unexpected errors
        # sourcery skip: raise-specific-error
        raise Exception(f"An error occurred while reading {dataset_fname}") from e


def report_header_mismatch(benchmark_data: List[str], test_data: List[str]) -> str:
    # recover mismatching fields to report
    # the report will consists of a list of tuples like
    # (2, benchmark[2] = "something" - test[2] = "something else")
    mmfields = [
        f'benchmark[{i}] = "{e}" - test[{i}] = "{test_data[i]}"'
        for i, e in enumerate(benchmark_data)
        if e != test_data[i]
    ]
    return "\n".join(mmfields)  # prepare the string to print to stderr


def check_headers(benchmark_header: List[str], test_header: List[str]) -> bool:
    """
    Compares two lists of headers to determine if they match exactly in both
    content and position.
    Returns True if there is any mismatch between the benchmark and test headers.

    Args:
        benchmark_header (List[str]): The list of benchmark headers to compare
            against.
        test_header (List[str]): The list of test headers to be compared.

    Returns:
        bool: True if there is a mismatch; otherwise, False.
    """

    # compare the benchmark and test headers -> each term must exactly match
    # in terms of position and content
    # if found a mismatch return true -> the parent will throw an exception
    return any(h != test_header[i] for i, h in enumerate(benchmark_header))


def format_values_list(values: str) -> str:
    values_ = values.split(",")  # retrieve individual values
    # sort individual values - old runs may show a different entries order than
    # more recent versions
    values_sorted = sorted(values_, reverse=False)
    return ",".join(values_sorted)


def compare_targets_number(
    benchmark_data: List[List[str]], test_data: List[List[str]]
) -> bool:
    # compare the number of reported targets; the number of reported targets must
    # match between the benchmark and test datasets
    return len(benchmark_data) != len(test_data)


def compare_targets_coords(
    benchmark_hash: Dict[str, Tuple[str, List[str]]],
    test_hash: Dict[str, Tuple[str, List[str]]],
    unreported_b: List[str],
    unreported_t: List[str],
) -> Tuple[List[str], List[str], Set[str]]:
    # check if there are mismatches in the targets coordinates
    benchmark_coords = set(benchmark_hash.keys())
    test_coords = set(test_hash.keys())
    mismatching_coords_b = benchmark_coords.difference(test_coords)
    mismatching_coords_t = test_coords.difference(benchmark_coords)
    if mismatching_coords_b:  # store targets reported in benchmark but not in test
        unreported_b += [benchmark_hash[k][0] for k in mismatching_coords_b]
    if mismatching_coords_t:  # store targets reported in test but not in benchmark
        unreported_t += [test_hash[k][0] for k in mismatching_coords_t]
    common_coords = benchmark_coords.intersection(
        test_coords
    )  # store common coordinates
    return unreported_b, unreported_t, common_coords


def hash_data(
    data: List[List[str]], idxstart: int, idxstop: int
) -> Dict[str, Tuple[str, List[str]]]:
    data_hash = {
        f"{line[0]}_{line[1]}_{line[idxstart]}": (
            "\t".join(
                line[:2]
                + [
                    format_values_list(e) if "," in e else e
                    for e in line[idxstart:idxstop]
                ]
            ),
            line[idxstart + 1 : idxstop],
        )
        for line in data
    }
    assert len(data) == len(data_hash)  # we shouldn't lose anything
    return data_hash


def compare_data_content(
    benchmark_hash: Dict[str, Tuple[str, List[str]]],
    test_hash: Dict[str, Tuple[str, List[str]]],
    coords: Set[str],
    unreported_b: List[str],
    unreported_t: List[str],
) -> Tuple[List[str], List[str]]:
    # iterate over target coordinates shared between benchamrk and test
    for k in coords:
        if benchmark_hash[k][0] != test_hash[k][0]:  # mismatch found -> report it
            unreported_b.append(benchmark_hash[k][0])
            unreported_t.append(test_hash[k][0])
    return unreported_b, unreported_t


def write_missing_targets_report(fname: str, header: str, targets: str):
    try:  # write missing targets report
        with open(fname, mode="w") as outfile:
            outfile.write(f"{header}\n{targets}")
    except PermissionError as e:  # catch permission errors
        raise PermissionError(f"Permission denied: cannot write {fname}") from e
    except OSError as e:  # catch os-related errors
        raise OSError(f"An error occurred while reporting targets in {fname}") from e
    except Exception as e:  # catch other unexpected errors
        # sourcery skip: raise-specific-error
        raise Exception(f"An error occurred while writing {fname}") from e


def report_missing_targets(
    unreported_b: List[str],
    unreported_t: List[str],
    header: List[str],
    fprefix: str,
    datatype: str,
) -> None:
    if unreported_b:
        fname = f"{fprefix}_{datatype}_benchmark.tsv"
        write_missing_targets_report(fname, "\t".join(header), "\n".join(unreported_b))
    if unreported_t:
        fname = f"{fprefix}_{datatype}_test.tsv"
        write_missing_targets_report(fname, "\t".join(header), "\n".join(unreported_t))


def compare_data_fields(
    benchamark_data: List[List[str]],
    test_data: List[List[str]],
    header: List[str],
    start: int,
    stop: int,
    datatype: str,
) -> bool:
    # compute hash maps for benchmark and test datasets
    benchmark_hash = hash_data(benchamark_data, start, stop)
    test_hash = hash_data(test_data, start, stop)
    # lists to store unreported targets in benchmark and test
    unreported_b, unreported_t = [], []
    # compare targets coords between benchmark and test datasets
    # recover coordinates shared between benchmark and test datasets
    unreported_b, unreported_t, common_coords = compare_targets_coords(
        benchmark_hash, test_hash, unreported_b, unreported_t
    )
    # compare targets data content
    unreported_b, unreported_t = compare_data_content(
        benchmark_hash, test_hash, common_coords, unreported_b, unreported_t
    )
    # report mismatching targets
    report_missing_targets(
        unreported_b, unreported_t, header, "missing_targets", datatype
    )
    return bool(unreported_b) or bool(unreported_t)


def compare_data(
    benchamark_data: List[List[str]], test_data: List[List[str]], header: List[str]
) -> None:
    # benchmark and test datasets must report the same exact number of targets
    if compare_targets_number(benchamark_data, test_data):
        raise ValueError("Mismatching number of reported targets")
    reportmsg = "Mismatches found in"
    # CFD score columns
    h = header[:2] + header[2:29]
    if compare_data_fields(benchamark_data, test_data, h, 2, 29, "CFD"):
        reportmsg = f"{reportmsg} CFD columns,"
    # fewest mm+b columns
    h = header[:2] + header[29:54]
    if compare_data_fields(benchamark_data, test_data, h, 29, 54, "mm+bulges"):
        reportmsg = f"{reportmsg} mm+bulges columns,"
    # CRISTA score columns
    h = header[:2] + header[54:79]
    if compare_data_fields(benchamark_data, test_data, h, 54, 79, "CRISTA"):
        reportmsg = f"{reportmsg} CRISTA columns,"
    if reportmsg != "Mismatches found in":
        raise ValueError(reportmsg)


def compare(benchmark_fname: str, test_fname: str):
    # load benchmark and test datasets
    benchmark_data = load_dataset(benchmark_fname)
    test_data = load_dataset(test_fname)
    # compare benchmark and test header content; the header content must match
    # by both position and content
    if check_headers(benchmark_data[0], test_data[0]):  # if true found mismatch
        mmfields = report_header_mismatch(benchmark_data[0], test_data[0])
        raise ValueError(f"Mismatches occurred on header:\n{mmfields}")
    header = benchmark_data[0]  # already ensured that headers match
    compare_data(benchmark_data[1:], test_data[1:], header)
    sys.stderr.write("The benchmark and test datasets fully match\n")


def main():
    # read input arguments -> expected arguments benchmark dataset and test
    # dataset to compare
    benchmark_fname, test_fname = parse_commandline(sys.argv[1:])
    start = time()
    sys.stderr.write("Starting dataset comparison\n")
    compare(benchmark_fname, test_fname)
    sys.stderr.write(f"Elapsed time {(time() - start):.2f}s\n")


if __name__ == "__main__":
    main()
