
"""Main script"""

import argparse
import sys
from time import perf_counter

import parameters
from lemmas import bernstein_bounds, set_verbose, substituting_estimates


def print_time(total_seconds):
    hours = int(total_seconds // 3600)
    minutes_left = total_seconds - 3600 * hours
    minutes = int(minutes_left // 60)
    seconds_left = minutes_left - 60 * minutes
    out = ""
    if hours > 0:
        out += f" {hours} hours,"
    if minutes > 0:
        out += f" {minutes} minutes,"
    out += f" {seconds_left:.2f} seconds."
    return out


def print_lemma(verified, lemma_label, bounds, elapsed):
    print("")
    status = "OK" if verified else "FAIL"
    print(f"{lemma_label}: {status}.")
    print(f"Execution time:{print_time(elapsed)}")
    if bounds is not None and parameters.VERBOSE:
        print("Used bounds were:")
        for bound, label, side_txt in bounds:
            if side_txt == "apprx":
                print(f" The number {label}: {bound}.")
            else:
                print(f" The {side_txt} bound {label}: {bound}.")
    print("", flush=True)


def print_subst(verified, description):
    if verified:
        print("Every substitution is OK.")
    else:
        print(description)


def run_block(block):
    failures = []
    for lemma in block:
        start = perf_counter()
        verified, lemma_label, bounds = lemma()
        elapsed = perf_counter() - start
        if not verified:
            failures.append(lemma_label)
        print_lemma(verified, lemma_label, bounds, elapsed)
    return failures


def run_substitutions():
    verified, description = substituting_estimates()
    print_subst(verified, description)
    return verified


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reproduce the interval-arithmetic checks for Table bb."
    )
    parser.add_argument(
        "-v", "--verbose", action="count", default=0,
        help=(
            "increase printed information. Use -v for row-level progress and "
            "the bounds read from bounds.txt; use -vv to also print every "
            "accepted subinterval and its local margins."
        ),
    )
    parser.add_argument(
        "--print-subintervals", action="store_true",
        help="print every accepted subinterval and its lower/upper relative margins",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    verbose_level = max(parameters.VERBOSE, int(args.verbose or 0))
    print_subintervals = bool(args.print_subintervals or verbose_level >= 2)
    set_verbose(verbose_level, print_subintervals=print_subintervals)

    blocks = [
        [bernstein_bounds],
    ]
    failures = []
    for block in blocks:
        failures.extend(run_block(block))

    if not run_substitutions():
        failures.append("substituting_estimates")

    if failures:
        print("FAILED:", failures)
        sys.exit(1)
    print(r"The bounds reported in Table~\ref{table:bb} have been verified.")


if __name__ == "__main__":
    main()
