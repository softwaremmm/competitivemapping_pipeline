"""Competitive Mapping"""

import sys

from competitivemapping.cli_args import Arguments


def cli_entry_point() -> None:
    """CLI entry point."""
    Arguments(sys.argv[1:])
