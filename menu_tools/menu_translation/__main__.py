"""Translate between MenuTools rate-table menus and CMSSW P2GT menus.

    python -m menu_tools.menu_translation yaml-to-cmssw MENU.yml -o seeds_cff.py
    python -m menu_tools.menu_translation cmssw-to-yaml MENU_DUMP.py -o menu.yml \\
        [--cmssw-src /path/to/CMSSW/src]

The CMSSW input is a structured config dump (dumpL1TGTMenu.py in CMSSW, or a
ConfDB export). Seeds that cannot be expressed in the target format are listed
and the exit code is 1; everything else is written.
"""

import argparse
import os
import sys

from menu_tools.menu_translation.cmssw_menu_reader import read_cmssw_menu
from menu_tools.menu_translation.cmssw_menu_writer import (
    DEFAULT_MENU_PACKAGE,
    write_cmssw_menu,
)
from menu_tools.menu_translation.object_mapping import ObjectMapping
from menu_tools.menu_translation.rate_table_yaml_io import (
    read_rate_table_menu,
    write_rate_table_menu,
)


def _report(skipped: dict[str, str], warnings: list[str]) -> int:
    for warning in warnings:
        print(f"WARNING {warning}", file=sys.stderr)
    for name, reason in skipped.items():
        print(f"SKIPPED {name}: {reason}", file=sys.stderr)
    return 1 if skipped else 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("direction", choices=["yaml-to-cmssw", "cmssw-to-yaml"])
    parser.add_argument("input")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--version", default="V43nano", help="MenuTools config version")
    parser.add_argument(
        "--object-map",
        help="object mapping (default: configs/<version>/cmssw_object_map.yaml)",
    )
    parser.add_argument("--cmssw-src", help="CMSSW source checkout (outside cmsenv)")
    parser.add_argument("--menu-package", default=DEFAULT_MENU_PACKAGE)
    args = parser.parse_args()

    map_path = args.object_map or f"configs/{args.version}/cmssw_object_map.yaml"
    mapping = ObjectMapping.from_file(map_path)
    warnings = [
        f"'{key}' is not defined in configs/{args.version}/objects"
        for key in mapping.undefined_menu_tools_keys(f"configs/{args.version}/objects")
    ]

    if args.direction == "yaml-to-cmssw":
        seeds, skipped, parse_warnings = read_rate_table_menu(
            args.input, mapping.primary_vertex
        )
        warnings += parse_warnings
        text, write_skipped = write_cmssw_menu(
            seeds, mapping, os.path.basename(args.input), args.menu_package
        )
        skipped.update(write_skipped)
        n_written = len(seeds) - len(write_skipped)
    else:
        seeds, skipped = read_cmssw_menu(args.input, mapping, args.cmssw_src)
        text = write_rate_table_menu(seeds, mapping.primary_vertex)
        n_written = len(seeds)

    with open(args.output, "w") as f:
        f.write(text)
    print(f"wrote {n_written} seeds to {args.output}")
    return _report(skipped, warnings)


if __name__ == "__main__":
    sys.exit(main())
