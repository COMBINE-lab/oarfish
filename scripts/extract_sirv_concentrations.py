#!/usr/bin/env python3
"""Extract SIRV E0/E1/E2 molar concentrations from Lexogen's XLSX file."""

import argparse
import re
import xml.etree.ElementTree as ET
import zipfile


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"


def column_number(cell_reference):
    letters = re.match(r"[A-Z]+", cell_reference).group(0)
    result = 0
    for letter in letters:
        result = result * 26 + ord(letter) - ord("A") + 1
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook")
    parser.add_argument("output")
    args = parser.parse_args()

    tag = f"{{{MAIN_NS}}}"
    with zipfile.ZipFile(args.workbook) as archive:
        shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        shared = [
            "".join(node.text or "" for node in item.iter(f"{tag}t"))
            for item in shared_root
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    rows = []
    for row in sheet.findall(f".//{tag}row"):
        values = {}
        for cell in row.findall(f"{tag}c"):
            value = cell.find(f"{tag}v")
            if value is None:
                continue
            parsed = value.text
            if cell.get("t") == "s":
                parsed = shared[int(parsed)]
            values[column_number(cell.get("r"))] = parsed
        name = values.get(7, "")
        concentrations = (values.get(11), values.get(12), values.get(13))
        if re.fullmatch(r"SIRV\d{3}", name) and all(concentrations):
            try:
                tuple(float(value) for value in concentrations)
            except ValueError:
                continue
            rows.append((name, *concentrations))

    with open(args.output, "w", encoding="utf-8") as output:
        output.write("transcript\tE0\tE1\tE2\n")
        for name, e0, e1, e2 in rows:
            output.write(f"{name}\t{e0}\t{e1}\t{e2}\n")


if __name__ == "__main__":
    main()
