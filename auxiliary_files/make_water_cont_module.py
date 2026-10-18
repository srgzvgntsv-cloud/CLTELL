import argparse
import re
import sys
import textwrap
from decimal import Decimal
from pathlib import Path
NUMBER = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
)
def read_values(data_section, name, expected_count):
    pattern = re.compile(
        r"(?m)^\s*" + re.escape(name) + r"\s*=\s*([^;]*);"
    )
    matches = list(pattern.finditer(data_section))
    if len(matches) != 1:
        raise ValueError(
            f"{name}: expected one assignment, found {len(matches)}"
        )
    body = matches[0].group(1).strip()
    tokens = [token.strip() for token in body.split(",")]
    if len(tokens) != expected_count:
        raise ValueError(
            f"{name}: expected {expected_count} values, found {len(tokens)}"
        )
    for index, token in enumerate(tokens, start=1):
        if NUMBER.fullmatch(token) is None:
            raise ValueError(
                f"{name}({index}): invalid numeric value {token!r}"
            )
    return tokens
def decimal_value(token):
    return Decimal(token.lower().replace("d", "e"))
def fortran_value(token):
    token = token.lower().replace("e", "d")
    if "d" in token:
        return token
    if "." not in token:
        token += ".0"
    return token + "d0"
def read_attribute(source, name):
    pattern = re.compile(
        r":" + re.escape(name)
        + r'\s*=\s*"((?:\\.|[^"\\])*)"\s*;',
        re.DOTALL,
    )
    match = pattern.search(source)
    if match is None:
        raise ValueError(f"Missing global attribute: {name}")
    return match.group(1).strip()
def append_comment(lines, label, value):
    lines.append(" ! " + label)
    paragraphs = value.splitlines()
    for paragraph in paragraphs:
        wrapped = textwrap.wrap(
            paragraph.strip(),
            width=110,
            break_long_words=False,
            break_on_hyphens=False,
        )
        if not wrapped:
            lines.append(" !")
        else:
            lines.extend(" ! " + line for line in wrapped)
def append_data(lines, name, values):
    block_size = 100
    values_per_line = 3
    for start in range(0, len(values), block_size):
        stop = min(start + block_size, len(values))
        lines.append(
            f" data ({name}(n),n={start + 1},{stop}) / &"
        )
        for offset in range(start, stop, values_per_line):
            last = min(offset + values_per_line, stop)
            row = ", ".join(
                fortran_value(value) for value in values[offset:last]
            )
            suffix = " /" if last == stop else ", &"
            lines.append("      " + row + suffix)
        lines.append("")
def make_module(source):
    dimension = re.search(
        r"(?m)^\s*wavenumbers\s*=\s*(\d+)\s*;",
        source,
    )
    if dimension is None:
        raise ValueError("Missing wavenumbers dimension")
    count = int(dimension.group(1))
    if count != 2003:
        raise ValueError(
            f"Expected 2003 wavenumbers for this table, found {count}"
        )
    section = re.search(r"(?m)^\s*data\s*:", source)
    if section is None:
        raise ValueError("Missing CDL data: section")
    data_section = source[section.end():]
    array_names = (
        "wavenumbers",
        "self_absco_ref",
        "for_absco_ref",
        "for_closure_absco_ref",
        "self_texp",
    )
    arrays = {
        name: read_values(data_section, name, count)
        for name in array_names
    }
    ref_press = read_values(data_section, "ref_press", 1)[0]
    ref_temp = read_values(data_section, "ref_temp", 1)[0]
    for index, token in enumerate(arrays["wavenumbers"]):
        expected = Decimal(-20 + 10 * index)
        if decimal_value(token) != expected:
            raise ValueError(
                f"Unexpected wavenumber at index {index + 1}: {token}"
            )
    for name in (
        "self_absco_ref",
        "for_absco_ref",
        "for_closure_absco_ref",
    ):
        if any(decimal_value(value) < 0 for value in arrays[name]):
            raise ValueError(f"{name}: negative continuum coefficient")
    if decimal_value(ref_press) <= 0:
        raise ValueError("ref_press must be positive")
    if decimal_value(ref_temp) <= 0:
        raise ValueError("ref_temp must be positive")
    lines = [
        "module WATER_VAROP_CONT",
        "",
    ]
    for name in ("Title", "Copyright", "Reference", "Notes"):
        append_comment(lines, name + ":", read_attribute(source, name))
        lines.append("")
    lines.extend([
        " implicit none",
        " private",
        "",
        " integer :: n",
        f" integer, parameter, public :: nwvcont = {count}",
        "",
        " ! Reference pressure from the table, in mbar.",
        " real(8), parameter, public :: water_cont_ref_press = &",
        "      " + fortran_value(ref_press),
        "",
        " ! Reference temperature, in K.",
        " real(8), parameter, public :: water_cont_ref_temp = &",
        "      " + fortran_value(ref_temp),
        "",
        " ! Conversion of table reference pressure from mbar to Pa.",
        " real(8), parameter, public :: water_cont_ref_press_pa = &",
        "      100.d0 * water_cont_ref_press",
        "",
        " real(8), protected, public, save :: wavenumbers(nwvcont)",
        " real(8), protected, public, save :: self_absco_ref(nwvcont)",
        " real(8), protected, public, save :: for_absco_ref(nwvcont)",
        " real(8), protected, public, save :: for_closure_absco_ref(nwvcont)",
        " real(8), protected, public, save :: self_texp(nwvcont)",
        "",
        " ! Dimensionless layer optical depths on the observed grid.",
        " real(8), allocatable, public :: water_cont_layer(:,:)",
        "",
        " ! Dimensionless total optical depth on the observed grid.",
        " real(8), allocatable, public :: water_cont_tau(:)",
        "",
    ])
    for name in array_names:
        append_data(lines, name, arrays[name])
    lines.append("end module WATER_VAROP_CONT")
    lines.append("")
    return "\n".join(lines)
def main():
    parser = argparse.ArgumentParser(
        description="Convert the MT_CKD CDL text dump to a Fortran DATA module."
    )
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    try:
        source = args.input.read_text(encoding="utf-8-sig")
        result = make_module(source)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")
    sys.stdout.write(result)
if __name__ == "__main__":
    main()