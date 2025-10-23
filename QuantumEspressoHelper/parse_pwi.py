#!/usr/bin/env python3
"""
qe_in_parser.py — Parse Quantum ESPRESSO pw.x input files to JSON.

Outputs a dict with:
{
  "num_atoms": int,
  "cell_data": [ax,ay,az, bx,by,bz, cx,cy,cz],   # Å
  "info": { "source": "qe_input", ... },
  "atoms_data": [
    {"symbols": "Mo", "positions": [x,y,z]},     # Å
    ...
  ]
}
"""
from __future__ import annotations
import sys, re, json, argparse
from pathlib import Path
from glob import glob
from typing import Iterator

# ---- constants (same as your out-parser) ----
kBohrToAngstrom = 0.529177

# -------------------- utils --------------------
_CAPS_LINE = re.compile(r"^[A-Z][A-Z0-9_]*(\s|$)")
_WS_ONLY = re.compile(r"^\s*$")

def _strip_comment(line: str) -> str:
    # QE treats '!' as comment; keep content before the first unquoted '!'
    # (inputs are simple—ignore quoted edge cases for brevity)
    idx = line.find("!")
    return line[:idx] if idx >= 0 else line

def _clean_lines(text: str) -> list[str]:
    out = []
    for raw in text.splitlines():
        line = _strip_comment(raw).rstrip()
        if line:
            out.append(line)
        else:
            out.append("")  # keep block separators
    return out

def _find_block_start(lines: list[str], key: str) -> int | None:
    for i, ln in enumerate(lines):
        if ln.upper().startswith(key):
            return i
    return None

# --- Replace the old card detection with this ---
_QE_CARD_NAMES = {
    # sections / cards commonly appearing in pw.x inputs
    "K_POINTS", "CELL_PARAMETERS", "ATOMIC_POSITIONS", "ATOMIC_SPECIES",
    "OCCUPATIONS", "CONSTRAINTS", "HUBBARD", "ATOMIC_FORCES",
    "STARTING_NS_EIGENVECTORS", "ATOMIC_VELOCITIES",
    # also common cards in tools/workflows
    "END", "EOF"
}

def _is_card_header(line: str) -> bool:
    s = line.strip()
    if not s:
        return False
    # Namelists like &control, &system, &electrons, &ions, &cell, &phonon
    if s.startswith("&"):
        return True
    # Known cards, possibly with unit suffixes (e.g., "ATOMIC_POSITIONS angstrom")
    first = s.split()[0]
    return first.upper() in _QE_CARD_NAMES

def _block_iter(lines: list[str], start_idx: int):
    """Yield lines after 'start_idx' until a blank line or the next QE card header."""
    i = start_idx + 1
    while i < len(lines):
        raw = lines[i]
        if not raw.strip():
            break
        if _is_card_header(raw):
            break
        yield raw
        i += 1

#def _block_iter(lines: list[str], start_idx: int) -> Iterator[str]:
#    """Yield lines after 'start_idx' until blank or next ALLCAPS card."""
#    i = start_idx + 1
#    while i < len(lines):
#        ln = lines[i].strip()
#        if not ln or _CAPS_LINE.match(lines[i]):  # blank or next card
#            break
#        yield lines[i]
#        i += 1

def _parse_named_value(lines: list[str], name_regex: str) -> float | None:
    """Find scalar like celldm(1) = 6.0  or  A = 5.43 in any &namelist."""
    pat = re.compile(rf"\b{name_regex}\b\s*=\s*([0-9.+\-Ee]+)")
    for ln in lines:
        m = pat.search(ln)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
    return None

def _parse_units(header_line: str, default_unit: str) -> str:
    # accept: "CELL_PARAMETERS angstrom" or "ATOMIC_POSITIONS (alat)"
    m = re.search(r"\b(CELL_PARAMETERS|ATOMIC_POSITIONS)\b\s*[\(\s]*([A-Za-z_]+)[\)]?\s*$",
                  header_line, flags=re.IGNORECASE)
    if m and m.group(2):
        return m.group(2).lower()
    return default_unit.lower()

# ---------------- lattice helpers ----------------
def _parse_cell(lines: list[str]) -> tuple[list[float], str] | None:
    """
    Return ([ax,ay,az, bx,by,bz, cx,cy,cz] in Å, unit_tag_used)
    Supported units: angstrom, bohr, alat
    """
    i = _find_block_start(lines, "CELL_PARAMETERS")
    if i is None:
        return None
    unit = _parse_units(lines[i], "alat")  # QE defaults CELL_PARAMETERS to 'alat'
    rows = []
    for j, ln in enumerate(_block_iter(lines, i)):
        parts = ln.split()
        if len(parts) < 3:
            break
        rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
        if len(rows) == 3:
            break
    if len(rows) != 3:
        raise ValueError("CELL_PARAMETERS: need 3 rows with 3 numbers each")

    # Convert to Å
    if unit == "angstrom" or unit == "ang":
        pass  # already Å
    elif unit == "bohr" or unit == "a.u.":
        rows = [[v * kBohrToAngstrom for v in r] for r in rows]
    elif unit == "alat":
        # Need alat in Å: derive from celldm(1) (Bohr) or A (Å)
        alat_bohr = _parse_named_value(lines, r"celldm\s*\(\s*1\s*\)")
        alat_ang = _parse_named_value(lines, r"A")
        if alat_ang is not None:
            scale = alat_ang
        elif alat_bohr is not None:
            scale = alat_bohr * kBohrToAngstrom
        else:
            raise ValueError("CELL_PARAMETERS (alat): missing celldm(1) or A to resolve alat")
        rows = [[v * scale for v in r] for r in rows]
    else:
        raise ValueError(f"CELL_PARAMETERS: unsupported unit '{unit}'")

    flat = [rows[0][0], rows[0][1], rows[0][2],
            rows[1][0], rows[1][1], rows[1][2],
            rows[2][0], rows[2][1], rows[2][2]]
    return flat, unit

def _matmul3(M: list[list[float]], v: list[float]) -> list[float]:
    return [
        M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
        M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
        M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2],
    ]

# ---------------- positions ----------------
def _parse_positions(lines: list[str], lattice_A: list[float] | None) -> tuple[list[dict], str] | None:
    """
    Return (atoms_data, unit_tag_used). Positions converted to Å.
    Supports units: angstrom, bohr, alat, crystal
    """
    i = _find_block_start(lines, "ATOMIC_POSITIONS")
    if i is None:
        return None
    unit = _parse_units(lines[i], "alat")  # QE default is alat
    # Prepare lattice matrix if needed
    L = None
    if lattice_A is not None:
        L = [
            [lattice_A[0], lattice_A[1], lattice_A[2]],
            [lattice_A[3], lattice_A[4], lattice_A[5]],
            [lattice_A[6], lattice_A[7], lattice_A[8]],
        ]

    # alat scaling in Å
    alat_bohr = _parse_named_value(lines, r"celldm\s*\(\s*1\s*\)")
    alat_ang = _parse_named_value(lines, r"A")
    alat_scale_A: float | None = None
    if alat_ang is not None:
        alat_scale_A = alat_ang
    elif alat_bohr is not None:
        alat_scale_A = alat_bohr * kBohrToAngstrom

    atoms: list[dict] = []
    for ln in _block_iter(lines, i):
        parts = ln.split()
        if len(parts) < 4:
            break
        sym = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            break  # stop at malformed line

        # Ignore any trailing if_pos flags or constraints
        posA: list[float]
        if unit in ("angstrom", "ang"):
            posA = [x, y, z]
        elif unit in ("bohr", "a.u."):
            posA = [x * kBohrToAngstrom, y * kBohrToAngstrom, z * kBohrToAngstrom]
        elif unit == "alat":
            if alat_scale_A is None:
                raise ValueError("ATOMIC_POSITIONS (alat): missing celldm(1) or A to resolve alat")
            posA = [x * alat_scale_A, y * alat_scale_A, z * alat_scale_A]
        elif unit == "crystal":
            if L is None:
                raise ValueError("ATOMIC_POSITIONS (crystal) requires CELL_PARAMETERS")
            posA = _matmul3(L, [x, y, z])
        else:
            raise ValueError(f"ATOMIC_POSITIONS: unsupported unit '{unit}'")

        atoms.append({"symbols": sym, "positions": posA})

    if not atoms:
        return None
    return atoms, unit

# ---------------- species (optional) ----------------
def _parse_species(lines: list[str]) -> list[tuple[str, float, str]] | None:
    """Return [(symbol, mass, pseudo)], if present. Not required for output shape."""
    i = _find_block_start(lines, "ATOMIC_SPECIES")
    if i is None:
        return None
    out = []
    for ln in _block_iter(lines, i):
        parts = ln.split()
        if len(parts) < 3:
            break
        sym = parts[0]
        try:
            mass = float(parts[1])
        except ValueError:
            break
        pseudo = parts[2]
        out.append((sym, mass, pseudo))
    return out or None

# ---------------- top-level parse ----------------
def qe_in_to_json(file_path: str | Path) -> dict:
    """
    Parse QE pw.x input and return json-compatible dict matching the output-parser shape.
    """
    text = Path(file_path).read_text(encoding="utf-8", errors="replace")
    lines = _clean_lines(text)

    # CELL_PARAMETERS (→ Å)
    cell_tuple = _parse_cell(lines)
    if cell_tuple is None:
        raise ValueError("Could not find CELL_PARAMETERS block")
    cell_A, cell_unit = cell_tuple

    # ATOMIC_POSITIONS (→ Å)
    pos_tuple = _parse_positions(lines, lattice_A=cell_A)
    if pos_tuple is None:
        raise ValueError("Could not find ATOMIC_POSITIONS block")
    atoms, pos_unit = pos_tuple

    # Shape compatible with your out-parser
    data: dict = {
        "num_atoms": len(atoms),
        "cell_data": cell_A,  # Å
        "info": {
            "source": "qe_input",
            "cell_unit": cell_unit,
            "positions_unit": pos_unit,
        },
        "atoms_data": atoms,
    }

    # Optional: species metadata (doesn't alter shape required for coords/box/symbols)
    species = _parse_species(lines)
    if species:
        data["info"]["species"] = [{"symbol": s, "mass": m, "pseudo": p} for (s, m, p) in species]

    return data

# ---------------- CLI (same style as your out-parser) ----------------
def _resolve_inputs(globs_or_paths: list[str]) -> list[Path]:
    files: list[Path] = []
    for pat in globs_or_paths:
        matches = [Path(p) for p in (glob(pat) if any(c in pat for c in "*?[]") else [pat])]
        files.extend(m for m in matches if m.is_file())
    # dedup
    seen, out = set(), []
    for p in files:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            out.append(p)
    return out

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Parse Quantum ESPRESSO pw.x input to JSON (Å).")
    ap.add_argument("input", nargs="+", help="QE input file(s) or globs (e.g., '*.in')")
    ap.add_argument("-o", "--out", default=None,
                    help="Output JSON file (single input) or output directory (multiple). "
                         "If omitted, prints to stdout.")
    ap.add_argument("--pretty", action="store_true", help="Pretty-print JSON.")
    ap.add_argument("--ndjson", action="store_true",
                    help="Stream per-file JSON lines to stdout (for multiple inputs).")
    args = ap.parse_args(argv)

    inputs = _resolve_inputs(args.input)
    if not inputs:
        print("No input files found.", file=sys.stderr)
        return 2

    indent = 2 if args.pretty else None

    if len(inputs) == 1 and not args.ndjson:
        data = qe_in_to_json(inputs[0])
        if args.out:
            outp = Path(args.out)
            if outp.exists() and outp.is_dir():
                outp = outp / (inputs[0].stem + ".json")
            outp.write_text(json.dumps(data, ensure_ascii=False, indent=indent) + ("\n" if indent else ""),
                            encoding="utf-8")
        else:
            json.dump(data, sys.stdout, ensure_ascii=False, indent=indent)
            if indent:
                sys.stdout.write("\n")
        return 0


    # multiple inputs
    if args.ndjson:
        for p in inputs:
            try:
                obj = qe_in_to_json(p)
                sys.stdout.write(json.dumps({"file": str(p), "data": obj}, ensure_ascii=False) + "\n")
            except Exception as e:
                sys.stdout.write(json.dumps({"file": str(p), "error": str(e)}, ensure_ascii=False) + "\n")
        return 0

    if args.out:
        out_dir = Path(args.out)
        out_dir.mkdir(parents=True, exist_ok=True)
        for p in inputs:
            try:
                obj = qe_in_to_json(p)
            except Exception as e:
                print(f"Skipping {p}: {e}", file=sys.stderr)
                continue
            (out_dir / (p.stem + ".json")).write_text(
                json.dumps(obj, ensure_ascii=False, indent=indent) + ("\n" if indent else ""),
                encoding="utf-8"
            )
        return 0
    else:
        # combined list to stdout
        all_objs = []
        for p in inputs:
            try:
                all_objs.append({"file": str(p), "data": qe_in_to_json(p)})
            except Exception as e:
                all_objs.append({"file": str(p), "error": str(e)})
        json.dump(all_objs, sys.stdout, ensure_ascii=False, indent=indent)
        if indent:
            sys.stdout.write("\n")
        return 0

if __name__ == "__main__":
    raise SystemExit(main())

