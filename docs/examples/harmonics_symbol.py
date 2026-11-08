# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: polar and axial symbols of harmonics of point groups, and the harmonics PDF.

    python docs/examples/harmonics_symbol.py

An axial harmonics has the same expression as the polar one; only the irrep, component,
and sign change. The PDF (others/pdf/PG/{no}-{tag}/harmonics.pdf) is created by
multipiegen/create_pdf/create_pdf_harmonics.py.
"""

from multipiegen.util.util_harmonics_symbol import harmonics_symbol_table, polar_to_axial, multipole_symbol

# ==================================================
# 1. table of Oh harmonics: (key, polar symbol, axial symbol, expression in LaTeX).
tbl = harmonics_symbol_table("Oh")
print("number of Oh harmonics (rank 0-11):", len(tbl))
for key, p_sym, a_sym, expr in tbl[:9]:
    print(f"    {str(key):22s} {p_sym:28s} {a_sym:30s} {expr}")

# ==================================================
# 2. for Oh and D6h, the axial irrep is obtained by g <-> u (same component, no sign).
for key, p_sym, a_sym, expr in harmonics_symbol_table("D6h"):
    if key[0] == 3:
        print(f"    D6h: {p_sym:28s} -> {a_sym}")

# ==================================================
# 3. polar -> axial for a single key (any of the 47 point groups).
a_key, sign = polar_to_axial("Oh", (2, "Eg", -1, 1))
print("\nOh  Q(Eg,2) ->", a_key, sign, multipole_symbol("G", a_key, sign))

#   Td: components are interchanged and a sign appears.
for key, p_sym, a_sym, expr in harmonics_symbol_table("Td")[4:9]:
    print(f"    Td: {p_sym:24s} -> {a_sym:26s} {expr}")

#   complex group: complex irreps E^(a), E^(b) and complex expressions.
for key, p_sym, a_sym, expr in harmonics_symbol_table("C3h-c")[1:4]:
    print(f"    C3h-c: {p_sym:30s} -> {a_sym:30s} {expr}")

# ==================================================
# 4. create the PDFs of all 47 point groups (needs ptex2pdf):
#    python -m multipiegen.create_pdf.create_pdf_harmonics
