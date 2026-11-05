# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Note:
definition of gv is changed. (it is multiplied by (-1)).
change name.
"""

"""
Misc. for response tensor.
"""

# list of multipoles in response tensor.
rt_multipole_list = {
    "cubic": [
        "s",
        "px",
        "py",
        "pz",
        "du",
        "dv",
        "dyz",
        "dxz",
        "dxy",
        "f3",
        "fax",
        "fay",
        "faz",
        "fbx",
        "fby",
        "fbz",
        "g",
        "gu",
        "gv",
        "gax",
        "gay",
        "gaz",
        "gbx",
        "gby",
        "gbz",
    ],
    "hexagonal": [
        "s",
        "pz",
        "du",
        "fs",
        "fc",
        "faz",
        "g0",
        "gc",
        "gs",
        "px",
        "py",
        "dyz",
        "dxz",
        "dv",
        "dxy",
        "fx",
        "fy",
        "f3",
        "fbz",
        "gau",
        "gav",
        "gc4",
        "gaz",
        "gv",
        "gbz",
    ],
}

# cartesian component in terms of multipoles.
rt_cartesian_multipole = {
    "cubic": {
        "s": [(1, "s")],
        "x": [(1, "px")],
        "y": [(1, "py")],
        "z": [(1, "pz")],
        "xx": [(-1, "du"), (1, "dv")],
        "xy": [(1, "dxy")],
        "xz": [(1, "dxz")],
        "yy": [(-1, "du"), (-1, "dv")],
        "yz": [(1, "dyz")],
        "zz": [(2, "du")],
        "xxx": [(2, "fax")],
        "xxy": [(-1, "fay"), (-1, "fby")],
        "xxz": [(-1, "faz"), (1, "fbz")],
        "xyy": [(-1, "fax"), (1, "fbx")],
        "xyz": [(1, "f3")],
        "xzz": [(-1, "fax"), (-1, "fbx")],
        "yyy": [(2, "fay")],
        "yyz": [(-1, "faz"), (-1, "fbz")],
        "yzz": [(-1, "fay"), (1, "fby")],
        "zzz": [(2, "faz")],
        "xxxx": [(2, "g"), (-1, "gu"), (-1, "gv")],
        "xxxy": [(1, "gaz"), (-1, "gbz")],
        "xxxz": [(-1, "gay"), (-1, "gby")],
        "xxyy": [(-1, "g"), (2, "gu")],
        "xxyz": [(2, "gbx")],
        "xxzz": [(-1, "g"), (-1, "gu"), (1, "gv")],
        "xyyy": [(-1, "gaz"), (-1, "gbz")],
        "xyyz": [(2, "gby")],
        "xyzz": [(2, "gbz")],
        "xzzz": [(1, "gay"), (-1, "gby")],
        "yyyy": [(2, "g"), (-1, "gu"), (1, "gv")],
        "yyyz": [(1, "gax"), (-1, "gbx")],
        "yyzz": [(-1, "g"), (-1, "gu"), (-1, "gv")],
        "yzzz": [(-1, "gax"), (-1, "gbx")],
        "zzzz": [(2, "g"), (2, "gu")],
    },
    "hexagonal": {
        "s": [(1, "s")],  # same as cubic.
        "x": [(1, "px")],  # same as cubic.
        "y": [(1, "py")],  # same as cubic.
        "z": [(1, "pz")],  # same as cubic.
        "xx": [(-1, "du"), (1, "dv")],  # same as cubic.
        "xy": [(1, "dxy")],  # same as cubic.
        "xz": [(1, "dxz")],  # same as cubic.
        "yy": [(-1, "du"), (-1, "dv")],  # same as cubic.
        "yz": [(1, "dyz")],  # same as cubic.
        "zz": [(2, "du")],  # same as cubic.
        "xxx": [(-3, "fx"), (1, "fc")],
        "xxy": [(-1, "fy"), (1, "fs")],
        "xxz": [(-1, "faz"), (1, "fbz")],  # same as cubic.
        "xyy": [(-1, "fx"), (-1, "fc")],
        "xyz": [(1, "f3")],  # same as cubic.
        "xzz": [(4, "fx")],
        "yyy": [(-3, "fy"), (-1, "fs")],
        "yyz": [(-1, "faz"), (-1, "fbz")],  # same as cubic.
        "yzz": [(4, "fy")],
        "zzz": [(2, "faz")],  # same as cubic.
        "xxxx": [(3, "g0"), (1, "gc4"), (-1, "gv")],  # fix incorrect expression.
        "xxxy": [(1, "gaz"), (-1, "gbz")],  # same as cubic.
        "xxxz": [(1, "gc"), (-3, "gav")],
        "xxyy": [(1, "g0"), (-1, "gc4")],
        "xxyz": [(-1, "gau"), (1, "gs")],
        "xxzz": [(-4, "g0"), (1, "gv")],  # fix incorrect expression.
        "xyyy": [(-1, "gaz"), (-1, "gbz")],  # same as cubic.
        "xyyz": [(-1, "gc"), (-1, "gav")],
        "xyzz": [(2, "gbz")],  # same as cubic.
        "xzzz": [(4, "gav")],
        "yyyy": [(3, "g0"), (1, "gc4"), (1, "gv")],  # fix incorrect expression.
        "yyyz": [(-3, "gau"), (-1, "gs")],
        "yyzz": [(-4, "g0"), (-1, "gv")],  # fix incorrect expression.
        "yzzz": [(4, "gau")],
        "zzzz": [(8, "g0")],
    },
}
