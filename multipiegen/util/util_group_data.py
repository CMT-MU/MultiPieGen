# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Group ids, tags, id sets, and names, derived from the group information (data_group_pg/sg/mpg/msg).

- pg_tag, sg_tag, mpg_tag, msg_tag: { id: tag } (info.tag).
- pg_id, sg_id, mpg_id, msg_id: { tag: id }.
- pg_id_set: { "all": [real PG], "crystal": {crystal: [PG]}, "complex": [complex PG ("-c")], "irrep": {"1d/2d/3d":
  [PG with the largest irrep of this dimension]} }.
- sg_id_set: { "all", "crystal", "PG" }; mpg_id_set: { "all", "crystal", "PG", "type" }; msg_id_set: { "all",
  "crystal", "PG", "SG", "MPG", "type" }. Each grouping is in the order of the group information.
- group_name_list: { "PG/SG/MPG/MSG": { id: (dir_name, md_name, no, tag, name, PG_id, SG_id, MPG_id, MSG_id) } }.
"""

from multipie.util.util import str_to_sympy

from multipiegen.data.data_character_table import character
from multipiegen.data.data_group_mpg import mpg_info
from multipiegen.data.data_group_msg import msg_info
from multipiegen.data.data_group_pg import pg_info
from multipiegen.data.data_group_sg import sg_info


# ==================================================
def _group_by(info, field):
    d = {}
    for k, v in info.items():
        d.setdefault(getattr(v, field), []).append(k)
    return d


def _max_dim(tag):
    return max(int(str_to_sympy(c)[0]) for c in character[tag][1].values())


# ==================================================
pg_tag = {k: v.tag for k, v in pg_info.items()}
sg_tag = {k: v.tag for k, v in sg_info.items()}
mpg_tag = {k: v.tag for k, v in mpg_info.items()}
msg_tag = {k: v.tag for k, v in msg_info.items()}

pg_id = {v: k for k, v in pg_tag.items()}
sg_id = {v: k for k, v in sg_tag.items()}
mpg_id = {v: k for k, v in mpg_tag.items()}
msg_id = {v: k for k, v in msg_tag.items()}

_irrep = {}
for _k, _v in pg_info.items():
    _irrep.setdefault(f"{_max_dim(_v.tag)}d", []).append(_k)
pg_id_set = {
    "all": [k for k, v in pg_info.items() if not v.tag.endswith("-c")],
    "crystal": _group_by(pg_info, "crystal"),
    "complex": [k for k, v in pg_info.items() if v.tag.endswith("-c")],
    "irrep": _irrep,
}
sg_id_set = {"all": list(sg_info), "crystal": _group_by(sg_info, "crystal"), "PG": _group_by(sg_info, "PG")}
mpg_id_set = {
    "all": list(mpg_info),
    "crystal": _group_by(mpg_info, "crystal"),
    "PG": _group_by(mpg_info, "PG"),
    "type": _group_by(mpg_info, "type"),
}
msg_id_set = {
    "all": list(msg_info),
    "crystal": _group_by(msg_info, "crystal"),
    "PG": _group_by(msg_info, "PG"),
    "SG": _group_by(msg_info, "SG"),
    "MPG": _group_by(msg_info, "MPG"),
    "type": _group_by(msg_info, "type"),
}


# ==================================================
def _name_pg(k, v):
    no = int(k.split(":")[1])
    return (f"{no:03d}-{v.tag}", f"No_{no}.md", no, v.tag, f"${v.schoenflies}$ (${v.international}$)", k, v.SG, v.MPG, v.MSG)


def _name_sg(k, v):
    no, t = int(k.split(":")[1]), v.tag.replace("^", "_")
    return (f"{no:03d}-{t}", f"No_{no}.md", no, t, f"${v.schoenflies}$ (${v.international}$)", v.PG, k, v.MPG, v.MSG)


def _name_mpg(k, v):
    return (v.tag, f"No_{v.tag}.md", v.tag, v.tag, f"${v.international}$", v.PG, v.SG, k, v.MSG)


def _name_msg(k, v):
    return (v.tag, f"No_{v.tag}.md", v.tag, v.tag, f"${v.BNS}$", v.PG, v.SG, v.MPG, k)


group_name_list = {
    "PG": {k: _name_pg(k, v) for k, v in pg_info.items()},
    "SG": {k: _name_sg(k, v) for k, v in sg_info.items()},
    "MPG": {k: _name_mpg(k, v) for k, v in mpg_info.items()},
    "MSG": {k: _name_msg(k, v) for k, v in msg_info.items()},
}
