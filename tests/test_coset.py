# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pickle

import numpy as np
import pytest
from multipiegen.create_binary.create_coset import coset_data
from multipiegen.util import util_coset
from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_cluster_samb_image import cluster_samb_transport, orbit_correspondence, _mat4
from multipiegen.util.util_coset import load_coset
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

info = load_info()
group = load_group()
data = load_coset()


def _transport(gid):
    f = data[gid]["first"]
    if f["kind"] != "inversion":
        return None
    pair = (gid, f["inversion"]["image"], _mat4(str(f["inversion"]["map"].tolist())), gid.split(":")[0])
    return cluster_samb_transport(load_cluster_samb(), orbit_correspondence(group, [pair]))[gid]


# ==================================================
def test_groups():
    # real PG and SG with a parent (no "-c", no top group).
    ids = info["id_set"]["PG"]["all"] + info["id_set"]["SG"]["all"]
    assert list(data.keys()) == [g for g in ids if util_coset.parent_chain(g)]
    assert len([g for g in data if g.startswith("PG")]) == 35 and len([g for g in data if g.startswith("SG")]) == 216
    kinds = [d["first"]["kind"] for d in data.values()]
    assert kinds.count("inversion") == 77  # 12 PG, 65 SG, all chiral.
    for gid, d in data.items():
        if d["first"]["kind"] == "inversion":
            assert all(
                np.linalg.det(np.array(W, dtype=float)[:3, :3]) > 0 for W in group[gid]["symmetry_operation"]["fractional"]
            )


# ==================================================
@pytest.mark.parametrize("gid", ["PG:8", "PG:15", "PG:1", "SG:152", "SG:198", "SG:225"])
def test_binary_is_up_to_date(gid):
    new = pickle.loads(pickle.dumps(coset_data(gid, _transport(gid) if gid in data else None)))
    if gid not in data:
        assert new is None
        return
    assert repr(new) == repr(data[gid])


# ==================================================
def test_first_coset():
    f = data["SG:152"]["first"]
    assert (f["parent"], f["rep"], f["order"], f["kind"], f["index"]) == ("SG:164", "-1:[0,0,0]", 2, "inversion", 2)
    assert f["inversion"]["image"] == "SG:154" and f["inversion"]["wyckoff"]["3a"] == ("3a", {"x": "-x"})
    assert f["parent_limit"]["3a"][0]["param"] == {"x": "0"}
    f = data[info["id"]["D2h"]]["first"]
    assert (f["parent"], f["rep"], f["order"], f["kind"]) == (info["id"]["D4h"], "4+[001]", 2, "rotation")
    assert f["frame"] == f["rep_matrix"] and f["inversion"] is None
    f = data[info["id"]["D4h"]]["first"]
    assert (f["rep"], f["order"]) == ("3+[111]", 3)
    assert data["SG:152"]["parent"] == ["SG:164", "SG:191"]
    assert data["SG:152"]["coset"]["SG:164"][0][0] == "-1:[0,0,0]"


# ==================================================
def test_cluster_samb_image():
    # 1847 Wyckoff positions of the 77 groups: label by label (eps = +-1) or U (unitary, within an irrep).
    n_diag = n_u = 0
    for gid, d in data.items():
        inv = d["first"]["inversion"]
        if inv is None:
            continue
        for kind in ["site", "bond"]:
            assert list(inv[kind]) == list(group[gid]["wyckoff"][kind])
            for w, t in inv[kind].items():
                assert sorted(t["perm"]) == list(range(len(group[gid]["wyckoff"][kind][w]["mapping"])))
                if t["diag"]:
                    n_diag += 1
                    assert all(complex(e) in (1, -1) for es in t["eps"].values() for e in es)
                else:
                    n_u += 1
                    assert t["image"] == (gid, w)  # only for a Wyckoff position mapped onto itself.
                    for rows, rows2, U in t["U"].values():
                        U = np.array(U, dtype=complex)
                        assert np.allclose(U @ U.conj().T, np.eye(len(rows)), atol=1e-8)
                        assert [k.Gamma for k, _ in rows] == [k.Gamma for k, _ in rows2]
    assert (n_diag, n_u) == (1809, 38)


# ==================================================
def test_inverse():
    # the map of an enantiomorphic pair and its reverse are inverse to each other (U^dagger, eps).
    for gid, d in data.items():
        inv = d["first"]["inversion"]
        if inv is None or inv["image"] == gid:
            continue
        rev = data[inv["image"]]["first"]["inversion"]
        assert rev["image"] == gid
        for kind in ["site", "bond"]:
            for w, t in inv[kind].items():
                r = rev[kind][t["image"][1]]
                assert r["image"] == (gid, w)
                assert [t["perm"][i] for i in r["perm"]] == list(range(len(t["perm"])))
                if t["diag"]:
                    assert r["diag"] and r["eps"] == t["eps"]
