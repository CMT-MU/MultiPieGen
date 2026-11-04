import numpy as np
from multipie import Group
from multipiegen.data.data_parent_pg import data_representative_parent_pg, data_parent_pg


def check_each_pg(c_group):
    p_group = data_representative_parent_pg[c_group]
    print(f"{p_group} => {c_group}")
    cg = Group(c_group)
    pg = Group(p_group)
    mapping = np.array(data_parent_pg[c_group][p_group][0][0]) - 1
    pg_wp = list(pg.wyckoff["site"].keys())[-1]
    cg_wp = list(cg.wyckoff["site"].keys())[-1]

    pw = pg.wyckoff["site"][pg_wp]["conventional"][mapping]
    cw = cg.wyckoff["site"][cg_wp]["conventional"]
    site_err = not bool(np.all(pw == cw))
    if site_err:
        print(c_group, pw, cw)

    pw = pg.wyckoff["bond"][pg.wyckoff["site"][pg_wp]["bond"][0]]["conventional"][mapping]
    cw = cg.wyckoff["bond"][cg.wyckoff["site"][cg_wp]["bond"][0]]["conventional"]
    bond_err = not bool(np.all(pw == cw))
    if bond_err:
        print(c_group, pw, cw)

    return site_err, bond_err


def check_pg():
    cnt_site = 0
    cnt_bond = 0
    for c_group in data_representative_parent_pg.keys():
        site_err, bond_err = check_each_pg(c_group)
        if site_err:
            cnt_site += 1
        if bond_err:
            bond_err += 1
    if cnt_site > 0:
        print("found error in site wyckoff.")
    if cnt_bond > 0:
        print("found error in bond wyckoff.")
    if cnt_site == 0 and cnt_bond == 0:
        print("no error.")


check_pg()
