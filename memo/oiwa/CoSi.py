"""
example of molecule and crystal models.
"""

models = []

# ==================================================
x_X = 0.143563986
x_Si = 0.843538046
# x_X = 0.25
# x_Si = 0.75
a = 4.4386000633000000

models = []
Nb_max = 0
for Nb in (
    3,
    50,
):
    models.append(
        {
            "model": f"B20_spinless_Co_spd_Si_sp_all_{Nb}_fixed",
            "group": 198,
            "cell": {"a": a, "b": a, "c": a},
            #
            "site": {
                "Co": (f"[{x_X},{x_X},{x_X}]", ["s", "p", "d"]),
                "Si": (f"[{x_Si},{x_Si},{x_Si}]", ["s", "p"]),
            },
            "bond": [
                ("Co", "Co", [i for i in range(1, Nb + 1)]),
                ("Si", "Si", [i for i in range(1, Nb + 1)]),
                ("Co", "Si", [i for i in range(1, Nb + 1)]),
            ],
            #
            "spinful": False,
            #
            "qtdraw": {"view": [0, 0, 1]},
            "SAMB_select": {
                # "X": ["T", "M"],  # type, "Q/G/T/M".
                "X": ["Q", "G"],  # type, "Q/G/T/M".
                "l": [],  # rank, []=all.
                # "Gamma": ["A", "T"],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
                "Gamma": [
                    "A"
                ],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
                "s": [],  # spin, []=all.
            },  # X, l, Gamma, s.
            "toroidal_priority": False,
            "max_neighbor": Nb,  # for DFT fitting, e.g. =200.
            "search_cell_range": (
                -10,
                10,
                -10,
                10,
                -10,
                10,
            ),  # (a1, a2, a3) range. for DFT fitting, e.g.=(-10, 10, -10, 10, -10, 10).
            "qtdraw": {
                "create": False,
                "mode": "standard",  # mode, "standard/detail".
                "background": True,  # create file or interactive.
                "view": None,  # [a,b,c], if None, default=[6,5,1].
                #
                "rep_site": True,  # show representative site ?
                "rep_bond": True,  # show representative bond ?
                #
                "max_neighbor": 3,  # max. neighbor to draw.
                "cell_mode": None,  # "off/single/all", if None, single for SG, off for PG.
                "scale": None,  # default or magnification scale.
                "site_radius": 0.07,  # base site radius.
                "bond_width": 0.02,  # base bond width.
            },
            "pdf": {
                "create": True,
                "common_samb": True,  # display common SAMB ?
                "harmonics": True,  # display harmonics ?
                "max_neighbor": Nb,  # max. neighbor to write (None: all).
            },
        }
    )


x_X = 0.856436014
x_Si = 0.156461954
# x_X = 0.25
# x_Si = 0.75
a = 4.4386000633000000

models = []
# Nb_max = 0
for Nb in (1, 3, 5, 50):
    models.append(
        {
            "model": f"B20_spinless_Co_spd_Si_sp_all_L_{Nb}_fixed",
            "group": 198,
            "cell": {"a": a, "b": a, "c": a},
            #
            "site": {
                "Co": (f"[{x_X},{x_X},{x_X}]", ["s", "p", "d"]),
                "Si": (f"[{x_Si},{x_Si},{x_Si}]", ["s", "p"]),
            },
            "bond": [
                ("Co", "Co", [i for i in range(1, Nb + 1)]),
                ("Si", "Si", [i for i in range(1, Nb + 1)]),
                ("Co", "Si", [i for i in range(1, Nb + 1)]),
            ],
            #
            "spinful": False,
            #
            "qtdraw": {"view": [0, 0, 1]},
            "SAMB_select": {
                "X": ["Q", "G"],  # type, "Q/G/T/M".
                "l": [],  # rank, []=all.
                "Gamma": [
                    "A"
                ],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
                "s": [],  # spin, []=all.
            },  # X, l, Gamma, s.
            "toroidal_priority": False,
            "max_neighbor": Nb,  # for DFT fitting, e.g. =200.
            "search_cell_range": (
                -10,
                10,
                -10,
                10,
                -10,
                10,
            ),  # (a1, a2, a3) range. for DFT fitting, e.g.=(-10, 10, -10, 10, -10, 10).
            "qtdraw": {
                "mode": "standard",  # mode, "standard/detail".
                "background": True,  # create file or interactive.
                "view": None,  # [a,b,c], if None, default=[6,5,1].
                #
                "rep_site": True,  # show representative site ?
                "rep_bond": True,  # show representative bond ?
                #
                "max_neighbor": 3,  # max. neighbor to draw.
                "cell_mode": None,  # "off/single/all", if None, single for SG, off for PG.
                "scale": None,  # default or magnification scale.
                "site_radius": 0.07,  # base site radius.
                "bond_width": 0.02,  # base bond width.
            },
            "pdf": {
                "common_samb": True,  # display common SAMB ?
                "harmonics": True,  # display harmonics ?
                "max_neighbor": Nb,  # max. neighbor to write (None: all).
            },
        }
    )
