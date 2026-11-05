"""
Create harmonics mapping from D6h and Oh to each PG.
    - harmonics_mapping: mapping and coeffcient to each PG.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import os
from multipiegen.data.data_top_harmonics import top_harmonics_tag
from multipiegen.data.data_compatibility_relation import compatibility_relation

# ==================================================
# polar (Q) and axial (G).
HARMONICS_TYPES = ("Q", "G")

# ==================================================
header = """
Mapping of harmonics key to each PG.

- data_harmonics_mapping.
    { PG: { "top": "Oh/D6h", "Q/G": { pg_harm_key: (coeff, top_harm_key) } } }.
    pg_harm_key = (rank, irrep, multiplicity (1-based), component (0-based)) [-1: no multiplicity/single component]
    top_harm_key = (rank, irrep, multiplicity, component)
"""


# ==================================================
# mapping (top group harmonics -> subgroup harmonics).
# ==================================================
def _top_irrep_comp(irrep, comp, sub_irrep):
    """
    Resolve (irrep, comp) of the top harmonics to be referred.

    If the irrep. of the subgroup is a complex one (contains "a"/"b"), the complex expression of the top group
    (irrep+"a" / irrep+"b", single component) is used instead of the real one.

    Args:
        irrep (str): irrep. of top group.
        comp (int): component of top group.
        sub_irrep (str): irrep. of subgroup (polar or axial).

    Returns:
        - (str, int) -- irrep. and component of top harmonics to be referred.
    """
    if "a" in sub_irrep or "b" in sub_irrep:
        if comp == 0:
            return irrep + "a", -1
        if comp == 1:
            return irrep + "b", -1
    return irrep, comp


# ==================================================
def create_mapping(group_tag):
    """
    Create mapping from top group harmonics to subgroup harmonics (no sympy data are needed).

    Args:
        group_tag (str): point group tag (with "-c" for complex expression).

    Returns:
        - (dict) -- mapping, {"Q": {...}, "G": {...}}.

    Notes:
        - mapping[X] = {(l, irrep_sub, comp_sub): [(coeff, (top, l, irrep_top, mul_top, comp_top)), ...]}
          ({(int,str,int): [(int, (str,int,str,int,int))]}).
        - the key is the harmonics of the subgroup, and the list is top group harmonics belonging to it.
        - the list order follows top_harmonics_tag.
        - multiplicity of the subgroup is not yet assigned (see renumber_multiplicity).
    """
    top, cr = compatibility_relation[group_tag]

    mapping = {X: {} for X in HARMONICS_TYPES}
    mapping["top"] = top
    for rank, irrep, mul, comp in top_harmonics_tag[top]:
        polar, axial = cr[(irrep, comp)]
        for X, (sub_irrep, sub_comp, coeff) in zip(HARMONICS_TYPES, (polar, axial)):
            top_irrep, top_comp = _top_irrep_comp(irrep, comp, sub_irrep)
            key = (rank, sub_irrep, sub_comp)
            mapping[X].setdefault(key, []).append((coeff, (rank, top_irrep, mul, top_comp)))

    return mapping


# ==================================================
def create_mappings(group_tags=None):
    """
    Create mappings for several point groups.

    Args:
        group_tags (list, optional): point group tags. If None, all groups in compatibility_relation.

    Returns:
        - (dict) -- {group_tag: mapping}.
    """
    if group_tags is None:
        group_tags = list(compatibility_relation)
    return {tag: create_mapping(tag) for tag in group_tags}


# ==================================================
def renumber_multiplicity(mapping):
    """
    Assign multiplicity of the subgroup harmonics.

    Args:
        mapping (dict): mapping created by create_mapping.

    Returns:
        - (dict) -- {"Q": {...}, "G": {...}}, each is {(l, irrep, n, comp): (coeff, top_key)}.

    Notes:
        - n = -1 if only one top harmonics belongs to the key, otherwise n = 1, 2, ... in the list order.
    """
    renumbered = {X: {} for X in HARMONICS_TYPES}
    renumbered["top"] = mapping["top"]
    del mapping["top"]
    for X, mapping_X in mapping.items():
        for (rank, irrep, comp), data in mapping_X.items():
            if len(data) == 1:
                renumbered[X][(rank, irrep, -1, comp)] = data[0]
            else:
                for n, v in enumerate(data, start=1):
                    renumbered[X][(rank, irrep, n, comp)] = v
    return renumbered


# ==================================================
def create_harmonics_mapping():
    group_tags = list(compatibility_relation)
    mappings = create_mappings(group_tags)
    mappings = {tag: renumber_multiplicity(m) for tag, m in mappings.items()}

    return mappings


# ==================================================
def create_file():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "data_harmonics_mapping.py")
    mapping = create_harmonics_mapping()
    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("data_harmonics_mapping =", mapping, file=f)


# ==================================================
create_file()
