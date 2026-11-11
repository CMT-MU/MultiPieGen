# Data tables

The data tables in `multipiegen/data` are the source of all generated data. They are curated by hand; some of them are first created by a script in `others/creation` (it writes next to itself, and the result is copied into `multipiegen/data`). "origin" is "MultiPie Ver.2.4.0" for the tables taken from MultiPie Ver.2.4.0 (some of them modified), and "new" otherwise. The tables taken from MultiPie Ver.2.4.0 were copied once and have been maintained here since then. The generation never reads the data of MultiPie Ver.2; it uses only its general functions and utilities.

Data that can be derived from these tables are computed in `multipiegen` instead of being stored:

| data | computed in | from |
|---|---|---|
| `*_tag`, `*_id`, `*_id_set`, `group_name_list` | `util/util_group_data` | `*_info` of `data_group_pg/sg/mpg/msg` |
| `top_harmonics_symmetry`, `top_harmonics_u_matrix` | `util/util_top_harmonics` | `data_top_harmonics` (tesseral part from the multipolar harmonics) |
| `data_harmonics_mapping` | `util/util_harmonics_mapping` | `data_compatibility_relation`, `data_top_harmonics.top_harmonics_tag` |
| `op_mapping_d6h_oh_cluster`, `op_mapping_parent_sg` | `util/util_op_mapping` | `data_parent_pg`, `data_parent_sg` |
| polar-axial conversion, symmetric and anti-symmetric products | `create_binary/create_group.character_plus` | `data_character_table` and the symmetry operations |

| file | content | origin | created by | used by |
|---|---|---|---|---|
| [`data_basis.py`](#data-basis-py) | orbital and basis definitions: atomic orbitals, spin, LM/JM/LMs/Lgs basis, bra-ket blocks, LaTeX names, Wannier90 aliases, cubic/hexagonal orbitals (internal basis of multipolar harmonics), harmonics keys, max ranks and rank names | new | - | `util/util_atomic_multipole`, `util/util_cluster_vector_basis`, `util/util_harmonics_multipole`, `create_binary/create_atomic_multipole_group`, `create_binary/create_info`, `create_pdf/create_pdf_group` |
| [`data_character_table.py`](#data-character-table-py) | character tables of all 47 point groups | MultiPie Ver.2.4.0 | - | `util/util_atomic_multipole`, `util/util_cluster_samb`, `util/util_group_data`, `util/util_harmonics`, `util/util_harmonics_multipole`, `util/util_harmonics_symbol`, `util/util_representation_matrix`, `create_binary/create_group` |
| [`data_compatibility_relation.py`](#data-compatibility-relation-py) | compatibility relations of the Oh and D6h subgroups | MultiPie Ver.2.4.0 | - | `util/util_harmonics_mapping`, `create_binary/create_info` |
| [`data_group_mpg.py`](#data-group-mpg-py) | group information of the magnetic point groups | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `create_binary/create_group` |
| [`data_group_msg.py`](#data-group-msg-py) | group information of the magnetic space groups | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `create_binary/create_group` |
| [`data_group_pg.py`](#data-group-pg-py) | group information of the point groups | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `util/util_symmetry_operation`, `create_binary/create_group`, `create_pdf/create_pdf_harmonics` |
| [`data_group_sg.py`](#data-group-sg-py) | group information of the space groups | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `util/util_symmetry_operation`, `create_binary/create_group` |
| [`data_parent_pg.py`](#data-parent-pg-py) | parent point groups of each point group | new | `others/creation/create_parent_pg_mapping.py` | `util/util_op_mapping` |
| [`data_parent_sg.py`](#data-parent-sg-py) | parent space groups of each space group | new | - | `util/util_op_mapping` |
| [`data_response_tensor.py`](#data-response-tensor-py) | multipoles in the response tensors and their cartesian components | MultiPie Ver.2.4.0 | - | `create_binary/create_info` |
| [`data_root_cluster.py`](#data-root-cluster-py) | Oh and D6h clusters with the representative site and vector | new | `others/creation/create_root_cluster.py` | `util/util_bond_cluster_basis`, `util/util_cluster_basis`, `create_binary/create_info` |
| [`data_top_harmonics.py`](#data-top-harmonics-py) | harmonics of the top groups (Oh, D6h) | MultiPie Ver.2.4.0 | - | `util/util_harmonics_mapping`, `util/util_top_harmonics` |
| [`data_wyckoff_bond_sym.py`](#data-wyckoff-bond-sym-py) | Wyckoff bonds of the point and space groups with the bond symmetry | MultiPie Ver.2.4.0 | `others/creation/create_wyckoff_bond.py` | `create_binary/create_group` |
| [`data_wyckoff_mpg.py`](#data-wyckoff-mpg-py) | Wyckoff sites of the magnetic point groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |
| [`data_wyckoff_msg.py`](#data-wyckoff-msg-py) | Wyckoff sites of the magnetic space groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |
| [`data_wyckoff_site.py`](#data-wyckoff-site-py) | Wyckoff sites of the point and space groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |

## data_basis.py

Orbital and basis definitions: atomic orbitals, spin, LM/JM/LMs/Lgs basis, bra-ket blocks, LaTeX names, Wannier90 aliases, cubic/hexagonal orbitals (internal basis of multipolar harmonics), harmonics keys, max ranks and rank names.

| variable | type (size) | description |
|---|---|---|
| `atomic_orbital_basis` | list (4) | regular atomic orbital basis order (same as tesseral functions up to rank 3). |
| `atomic_spin_basis` | list (2) | regular atomic spin basis order. |
| `atomic_lm_basis` | list (4) | regular atomic LM basis order (M: descending). |
| `atomic_jm_basis` | list (4) | regular atomic JM basis order (J: ascending, M: descending). |
| `atomic_lms_basis` | list (4) | regular atomic spinful basis order. |
| `atomic_lgs_basis` | list (4) | regular atomic spinful basis order of the orbitals (Lgs) |
| `atomic_rank_block` | list (10) | bra-ket rank block, (L1,L2) for matrix element of atomic SAMB (uppper triangular block only). |
| `atomic_latex_name` | dict (18) | convert to regular atomic latex name. |
| `wannier90_orbital` | dict (9) | alias for Wannier90: wannier90 -> regular atomic orbital. |
| `cubic_orbital_basis` | list (4) | cubic basis order (multipolar PG harmonics internal basis). |
| `hexagonal_orbital_basis` | list (4) | hexagonal basis order (multipolar PG harmonics internal basis). |
| `cubic_basis` | list (4) | spinful basis order of cubic_orbital_basis |
| `hexagonal_basis` | list (4) | spinful basis order of hexagonal_orbital_basis |
| `orbital_latex_name` | dict (41) | convert to latex name up to for internal basis of multipolar PG harmonics and response tensor (add prefix Q/G/T/M). |
| `harmonics_key` | dict (41) | convert from name to key, ("tesseral", (l,m)), ("Oh", (rank, irrep, mul, comp)), ("D6h", (rank, irrep, mul, comp)). negative m represents S, otherwise C. |
| `max_rank` | int | max. rank of spherical harmonics. |
| `max_s_rank` | int | max. internal rank of spherical harmonics. |
| `rank_name` | dict (14) | rank name. |

## data_character_table.py

Character tables of all 47 point groups.

| variable | type (size) | description |
|---|---|---|
| `character` | dict (47) | { PG: ( [[SO in conjugacy class]], { irrep: [character] } ) } |

```text
Character table for all point groups.
    { point_group: ( [[symmetry operations in conjugacy class]], { irrep: [character] } ) }

- wp = exp(+2pi i/3), wm = exp(-2pi i/3).
- polar-axial conversion and symmetric/anti-symmetric products are computed from the characters
  (create_binary.create_group.character_plus).
```

## data_compatibility_relation.py

Compatibility relations of the Oh and D6h subgroups.

| variable | type (size) | description |
|---|---|---|
| `compatibility_relation` | dict (47) | { PG: (top PG, { (top irrep, comp): ((polar irrep, comp, coeff), (axial irrep, comp, coeff)) }) } |

```text
Compatibility relation for Oh and D6h subgroups.
    { point group: (top_pg, compatibility_relation }.
    compatibility_relation = { (top_pg_irrep, top_pg_comp): ((irrep_polar, comp, coeff.), (irrep_axial, comp, coeff.)) }.

- real expressions are used for complex group, and corresponding complex version is the point group with "-c".
```

## data_group_mpg.py

Group information of the magnetic point groups.

| variable | type (size) | description |
|---|---|---|
| `mpg_info` | dict (157) | { MPG_id: field } |

```text
Group information for all magnetic point groups.

- mpg_info.
    { MPG_id: field }
    MPG_id: MPG:{PG}.{no}.{ID} (PG:1-37, ID:1-155).
- mpg_tag, mpg_id, mpg_id_set: derived from mpg_info in multipiegen.util.util_group_data.

- field.
    - tag: MPG_id (PG.no.ID).
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, MPG -> MSG -> SG.
    - MPG: = id.
    - MSG: associated MSG, 1st MSG in the same MPGs.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III.
    - SO: symmetry operations.
```

## data_group_msg.py

Group information of the magnetic space groups.

| variable | type (size) | description |
|---|---|---|
| `msg_info` | dict (1651) | { MSG_id: field } |

```text
Group information for all magnetic space groups.

- msg_info.
    { MSG_id: field }
    MSG_id: MSG:{SG}.{no} (SG:1-230).
- msg_tag, msg_id, msg_id_set: derived from msg_info in multipiegen.util.util_group_data.

- field.
    - tag: BNS_id (SG.no).
    - BNS: Belov-Neronova-Smirnova (international) notation in LaTeX.
    - OG: Opechowski-Guccione notation in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, unique.
    - MPG: associated MPG, unique.
    - MSG: = id.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III, IV.
    - SO: symmetry operations.
```

## data_group_pg.py

Group information of the point groups.

| variable | type (size) | description |
|---|---|---|
| `pg_info` | dict (47) | { PG_id: field } |

```text
Group information for all point groups.

- pg_info.
    { PG_id: field }
    PG_id: PG:{ID} (ID:1-37; complex:38-47).
- pg_tag, pg_id, pg_id_set: derived from pg_info in multipiegen.util.util_group_data.

- field.
    - tag: Schoenflies in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: = id.
    - SG: associated SG, 1st SG in the same PGs.
    - MPG: associated MPG, 1st MPG with type II in the same PGs.
    - MSG: associated MSG, PG -> MPG -> MSG.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
```

## data_group_sg.py

Group information of the space groups.

| variable | type (size) | description |
|---|---|---|
| `sg_info` | dict (230) | { SG_id: field } |

```text
Group information for all space groups.

- sg_info.
    { SG_id: field }
    SG_id: SG:{ID} (ID:1-230).
- sg_tag, sg_id, sg_id_set: derived from sg_info in multipiegen.util.util_group_data.

- field.
    - tag: Schoenflies symbol in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting comment.
    - PG: associated PG, unique.
    - SG: = id.
    - MPG: associated MPG, SG -> MSG -> MPG.
    - MSG: associated MSG, 1st MSG with type II in the same SGs.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
```

## data_parent_pg.py

Parent point groups of each point group.

Created by `others/creation/create_parent_pg_mapping.py`.

| variable | type (size) | description |
|---|---|---|
| `data_parent_pg` | dict (35) | { child PG_id: { parent PG_id: ([[parent SO]], R) } } |
| `data_representative_parent_pg` | dict (35) | { child PG_id: parent_PG_id } |

```text
Parent group information for all point groups.

- data_parent_pg.
    { child PG_id: { parent PG_id: ([[parent SO]], R) } }.
    SO: 1-based, [parent SO] in order of child-PG SO.
    R: rotation matrix for converting principal axis (fractional).

- data_representative_parent_pg (first parent PG in data_parent_pg).
    { child PG_id: parent_PG_id }.
```

## data_parent_sg.py

Parent space groups of each space group.

| variable | type (size) | description |
|---|---|---|
| `data_parent_sg` | dict (216) | { child SG_id: { parent SG_id: ([[parent SO]], R) } }. [only one parent group compatible with representative parent PG] |
| `data_representative_parent_sg` | dict (216) | { child SG_id: parent_SG_id } |

```text
Parent group information for all space groups.

- data_parent_sg.
    { child SG_id: { parent SG_id: ([[parent SO]], R) } }. [only one parent group compatible with representative parent PG]
    SO: 1-based, [parent SO] in order of child-SG SO.
    R: rotation + translation matrix for converting principal axis (fractional).

- data_representative_parent_sg.
    { child SG_id: parent_SG_id }.
```

## data_response_tensor.py

Multipoles in the response tensors and their cartesian components.

| variable | type (size) | description |
|---|---|---|
| `rt_multipole_list` | dict (2) | list of multipoles in response tensor. |
| `rt_cartesian_multipole` | dict (2) | cartesian component in terms of multipoles. |

```text
Note:
definition of gv is changed. (it is multiplied by (-1)).
change name.
```

## data_root_cluster.py

Oh and D6h clusters with the representative site and vector.

Created by `others/creation/create_root_cluster.py`.

| variable | type (size) | description |
|---|---|---|
| `oh_cluster_cart_site0` | str | representative site in cartesian coordinate (normalized when necessary), radius of scaled one is given by 1/scale. d6h with principal axis [111] is equivanlent to oh. |
| `d6h_cluster_cart_site0` | str | representative site of D6h cluster. (same as oh for principal axis [111]) |
| `oh_cluster_cart_vector0` | str | representative vector in cartesian coordinate (normalized when necessary). it should not be parallel nor perpendicular to site0 vector. d6h with principal axis [111] is equivanlent to oh. |
| `d6h_cluster_cart_vector0` | str | representative vector of D6h cluster. (same as oh for principal axis [111]) |
| `cluster_scale` | int | scale factor for space group purpose. site0 is scaled as site0/scale. |
| `oh_cluster_cart` | str | positions of oh0 cluster in cartesian coordinate (48 points). principal axis [001]. (same as general Wyckoff points) |
| `d6h_cluster_cart` | str | positions of D6h cluster in cartesian coordinate (24 points). principal axis [001]. (same as general Wyckoff points) |

```text
Definition of root cluseter.

- oh_cluster_cart_site0 (str): representative site of Oh cluster.
- d6h_cluster_cart_site0 (str): representative site of D6h cluster. (same as oh for principal axis [111])
- oh_cluster_cart_vector0 (str): representative vector of Oh cluster.
- d6h_cluster_cart_vector0 (str): representative vector of D6h cluster. (same as oh for principal axis [111])
- cluster_scale (int): scale factor (site -> site/scale).
- oh_cluster_cart (str): sites of Oh cluster (first=[x,y,z]).
- d6h_cluster_cart (str): sites of D6h cluster (first=[x,y,z]).

note:
- general points of Oh (oh0).
    - oh0 : gs @ [x,y,z]

- general points of D6h (d6h0).
    - d6h0 : g0s @ [x,y,z], where g0s in principal axis [111] ([-110],[-1-12])
```

## data_top_harmonics.py

Harmonics of the top groups (Oh, D6h).

| variable | type (size) | description |
|---|---|---|
| `top_harmonics_tag` | dict (2) | tags of the top harmonics |
| `top_harmonics` | dict (480) | top harmonics (cartesian expression) of Oh and D6h |

```text
Note:
definitions of ("D6h", 4, "E2g", 2, 0) and ("D6h", 4, "E2g", 2, 1) are changed. (both is multiplied by (-1))
```

## data_wyckoff_bond_sym.py

Wyckoff bonds of the point and space groups with the bond symmetry.

Created by `others/creation/create_wyckoff_bond.py`.

| variable | type (size) | description |
|---|---|---|
| `wyckoff_bond_data` | dict (277) | { group: { Wyckoff: (vector, center, mapping, bond_symmetry) } } |

```text
Wyckoff bond for all space and point groups (fractional coordinate, conventional cell, no plus set).
- wyckoff_bond_data.
    { group: { Wyckoff: (vector, center, mapping, bond_symmetry) } }
```

## data_wyckoff_mpg.py

Wyckoff sites of the magnetic point groups.

| variable | type (size) | description |
|---|---|---|
| `wyckoff_mpg_data` | dict (157) | { MPG_id: { Wyckoff: (position, site-symmetry) } } |

```text
Wyckoff positions for all magnetic point groups (fractional coordinate).
    MPG_id: { Wyckoff: (position, site-symmetry) }
```

## data_wyckoff_msg.py

Wyckoff sites of the magnetic space groups.

| variable | type (size) | description |
|---|---|---|
| `wyckoff_msg_data` | dict (1651) | { MSG_id: { Wyckoff: (position, site-symmetry) } } |

```text
Wyckoff positions for all magnetic space groups (fractional coordinate).
    MSG_id: { Wyckoff: (position, site-symmetry) }

- alpha in wyckoff letter is replaced by A.
```

## data_wyckoff_site.py

Wyckoff sites of the point and space groups.

| variable | type (size) | description |
|---|---|---|
| `wyckoff_site_data` | dict (267) | { PG/SG: { Wyckoff: (position, site-symmetry) } } (point groups first) |

```text
Wyckoff positions for all point groups (fractional coordinate) and space groups (fractional coordinate, conventional
cell).
    wyckoff_site_data = { point_group/space_group: { Wyckoff: (position, site-symmetry) } } (point groups first)

- only for standard setting, (unique b axes; for space groups origin choice 2, hexagonal axes, abc).
- positions are in the same order of ITA (BCS).
```
