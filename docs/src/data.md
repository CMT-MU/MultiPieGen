# Data tables

The data tables in `multipiegen/data` are the source of all generated data. They are curated by hand; some of them are first created by a script in `others/creation` (it writes next to itself, and the result is copied into `multipiegen/data`). "origin" is "MultiPie Ver.2.4.0" for the tables taken from MultiPie Ver.2.4.0 (some of them modified), and "new" otherwise.

| file | content | origin | created by | used by |
|---|---|---|---|---|
| [`data_basis.py`](#data-basis-py) | orbital and basis definitions: atomic orbitals, spin, LM/JM/LMs/Lgs basis, bra-ket blocks, LaTeX names, Wannier90 aliases, cubic/hexagonal orbitals (internal basis of multipolar harmonics), harmonics keys | new | - | `util/util_atomic_multipole`, `util/util_cluster_vector_basis`, `util/util_harmonics_multipole`, `create_binary/create_atomic_multipole_group`, `create_binary/create_info`, `create_pdf/create_pdf_group` |
| [`data_character_table.py`](#data-character-table-py) | character tables of all 47 point groups | MultiPie Ver.2.4.0 | - | `util/util_atomic_multipole`, `util/util_cluster_samb`, `util/util_harmonics`, `util/util_harmonics_multipole`, `util/util_harmonics_symbol`, `util/util_representation_matrix`, `create_binary/create_group` |
| [`data_compatibility_relation.py`](#data-compatibility-relation-py) | compatibility relations of the Oh and D6h subgroups | MultiPie Ver.2.4.0 | - | `create_binary/create_info` |
| [`data_coset_so.py`](#data-coset-so-py) | symmetry operations of the cosets to the parent group | new | - | - |
| [`data_group_mpg.py`](#data-group-mpg-py) | group information of the magnetic point groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group`, `create_binary/create_info` |
| [`data_group_msg.py`](#data-group-msg-py) | group information of the magnetic space groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group`, `create_binary/create_info` |
| [`data_group_name_list.py`](#data-group-name-list-py) | names of all groups (directory, markdown, LaTeX names) | MultiPie Ver.2.4.0 | - | `create_binary/create_info` |
| [`data_group_pg.py`](#data-group-pg-py) | group information of the point groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group`, `create_binary/create_info`, `create_pdf/create_pdf_harmonics` |
| [`data_group_sg.py`](#data-group-sg-py) | group information of the space groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group`, `create_binary/create_info` |
| [`data_harmonics.py`](#data-harmonics-py) | definition of harmonics (old orbital names, use data_basis for names) | new | - | `create_binary/create_info` |
| [`data_harmonics_mapping.py`](#data-harmonics-mapping-py) | mapping of the harmonics of each point group to the top harmonics (Oh, D6h) | new | `others/creation/create_harmonics_mapping.py` | `util/util_atomic_multipole`, `util/util_cluster_samb`, `util/util_cluster_vector_basis`, `util/util_harmonics_multipole`, `util/util_harmonics_symbol`, `util/util_representation_matrix`, `create_binary/create_atomic_multipole_group`, `create_pdf/create_pdf_info` |
| [`data_op_mapping.py`](#data-op-mapping-py) | mapping of symmetry operations to the Oh/D6h cluster (PG), the parent SG (SG), and the root cluster | new | `others/creation/create_group_mapping.py` | `util/util_cluster_samb`, `create_binary/create_group` |
| [`data_parent_pg.py`](#data-parent-pg-py) | parent point groups of each point group | new | `others/creation/create_parent_pg_mapping.py` | - |
| [`data_parent_sg.py`](#data-parent-sg-py) | parent space groups of each space group | new | - | - |
| [`data_response_tensor.py`](#data-response-tensor-py) | multipoles in the response tensors and their cartesian components | MultiPie Ver.2.4.0 | - | `create_binary/create_info` |
| [`data_root_cluster.py`](#data-root-cluster-py) | root cluster, Oh and D6h clusters with the representative site and vector | new | `others/creation/create_root_cluster.py` | `util/util_bond_cluster_basis`, `util/util_cluster_basis`, `create_binary/create_info` |
| [`data_top_harmonics.py`](#data-top-harmonics-py) | harmonics of the top groups (Oh, D6h) | MultiPie Ver.2.4.0 | - | - |
| [`data_top_harmonics_info.py`](#data-top-harmonics-info-py) | symmetry (cartesian expression) and u-matrix of the top harmonics | new | `others/creation/create_top_harm.py` | `util/util_atomic_multipole`, `util/util_cluster_basis`, `util/util_cluster_vector_basis`, `util/util_harmonics_multipole`, `util/util_harmonics_symbol`, `util/util_representation_matrix`, `create_binary/create_info`, `create_pdf/create_pdf_info` |
| [`data_wyckoff_bond_sym.py`](#data-wyckoff-bond-sym-py) | Wyckoff bonds of the point and space groups with the bond symmetry | MultiPie Ver.2.4.0 | `others/creation/create_wyckoff_bond.py` | `create_binary/create_group` |
| [`data_wyckoff_mpg.py`](#data-wyckoff-mpg-py) | Wyckoff sites of the magnetic point groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |
| [`data_wyckoff_msg.py`](#data-wyckoff-msg-py) | Wyckoff sites of the magnetic space groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |
| [`data_wyckoff_site.py`](#data-wyckoff-site-py) | Wyckoff sites of the point and space groups | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |

## data_basis.py

Orbital and basis definitions: atomic orbitals, spin, LM/JM/LMs/Lgs basis, bra-ket blocks, LaTeX names, Wannier90 aliases, cubic/hexagonal orbitals (internal basis of multipolar harmonics), harmonics keys.

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

## data_character_table.py

Character tables of all 47 point groups.

| variable | type (size) | description |
|---|---|---|
| `character` | dict (47) | { PG: ( [[SO in conjugacy class]], { irrep: [character] } ) } |
| `character_plus` | dict (47) | { PG: { "polar_axial": { irrep: irrep }, "symmetric_product": { (irrep1, irrep2): [(n, irrep)] }, "anti_symmetric_product": { irrep: [(n, irrep)] } } } |

```text
Character table for all point groups.
    { point_group: ( [[symmetry operations in conjugacy class]], { irrep: [character] } ) }

- wp = exp(+2pi i/3), wm = exp(-2pi i/3).
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

## data_coset_so.py

Symmetry operations of the cosets to the parent group.

| variable | type (size) | description |
|---|---|---|
| `data_coset_so` | dict (251) | (parent, [SO]) } |

```text
Symmetry operations of coset.

- data_coset_so = { child: (parent, [SO]) }
```

## data_group_mpg.py

Group information of the magnetic point groups.

| variable | type (size) | description |
|---|---|---|
| `mpg_info` | dict (157) | { MPG_id: field } |
| `mpg_id_set` | dict (4) | { all/crystal/PG/type: [MPG_id] } |
| `mpg_tag` | dict (157) | { MPG_id: tag } |
| `mpg_id` | dict (157) | { tag: MPG_id } |

```text
Group information for all magnetic point groups.

- mpg_info.
    { MPG_id: field }
    MPG_id: MPG:{PG}.{no}.{ID} (PG:1-37, ID:1-155).
- mpg_id_set
    { all/crystal/PG/type: [MPG_id] }.
- mpg_tag.
    { MPG_id: tag }
- mpg_id.
    { tag: MPG_id }

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
| `msg_id_set` | dict (6) | { all/crystal/PG/SG/MPG/type: [MPG_id] } |
| `msg_tag` | dict (1651) | { MSG_id: tag } |
| `msg_id` | dict (1651) | { tag: MSG_id } |

```text
Group information for all magnetic space groups.

- msg_info.
    { MSG_id: field }
    MSG_id: MSG:{SG}.{no} (SG:1-230).
- msg_id_set
    { all/crystal/PG/SG/MPG/type: [MPG_id] }.
- msg_tag.
    { MSG_id: tag }
- msg_id.
    { tag: MSG_id }

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

## data_group_name_list.py

Names of all groups (directory, markdown, LaTeX names).

| variable | type (size) | description |
|---|---|---|
| `group_name_list` | dict (4) | { PG/SG/MPG/MSG: { id: (dir_name, md_name, no, tag, name, PG_id, SG_id, MPG_id, MSG_id) } } |

```text
Group name list for all magnetic space groups.

id : (dir_name, md_name, no, tag, name, PG_id, SG_id, MPG_id, MSG_id).
```

## data_group_pg.py

Group information of the point groups.

| variable | type (size) | description |
|---|---|---|
| `pg_info` | dict (47) | { PG_id: field } |
| `pg_id_set` | dict (4) | { all/crystal/complex/irrep: [PG_id]} |
| `pg_tag` | dict (47) | { PG_id: tag } |
| `pg_id` | dict (47) | { tag: PG_id } |

```text
Group information for all point groups.

- pg_info.
    { PG_id: field }
    PG_id: PG:{ID} (ID:1-37; complex:38-47).
- pg_id_set.
    { all/crystal/complex/irrep: [PG_id]}
    irrep: 1d, 2d, 3d.
- pg_tag.
    { PG_id: tag }
- pg_id.
    { tag: PG_id }

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
| `sg_id_set` | dict (3) | { all/crystal/PG: [SG_id]} |
| `sg_tag` | dict (230) | { SG_id: tag } |
| `sg_id` | dict (230) | { tag: SG_id } |

```text
Group information for all space groups.

- sg_info.
    { SG_id: field }
    SG_id: SG:{ID} (ID:1-230).
- sg_id_set.
    { all/crystal/PG: [SG_id]}
- sg_tag.
    { SG_id: tag }
- sg_id.
    { tag: SG_id }

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

## data_harmonics.py

Definition of harmonics (old orbital names, use data_basis for names).

| variable | type (size) | description |
|---|---|---|
| `max_rank` | int | max. rank of spherical harmonics. |
| `max_s_rank` | int | max. internal rank of spherical harmonics. |
| `tesseral_harmonics` | dict (144) | tesseral harmonics. |
| `rank_name` | dict (14) | rank name. |
| `orbital_rank` | dict (8) | rank and pre-letter. |
| `orbital_set` | dict (35) | definition of orbitals up to rank 4, (cartesian expression, <m|U|gamma>.T, name (LaTeX), tensor_name (LaTeX)). |
| `wannier90_orbital` | dict (9) | alias for Wannier90: wannier -> orbital name. |

```text
Definition of harmonics.

- max_rank (int): max rank of spherical harmonics.
- max_s_rank (int): max rank of internal spherical harmonics.
- tesseral_harmonics (dict): "Cl,m"/"Sl,m" -> expression.
- rank_name (dict): rank <-> name.
- orbital_rank (dict): rank <-> pre-letter.
- orbital_set (dict): name -> (cartesian expression, <m|U|gamma>.T, name, tensor_name).
- wannier90_orbital (dict): wannier90 name -> name.
```

## data_harmonics_mapping.py

Mapping of the harmonics of each point group to the top harmonics (Oh, D6h).

Created by `others/creation/create_harmonics_mapping.py`.

| variable | type (size) | description |
|---|---|---|
| `data_harmonics_mapping` | dict (47) | { PG: { "top": "Oh/D6h", "Q/G": { pg_harm_key: (coeff, top_harm_key) } } } |

```text
Mapping of harmonics key to each PG.

- data_harmonics_mapping.
    { PG: { "top": "Oh/D6h", "Q/G": { pg_harm_key: (coeff, top_harm_key) } } }.
    pg_harm_key = (rank, irrep, multiplicity (1-based), component (0-based)) [-1: no multiplicity/single component]
    top_harm_key = (rank, irrep, multiplicity, component)
```

## data_op_mapping.py

Mapping of symmetry operations to the Oh/D6h cluster (PG), the parent SG (SG), and the root cluster.

Created by `others/creation/create_group_mapping.py`.

| variable | type (size) | description |
|---|---|---|
| `op_mapping_d6h_oh_cluster` | dict (37) | operation mapping of general points of point groups to D6h and Oh cluster (1-based). |
| `op_mapping_parent_sg` | dict (230) | operation mapping of space groups to the top parent space group (221-230 or 191-194) through the chain of maximal subgroups (data_parent_sg), with the change of axes (1-based). |
| `op_mapping_root_cluster` | dict (267) | operation mapping of general points of point and space groups to root cluster (1-based). |

```text
Mapping of symmetry operations to top parent group.

- op_mapping_d6h_oh_cluster (dict): mapping of point groups to top Oh/D6h cluster (operation no. 1-based).
    - PG_tag : (top group, mapping), mapping[coset][parent op no. of child op. from 1 to g].
    - top group is child when child belongs to Oh/D6h.
    - the mapping of a space group to the Oh/D6h cluster is that of its point group (same point parts of operations).
- op_mapping_parent_sg (dict): mapping of space groups to top parent space group (operation no. 1-based).
    - SG_no : (top SG_no, mapping), mapping[coset][parent op no. of child op. from 1 to g].
    - top SG is 221-230 (Oh) or 191-194 (D6h), through the chain of maximal subgroups (data_parent_sg).
    - the axes may be changed along the chain, so the point part of a child operation may differ from that of the
      parent operation (e.g. 2[010] -> 2[100]); not a mapping to the Oh/D6h cluster.
    - top SG is child when child is the top SG.
- op_mapping_root_cluster (dict): mapping to root cluster (operation no. 1-based). [future extension]
    - PG_tag/SG_no : ("RC", mapping), mapping[coset][parent op no. of child op. from 1 to g].
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

Root cluster, Oh and D6h clusters with the representative site and vector.

Created by `others/creation/create_root_cluster.py`.

| variable | type (size) | description |
|---|---|---|
| `root_cluster_cart_site0` | str | representative site in cartesian coordinate (normalized when necessary), radius of scaled one is given by 1/scale. d6h with principal axis [111] is equivanlent to oh. |
| `oh_cluster_cart_site0` | str | representative site of Oh cluster |
| `d6h_cluster_cart_site0` | str | representative site of D6h cluster. (same as oh for principal axis [111]) |
| `root_cluster_cart_vector0` | str | representative vector in cartesian coordinate (normalized when necessary). it should not be parallel nor perpendicular to site0 vector. d6h with principal axis [111] is equivanlent to oh. |
| `oh_cluster_cart_vector0` | str | representative vector of Oh cluster |
| `d6h_cluster_cart_vector0` | str | representative vector of D6h cluster. (same as oh for principal axis [111]) |
| `cluster_scale` | int | scale factor for space group purpose. site0 is scaled as site0/scale. |
| `oh_cluster_cart` | str | positions of oh0 cluster in cartesian coordinate (48 points). principal axis [001]. (same as general Wyckoff points) |
| `d6h_cluster_cart` | str | positions of D6h cluster in cartesian coordinate (24 points). principal axis [001]. (same as general Wyckoff points) |
| `root_cluster_cart` | str | positions of root cluster in cartesian coordinate (24x4=96 points). contain general points of D6h (d6h0, d6h1, d6h2, d6h3) and Oh (oh0, oh1). |

```text
Definition of root cluseter.

- root_cluster_cart_site0 (str): representative site of root cluster.
- oh_cluster_cart_site0 (str): representative site of Oh cluster.
- d6h_cluster_cart_site0 (str): representative site of D6h cluster. (same as oh for principal axis [111])
- root_cluster_cart_vector0 (str): representative vector of root cluster.
- oh_cluster_cart_vector0 (str): representative vector of Oh cluster.
- d6h_cluster_cart_vector0 (str): representative vector of D6h cluster. (same as oh for principal axis [111])
- cluster_scale (int): scale factor (site -> site/scale).
- oh_cluster_cart (str): sites of Oh cluster (first=[x,y,z]).
- d6h_cluster_cart (str): sites of D6h cluster (first=[x,y,z]).
- root_cluster_cart (str): sites of root cluster (first=[x,y,z]). [future extension]

note:
- general points of Oh [oh0, oh1].
    - oh0 : gs @ [x,y,z]
    - oh1 : gs @ m[-1-12] @ [x,y,z], where m[-1-12] = [[2/3,-1/3,2/3],[-1/3,2/3,2/3],[2/3,2/3,-1/3]].

- general points of D6h with .
    - d6h0 : g0s @ [x,y,z], where g0s in principal axis [111] ([-110],[-1-12])
    - d6h1 : Rx @ d6h0, where Rx = [[1,0,0],[0,0,-1],[0,-1,0]] (principal axis [1-1-1] ([-1-10],[-11-2])).
    - d6h2 : Ry @ d6h0, where Ry = [[0,0,-1],[0,1,0],[-1,0,0]] (principal axis [-11-1] ([110],[1-1-2])).
    - d6h3 : Rz @ d6h0, where Rz = [[0,-1,0],[-1,0,0],[0,0,1]] (principal axis [-1-11] ([1-10],[112])).
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

## data_top_harmonics_info.py

Symmetry (cartesian expression) and u-matrix of the top harmonics.

Created by `others/creation/create_top_harm.py`.

| variable | type (size) | description |
|---|---|---|
| `top_harmonics_symmetry` | dict (3) | { "tesseral/Oh/D6h": { tess_key/top_key: symmetry (cartesian expression) } } |
| `top_harmonics_u_matrix` | dict (3) | { "tesseral/Oh/D6h": { tess_key/top_key: u-matrix <m|gamma> } } |

```text
Top harmonics infomation (polar) [axial can be obtained by interchanging g and u in irrep].

- top_harmonics_symmetry.
    { "tesseral/Oh/D6h": { tess_key/top_key: symmetry } }.
- top_harmonics_u_matrix, <m|gamma>.
    { "tesseral/Oh/D6h": { tess_key/top_key: u-matrix } }.

tess_key = (l, m), negative m is for S, otherwise for C.
top_key = (rank, irrep, multiplicity [1-based], component [0-based]).
-1 in multiplicity and component represent unique value.
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
| `wyckoff_site_pg` | dict (37) | { PG: { Wyckoff: (position, site-symmetry) } } |
| `wyckoff_site_sg` | dict (230) | { SG: { Wyckoff: (position, site-symmetry) } } |
| `wyckoff_site_data` | dict (267) | wyckoff_site_pg and wyckoff_site_sg together |

```text
Wyckoff positions for all point groups (fractional coordinate).
    point_group: { Wyckoff: (position, site-symmetry) }

- only for standard setting, (unique b axes).
- positions are in the same order of ITA (BCS).
```
