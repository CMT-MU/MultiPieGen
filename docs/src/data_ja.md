---
orphan: true
---

# データテーブル

英語版は[こちら](data.md)です。

`multipiegen/data` のデータテーブルは、生成するすべてのデータの元です。データテーブルは手で管理しています。一部は、最初に `others/creation` のスクリプトで作ります（スクリプトは自分と同じ場所に書き出し、その結果を `multipiegen/data` にコピーします）。「由来」は、MultiPie Ver.2.4.0 から取ったテーブル（一部は変更あり）では「MultiPie Ver.2.4.0」、それ以外では「新規」です。MultiPie Ver.2.4.0 から取ったテーブルは一度だけコピーし、それ以降はここで管理しています。生成では MultiPie Ver.2 のデータを読むことはなく、その一般的な関数とユーティリティだけを使います。

これらのテーブルから導けるデータは、保存せずに `multipiegen` の中で計算します：

| データ | 計算する場所 | 元のデータ |
|---|---|---|
| `*_tag`, `*_id`, `*_id_set`, `group_name_list` | `util/util_group_data` | `data_group_pg/sg/mpg/msg` の `*_info` |
| `top_harmonics_symmetry`, `top_harmonics_u_matrix` | `util/util_top_harmonics` | `data_top_harmonics`（方域の部分は点群多極子調和関数から） |
| `data_harmonics_mapping` | `util/util_harmonics_mapping` | `data_compatibility_relation`, `data_top_harmonics.top_harmonics_tag` |
| `op_mapping_d6h_oh_cluster`, `op_mapping_parent_sg` | `util/util_op_mapping` | `data_parent_pg`, `data_parent_sg` |
| 極性と軸性の変換、対称積と反対称積 | `create_binary/create_group.character_plus` | `data_character_table` と対称操作 |

| ファイル | 内容 | 由来 | 作成 | 使用 |
|---|---|---|---|---|
| [`data_basis.py`](#data-basis-py) | 軌道と基底の定義：原子軌道、スピン、LM/JM/LMs/Lgs 基底、bra-ket ブロック、LaTeX の名前、Wannier90 の別名、立方晶・六方晶の軌道（点群多極子調和関数の内部基底）、調和関数のキー、最大ランクとランクの名前 | 新規 | - | `util/util_atomic_multipole`, `util/util_cluster_vector_basis`, `util/util_harmonics_multipole`, `create_binary/create_atomic_multipole_group`, `create_binary/create_info`, `create_pdf/create_pdf_group` |
| [`data_character_table.py`](#data-character-table-py) | 47 点群すべての指標表 | MultiPie Ver.2.4.0 | - | `util/util_atomic_multipole`, `util/util_cluster_samb`, `util/util_group_data`, `util/util_harmonics`, `util/util_harmonics_multipole`, `util/util_harmonics_symbol`, `util/util_representation_matrix`, `create_binary/create_group` |
| [`data_compatibility_relation.py`](#data-compatibility-relation-py) | Oh と D6h の部分群の適合関係 | MultiPie Ver.2.4.0 | - | `util/util_harmonics_mapping`, `create_binary/create_info` |
| [`data_group_mpg.py`](#data-group-mpg-py) | 磁気点群の情報 | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `create_binary/create_group` |
| [`data_group_msg.py`](#data-group-msg-py) | 磁気空間群の情報 | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `create_binary/create_group` |
| [`data_group_pg.py`](#data-group-pg-py) | 点群の情報 | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `util/util_symmetry_operation`, `create_binary/create_group`, `create_pdf/create_pdf_harmonics` |
| [`data_group_sg.py`](#data-group-sg-py) | 空間群の情報 | MultiPie Ver.2.4.0 | - | `util/util_group_data`, `util/util_symmetry_operation`, `create_binary/create_group` |
| [`data_parent_pg.py`](#data-parent-pg-py) | 各点群の親の点群 | 新規 | `others/creation/create_parent_pg_mapping.py` | `util/util_op_mapping` |
| [`data_parent_sg.py`](#data-parent-sg-py) | 各空間群の親の空間群 | 新規 | - | `util/util_op_mapping` |
| [`data_response_tensor.py`](#data-response-tensor-py) | 応答テンソルに現れる多極子と、そのデカルト成分 | MultiPie Ver.2.4.0 | - | `create_binary/create_info` |
| [`data_root_cluster.py`](#data-root-cluster-py) | 代表サイトと代表ベクトルをもつ Oh と D6h のクラスター | 新規 | `others/creation/create_root_cluster.py` | `util/util_bond_cluster_basis`, `util/util_cluster_basis`, `create_binary/create_info` |
| [`data_top_harmonics.py`](#data-top-harmonics-py) | トップの群（Oh、D6h）の調和関数 | MultiPie Ver.2.4.0 | - | `util/util_harmonics_mapping`, `util/util_top_harmonics` |
| [`data_wyckoff_bond_sym.py`](#data-wyckoff-bond-sym-py) | 点群と空間群の Wyckoff ボンドと、ボンドの対称性 | MultiPie Ver.2.4.0 | `others/creation/create_wyckoff_bond.py` | `create_binary/create_group` |
| [`data_wyckoff_mpg.py`](#data-wyckoff-mpg-py) | 磁気点群の Wyckoff サイト | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |
| [`data_wyckoff_msg.py`](#data-wyckoff-msg-py) | 磁気空間群の Wyckoff サイト | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |
| [`data_wyckoff_site.py`](#data-wyckoff-site-py) | 点群と空間群の Wyckoff サイト | MultiPie Ver.2.4.0 | - | `create_binary/create_group` |

## data_basis.py

軌道と基底の定義：原子軌道、スピン、LM/JM/LMs/Lgs 基底、bra-ket ブロック、LaTeX の名前、Wannier90 の別名、立方晶・六方晶の軌道（点群多極子調和関数の内部基底）、調和関数のキー、最大ランクとランクの名前です。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `atomic_orbital_basis` | list (4) | regular atomic orbital basis order (same as tesseral functions up to rank 3). |
| `atomic_spin_basis` | list (2) | regular atomic spin basis order. |
| `atomic_lm_basis` | list (4) | regular atomic LM basis order (M: descending). |
| `atomic_jm_basis` | list (4) | regular atomic JM basis order (J: ascending, M: descending). |
| `atomic_lms_basis` | list (4) | regular atomic spinful basis order. |
| `atomic_lgs_basis` | list (4) | 原子軌道の標準のスピンあり基底の並び（Lgs） |
| `atomic_rank_block` | list (10) | bra-ket rank block, (L1,L2) for matrix element of atomic SAMB (uppper triangular block only). |
| `atomic_latex_name` | dict (18) | convert to regular atomic latex name. |
| `wannier90_orbital` | dict (9) | alias for Wannier90: wannier90 -> regular atomic orbital. |
| `cubic_orbital_basis` | list (4) | cubic basis order (multipolar PG harmonics internal basis). |
| `hexagonal_orbital_basis` | list (4) | hexagonal basis order (multipolar PG harmonics internal basis). |
| `cubic_basis` | list (4) | cubic_orbital_basis のスピンあり基底の並び |
| `hexagonal_basis` | list (4) | hexagonal_orbital_basis のスピンあり基底の並び |
| `orbital_latex_name` | dict (41) | convert to latex name up to for internal basis of multipolar PG harmonics and response tensor (add prefix Q/G/T/M). |
| `harmonics_key` | dict (41) | convert from name to key, ("tesseral", (l,m)), ("Oh", (rank, irrep, mul, comp)), ("D6h", (rank, irrep, mul, comp)). negative m represents S, otherwise C. |
| `max_rank` | int | max. rank of spherical harmonics. |
| `max_s_rank` | int | max. internal rank of spherical harmonics. |
| `rank_name` | dict (14) | rank name. |

## data_character_table.py

47 点群すべての指標表です。

| 変数 | 型（サイズ） | 説明 |
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

Oh と D6h の部分群の適合関係です。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `compatibility_relation` | dict (47) | { PG: (top PG, { (top irrep, comp): ((polar irrep, comp, coeff), (axial irrep, comp, coeff)) }) } |

```text
Compatibility relation for Oh and D6h subgroups.
    { point group: (top_pg, compatibility_relation }.
    compatibility_relation = { (top_pg_irrep, top_pg_comp): ((irrep_polar, comp, coeff.), (irrep_axial, comp, coeff.)) }.

- real expressions are used for complex group, and corresponding complex version is the point group with "-c".
```

## data_group_mpg.py

磁気点群の情報です。

| 変数 | 型（サイズ） | 説明 |
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

磁気空間群の情報です。

| 変数 | 型（サイズ） | 説明 |
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

点群の情報です。

| 変数 | 型（サイズ） | 説明 |
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

空間群の情報です。

| 変数 | 型（サイズ） | 説明 |
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

各点群の親の点群です。

`others/creation/create_parent_pg_mapping.py` で作ります。

| 変数 | 型（サイズ） | 説明 |
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

各空間群の親の空間群です。

| 変数 | 型（サイズ） | 説明 |
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

応答テンソルに現れる多極子と、そのデカルト成分です。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `rt_multipole_list` | dict (2) | list of multipoles in response tensor. |
| `rt_cartesian_multipole` | dict (2) | cartesian component in terms of multipoles. |

```text
Note:
definition of gv is changed. (it is multiplied by (-1)).
change name.
```

## data_root_cluster.py

代表サイトと代表ベクトルをもつ Oh と D6h のクラスターです。

`others/creation/create_root_cluster.py` で作ります。

| 変数 | 型（サイズ） | 説明 |
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

トップの群（Oh、D6h）の調和関数です。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `top_harmonics_tag` | dict (2) | トップの調和関数のタグ |
| `top_harmonics` | dict (480) | Oh と D6h のトップの調和関数（デカルト座標の式） |

```text
Note:
definitions of ("D6h", 4, "E2g", 2, 0) and ("D6h", 4, "E2g", 2, 1) are changed. (both is multiplied by (-1))
```

## data_wyckoff_bond_sym.py

点群と空間群の Wyckoff ボンドと、ボンドの対称性です。

`others/creation/create_wyckoff_bond.py` で作ります。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `wyckoff_bond_data` | dict (277) | { group: { Wyckoff: (vector, center, mapping, bond_symmetry) } } |

```text
Wyckoff bond for all space and point groups (fractional coordinate, conventional cell, no plus set).
- wyckoff_bond_data.
    { group: { Wyckoff: (vector, center, mapping, bond_symmetry) } }
```

## data_wyckoff_mpg.py

磁気点群の Wyckoff サイトです。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `wyckoff_mpg_data` | dict (157) | { MPG_id: { Wyckoff: (position, site-symmetry) } } |

```text
Wyckoff positions for all magnetic point groups (fractional coordinate).
    MPG_id: { Wyckoff: (position, site-symmetry) }
```

## data_wyckoff_msg.py

磁気空間群の Wyckoff サイトです。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `wyckoff_msg_data` | dict (1651) | { MSG_id: { Wyckoff: (position, site-symmetry) } } |

```text
Wyckoff positions for all magnetic space groups (fractional coordinate).
    MSG_id: { Wyckoff: (position, site-symmetry) }

- alpha in wyckoff letter is replaced by A.
```

## data_wyckoff_site.py

点群と空間群の Wyckoff サイトです。

| 変数 | 型（サイズ） | 説明 |
|---|---|---|
| `wyckoff_site_data` | dict (267) | { PG/SG: { Wyckoff: (position, site-symmetry) } }（点群が先） |

```text
Wyckoff positions for all point groups (fractional coordinate) and space groups (fractional coordinate, conventional
cell).
    wyckoff_site_data = { point_group/space_group: { Wyckoff: (position, site-symmetry) } } (point groups first)

- only for standard setting, (unique b axes; for space groups origin choice 2, hexagonal axes, abc).
- positions are in the same order of ITA (BCS).
```
