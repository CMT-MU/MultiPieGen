# Binary data

The binary data in `multipiegen/binary` are created from the data tables by the modules in `multipiegen/create_binary`. They use the format of MultiPie Ver.2's `BinaryManager` (gzip-compressed pickle with a header), so they can be read without MultiPieGen, e.g. `BinaryManager("info", topdir=...)`. `str(data)` shows the header, the description of the dictionary keys given in each section below.

| file | content | top-level keys | create | load |
|---|---|---|---|---|
| [`harmonics.pkl`](#harmonics-pkl) | expression and u-matrix of the harmonics of all PG (rank <= 11), looked up by a multipole key | 47 keys: `PG:1`, `PG:2`, ..., `PG:47` | `python -m multipiegen.create_binary.create_harmonics` | `multipiegen.util.util_harmonics.load_harmonics()` |
| [`info.pkl`](#info-pkl) | overall information of all groups, harmonics, basis, root cluster, response tensor | `id_set`, `tag`, `id`, `group_name`, `character`, `root_cluster`, `harmonics`, `response_tensor` | `python -m multipiegen.create_binary.create_info` | `multipiegen.util.util_info.load_info()` |
| [`group.pkl`](#group-pkl) | data of each PG, SG, MPG, MSG: symmetry operations, character, Wyckoff sites and bonds | 2085 keys: `PG:1`, `PG:2`, ..., `MSG:230.149` | `python -m multipiegen.create_binary.create_group` | `multipiegen.util.util_group.load_group()` |
| [`response_tensor.pkl`](#response-tensor-pkl) | response tensors (rank 0-4) of all MPG | 157 keys: `MPG:1.1.1`, `MPG:1.2.2`, ..., `MPG:37.5.155` | `python -m multipiegen.create_binary.create_response_tensor` | `multipiegen.util.util_response_tensor.load_response_tensor()` |
| [`multipolar_harmonics.pkl`](#multipolar-harmonics-pkl) | multipolar harmonics (s=1-3, l<=3) of all PG | 47 keys: `PG:1`, `PG:2`, ..., `PG:47` | `python -m multipiegen.create_binary.create_multipolar_harmonics` | `multipiegen.util.util_harmonics_multipole.load_multipolar_harmonics()` |
| [`atomic_multipole.pkl`](#atomic-multipole-pkl) | matrices of atomic multipoles in the LMs/JML/LM basis | `lms`, `jml`, `lm` | `python -m multipiegen.create_binary.create_atomic_multipole` | `multipiegen.util.util_atomic_multipole.load_atomic_multipole()` |
| [`atomic_multipole_group.pkl`](#atomic-multipole-group-pkl) | atomic multipole SAMB of all PG in the JML/Lgs/Lg basis | 47 keys: `PG:1`, `PG:2`, ..., `PG:47` | `python -m multipiegen.create_binary.create_atomic_multipole_group` | `multipiegen.util.util_atomic_multipole.load_atomic_multipole_group()` |
| [`cluster_scalar_basis.pkl`](#cluster-scalar-basis-pkl) | orthonormal (scalar) cluster basis of the Oh and D6h clusters | `Oh`, `D6h` | `python -m multipiegen.create_binary.create_cluster_scalar_basis` | `multipiegen.util.util_cluster_basis.load_cluster_basis()` |
| [`cluster_vector_basis.pkl`](#cluster-vector-basis-pkl) | orthonormal vector cluster basis (s=1) of the Oh and D6h clusters | `Oh`, `D6h` | `python -m multipiegen.create_binary.create_cluster_vector_basis` | `multipiegen.util.util_cluster_vector_basis.load_cluster_vector_basis()` |
| [`cluster_bond_basis.pkl`](#cluster-bond-basis-pkl) | orthonormal bond cluster basis of the Oh and D6h clusters | `Oh`, `D6h` | `python -m multipiegen.create_binary.create_cluster_bond_basis` | `multipiegen.util.util_bond_cluster_basis.load_bond_cluster_basis()` |
| [`cluster_samb.pkl`](#cluster-samb-pkl) | cluster SAMB of each Wyckoff site and bond of the real PG and SG | 267 keys: `PG:1`, `PG:2`, ..., `SG:230` | `python -m multipiegen.create_binary.create_cluster_samb` | `multipiegen.util.util_cluster_samb.load_cluster_samb()` |

## harmonics.pkl

Expression and u-matrix of the harmonics of all PG (rank <= 11), looked up by a multipole key.

- create: `python -m multipiegen.create_binary.create_harmonics`
- load: `multipiegen.util.util_harmonics.load_harmonics()`
- top-level keys: 47 keys: `PG:1`, `PG:2`, ..., `PG:47`
- size: 604 KB

Description of the dictionary keys (header of the binary):

```text
* Harmonics of all PG (rank l <= 11): cartesian expression and u-matrix for each multipole key.
- PG_id (str): (dict) harmonics.
  - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,sympy), ndarray(2l+1,dim,sympy)) [expression] and [u-matrix].
NOTE:
  - X = Q/G, p = -1, s = k = 0, x = "q" in the key. Look up the harmonics of a key of the other binaries by its
    (X,l,Gamma,n), with T -> Q and M -> G (multipiegen.util.util_harmonics.harmonics_expression, harmonics_u_matrix).
  - n is the multiplicity of the harmonics (X,l,Gamma), 1-based; "-1" represents no multiplicity.
  - expression: polynomial in x, y, z of each component; u-matrix: <m|component> for m = l, l-1, ..., -l.
  - the harmonics is c times a top harmonics of Oh or D6h (data_harmonics_mapping).
  - sorted in order of [Gamma, l, Q/G, n].
```

## info.pkl

Overall information of all groups, harmonics, basis, root cluster, response tensor.

- create: `python -m multipiegen.create_binary.create_info`
- load: `multipiegen.util.util_info.load_info()`
- top-level keys: `id_set`, `tag`, `id`, `group_name`, `character`, `root_cluster`, `harmonics`, `response_tensor`
- size: 137 KB

Description of the dictionary keys (header of the binary):

```text
* Information of all groups.
- "id_set" (str): (dict) info. of group's IDs.
  - "PG/SG/MPG/MSG" (str):  (dict) ID set.
    - PG: all/crystal/complex/irrep.
    - SG: all/crystal/PG.
    - MPG: all/crystal/PG/type.
    - MSG: all/crystal/PG/SG/MPG/type.
- "tag" (str): (dict) from id to tag.
- "id" (str):  (dict) from tag to id.
- "group_name" (str): (dict) name list of groups (data_group_name_list).
  - "PG/SG/MPG/MSG" (str): (dict) name info. of each group.
    - id (str): (dict) "dir" (directory name), "md" (markdown file name), "no" (number), "tag" (tag),
                       "name" (name in LaTeX), "PG/SG/MPG/MSG" (associated group id).
- "character" (str): (dict) info. of character.
  - "alias" (str): (dict) alias.
    - "wp/wm" (str): (sympy) definition of wp and wm, exp(2 pi i/3) or exp(-2 pi i/3).
  - "compatibility_relation" (str): (dict) compatibility relation.
    - (PG_tag1,PG_tag2) (str,str): ([(str,str)]) list of compatibility relation, (irrep1,irrep2).
- "root_cluster" (str): (dict) info. of root cluster.
  - "point_group" (str): ([str]) list of PG_tag for root cluster, ["Oh", "D6h"].
  - "scale" (str): (int) scale factor, scaled site = site/scale.
  - "oh_site0_cart" (str): (ndarray(3,sympy)) representative site of Oh cluster (cartesian = fractional for cubic, r=1).
  - "oh_scaled_site0_cart" (str): (ndarray(3,sympy)) scaled representative site of Oh cluster (cartesian, r=1/scale).
  - "oh_vector0_cart" (str): (ndarray(3,sympy)) representative vector of Oh cluster (cartesian = fractional for cubic, v=1).
  - "oh_site0_frac" (str): (ndarray(3,sympy)) representative site of Oh cluster (fractional for cubic, same as cartesian).
  - "oh_scaled_site0_frac" (str): (ndarray(3,sympy)) scaled representative site of Oh cluster (fractional for cubic, same as cartesian).
  - "oh_vector0_frac" (str): (ndarray(3,sympy)) representative vector of Oh cluster (fractional for cubic, same as cartesian).
  - "d6h_site0_cart" (str): (ndarray(3,sympy)) representative site of D6h cluster (cartesian, r=1).
  - "d6h_scaled_site0_cart" (str): (ndarray(3,sympy)) scaled representative site of D6h cluster (cartesian, r=1/scale).
  - "d6h_vector0_cart" (str): (ndarray(3,sympy)) representative vector of D6h cluster (cartesian, v=1).
  - "d6h_site0_frac" (str): (ndarray(3,sympy)) representative site of D6h cluster (fractional for hexagonal).
  - "d6h_scaled_site0_frac" (str): (ndarray(3,sympy)) scaled representative site of D6h cluster (fractional for hexagonal).
  - "d6h_vector0_frac" (str): (ndarray(3,sympy)) representative vector of D6h cluster (fractional for hexagonal).
- "harmonics" (str): (dict) info. of harmonics.
  - "max_rank" (str): (int) max rank of harmonics.
  - "max_s_rank" (str): (int) max internal rank.
  - "variable_site" (str): (ndarray(3,sympy)) variable of site, [x,y,z].
  - "variable_vector" (str): (ndarray(3,sympy)) variable of (bond) vector, [X,Y,Z].
  - "internal_basis" (str): (dict) internal basis of multipolar harmonics.
    - "cubic/hexagonal" (str): (dict) basis (data_basis.cubic_orbital_basis/hexagonal_orbital_basis).
      - (s,x) (int,str): ([sympy]) list of basis, x = "q" (polar)/"g" (axial).
  - "atomic_basis" (str): (dict) atomic basis.
    - "spinless" (str): (dict) spinless basis.
      - "lm/lg/cubic/hexagonal" (str): (dict) basis.
        - L (int): ([str]) list of basis, M or gamma.
    - "spinful" (str): (dict) spinful basis.
      - "jml/lms/lgs/cubic/hexagonal" (str): (dict) basis.
        - L (int): ([str]) list of basis, (J,M) or (M,s) or (gamma,s).
    - "spin" (str): ([str]) spin basis, ["u","d"].
    - "rank_block" (str): ([(int,int)]) bra-ket rank block (L1,L2) (upper triangular block only).
  - "basis_function" (str): (dict) orbital basis function.
    - name (str): ((sympy,ndarray(2L+1,sympy))) cartesian expression and u-matrix, <m|gamma>.T (m: descending).
  - "harmonics_key" (str): (dict) key of top harmonics of orbital.
    - name (str): ((str,tuple)) ("tesseral", (l,m)) or ("Oh"/"D6h", (rank, irrep, mul, comp)), m<0 for S.
  - "wannier90" (str): (dict) alias for Wannier90.
    - tag (str): (str) name corresponding to Wannier90 tag.
  - "tesseral" (str): (dict) tesseral tag.
    - name (str): ((int,int,str)) tesseral tag (L,M,"c/s").
  - "rank_name" (str): (dict) conversion between rank and name.
    - name/rank (str/int): (int/str) rank/name.
  - "latex" (str): (dict) LaTeX name.
    - "atomic" (str): (dict) name -> LaTeX of atomic orbital and spin.
    - "orbital" (str): (dict) name -> LaTeX suffix of internal basis (add prefix Q/G/T/M).
  NOTE:
    - see "multipiegen/data/data_basis.py" for detail.
- "response_tensor" (str): (dict) info. of response tensor.
  - "multipole" (str): (dict) multipole info.
    - "cubic/hexagonal" (str): ([str]) list of multipoles.
  - "cartesian_multipole" (str): (dict) linear combination info.
    - "cubic/hexagonal" (str): { comp (str): ([(int,str)]) }, relation between cartesian component and multipole.
```

## group.pkl

Data of each PG, SG, MPG, MSG: symmetry operations, character, Wyckoff sites and bonds.

- create: `python -m multipiegen.create_binary.create_group`
- load: `multipiegen.util.util_group.load_group()`
- top-level keys: 2085 keys: `PG:1`, `PG:2`, ..., `MSG:230.149`
- size: 2,114 KB

Description of the dictionary keys (header of the binary):

```text
* Information for all PG ("PG").
- PG_id (str): (dict) id, "PG:{ID}".
  - "info" (str): (namedtuple) info.
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
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "generator" (str): ([str]) list of generator SO tag.
    - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
    - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
      - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
    - "product" (str): (dict) product table.
      - (SO_tag1,SO_tag2) (str,str): (str) SO_tag of product, SO1 SO2.
    - "mapping" (str): (dict) mapping of SO to the Oh/D6h cluster (data_op_mapping).
      - "top" (str): (str) top group, "Oh" or "D6h".
      - "mapping" (str): ([[int]]) mapping[coset][i] = SO number (from 1) of top group for i-th SO of group.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SOs of SG and PG (without translational part) are in the same order.
      - SO matrices are in the same order of "tag".
      - monopole (s=0) and dipole (s=1) basis are [1] and [x,y,z].
      - mapping[0] is the group itself, and mapping[k] (k>0) is the left coset g_k * group.
  - "wyckoff" (str): (dict) Wyckoff site/bond of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
        - "bond" (str): ([str]) list of bw_tag.
    - "bond" (str): (dict) info. of Wyckoff bond.
      - bw_tag (str): (dict) bond Wyckoff tag.
        - "symmetry" (str): (str) bond symmetry.
        - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional).
        - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
    NOTE:
      - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
      - Wyckoff sites are in the same order of ITA (BCS).
      - representative site/vector are scaled site0/vector0 of the root cluster (info["root_cluster"]).
  - "character" (str): (dict) character table of group.
    - "conjugacy" (str): ([[str]]) list of SO_tag in each conjugacy class.
    - "dimension" (str): (dict) dimension of irrep.
      - irrep_tag (str): (int) dimension.
    - "table" (str): (dict) character table (first SO in each conjugacy class).
      - irrep_tag (str): (ndarray(n,sympy) list of character.
    - "table_full" (str): (dict) character table for all SO.
      - irrep_tag (str): (ndarray(n,sympy) list of character.
    - "polar_axial_conversion" (str): (dict) conversion between polar and axial.
      - irrep_tag (str): (str) converted irrep_tag.
    - "symmetric_product" (str): (dict) symmetric product decomposition.
      - (irrep_tag,irrep_tag) (str,str): [(int,str)] list of (n, irrep_tag).
    - "anti_symmetric_product" (str): (dict) anti-symmetric product decomposition.
      - irrep_tag (str): ([(int,str)]) list of (n,irrep_tag).
    NOTE:
      - wp = exp(+2pi i/3), wm = exp(-2pi i/3) in "table", and they are substituted in "table_full".
      - for complex character, use PG_tag with "-c", otherwise without "-c", real version is given.
  NOTE:
    - complex PG ("-c") has the same "symmetry_operation" and "wyckoff" as the real PG.

* Information for all SG ("SG").
- SG_id (str): (dict) id, "SG:{ID}".
  - "info" (str): (namedtuple) info.
    - tag: Schoenflies in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting comment.
    - PG: associated PG, unique.
    - SG: = id.
    - MPG: associated MPG, SG -> PG -> MPG.
    - MSG: associated MSG, 1st MSG with type II in the same SGs.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "generator" (str): ([str]) list of generator SO tag.
    - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
    - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "cartesian_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, primitive).
    - "fractional_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, primitive).
    - "plus_set" (str): (ndarray(n,3,sympy)) list of plus set vector (fractional, conventional).
    - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
      - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
    - "mapping" (str): (dict) mapping of SO to the Oh/D6h cluster (that of the PG, data_op_mapping).
      - "top" (str): (str) top group, "Oh" or "D6h".
      - "mapping" (str): ([[int]]) mapping[coset][i] = SO number (from 1) of top group for i-th SO of group.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SOs of SG and PG (without translational part) are in the same order.
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site/bond of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
        - "primitive" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, primitive).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
        - "bond" (str): ([str]) list of bw_tag.
    - "bond" (str): (dict) info. of Wyckoff bond.
      - bw_tag (str): (dict) bond Wyckoff tag.
        - "symmetry" (str): (str) bond symmetry.
        - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, conventional, no plus set).
        - "primitive" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, primitive).
        - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
    NOTE:
      - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
      - Wyckoff sites are in the same order of ITA (BCS).

* Information for all MPG ("MPG").
- MPG_id (str): (dict) id, "MPG:{PG}.{no}.{ID}".
  - "info" (str): (namedtuple) info.
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
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
    - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
  - "active_multipole" (list): active MPs list, [str] (multipoles in info["response_tensor"]["multipole"]).

* Information for all MSG ("MSG").
- MSG_id (str): (dict) id, "MSG:{SG}.{no}".
  - "info" (str): (namedtuple) info.
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
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
    - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.

* all sympy values are expanded with rationalized denominators.
```

## response_tensor.pkl

Response tensors (rank 0-4) of all MPG.

- create: `python -m multipiegen.create_binary.create_response_tensor`
- load: `multipiegen.util.util_response_tensor.load_response_tensor()`
- top-level keys: 157 keys: `MPG:1.1.1`, `MPG:1.2.2`, ..., `MPG:37.5.155`
- size: 159 KB

Description of the dictionary keys (header of the binary):

```text
* Response tensors of all MPG (up to rank 4).
- MPG_id (str): (dict) response tensors, "MPG:{PG}.{no}.{ID}".
  - (X, rank, opt) (str,int,str): ((ndarray(sympy), dict)) tensor matrix and expression.
    - X: tensor type, "Q" (polar), "G" (axial), "T" (magnetic polar), "M" (magnetic axial).
    - rank: tensor rank, 0-4.
    - opt: "" (rank 0,1), "s/a" (rank 2), "ss/sa/a" (rank 3), "sss/ssa/aas/aaa/sa/as/s/a/t" (rank 4).
    - matrix: tensor components, e.g. Q_{xyz} (sympy symbol), 0 for vanishing components.
      - rank 0: (1,1), rank 1: (1,3) [x,y,z].
      - rank 2 (s/a): (3,3) [x,y,z] x [x,y,z].
      - rank 3: (6,3) [xx,yy,zz,yz,zx,xy] x [x,y,z] (ss/sa), (3,3) [yz,zx,xy] x [x,y,z] (a).
      - rank 4: (6,6) (sss/ssa), (3,3) (aas/aaa), (6,3) (sa), (3,6) (as), (6,9) (s), (3,9) (a), (10,3) (t),
                with pairs v=[xx,yy,zz,yz,zx,xy], a=[yz,zx,xy], u=[xx,yy,zz,yz,zx,xy,zy,xz,yx],
                and triples t=[xxx,yyy,zzz,yyz,zzx,xxy,yzz,zxx,xyy,xyz].
    - expression: {component symbol: expression in terms of active multipoles (LaTeX symbols)}.
NOTE:
  - multipole symbol: \mathbb{X}_{sub}^{sup(n)} with the LaTeX names of info["harmonics"]["latex"]["orbital"],
    (n) labels the independent cartesian decomposition.
  - same as MultiPie Ver.2 Group.response_tensor_all except for the multipole names.
```

## multipolar_harmonics.pkl

Multipolar harmonics (s=1-3, l<=3) of all PG.

- create: `python -m multipiegen.create_binary.create_multipolar_harmonics`
- load: `multipiegen.util.util_harmonics_multipole.load_multipolar_harmonics()`
- top-level keys: 47 keys: `PG:1`, `PG:2`, ..., `PG:47`
- size: 1,187 KB

Description of the dictionary keys (header of the binary):

```text
* Multipolar harmonics of all PG (s = 1..3, rank l <= 3) with polar (q) and axial (g) internal basis.
- PG_id (str): (dict) multipolar harmonics.
  - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,2s+1,sympy), ndarray(dim,sympy))
    [internal coefficient] (polynomial in x,y,z multiplying each internal basis function) and [multipole expression]
    (with the internal symbols), for each component.
NOTE:
  - the expression of the harmonics (X,l,Gamma,n) of a key is in harmonics.pkl.
  - internal basis: orbitals of data_basis.cubic_orbital_basis[s] (cubic) or hexagonal_orbital_basis[s] (trigonal,
    hexagonal), with the symbols Q/G + orbital_latex_name (info["harmonics"]["internal_basis"]).
  - x = q if (X = Q and s+k even) or (X = G and s+k odd), otherwise g.
  - n is the multiplicity of the harmonics (X,l,Gamma) (its xyz polynomial, not renumbered); p numbers the same
    (x, l, X, Gamma, n, k) as 1, 2, ...; "-1" represents no multiplicity / no additional multiplicity.
  - sorted in order of [q/g, s, Gamma, l, Q/G, k, n, p].
```

## atomic_multipole.pkl

Matrices of atomic multipoles in the LMs/JML/LM basis.

- create: `python -m multipiegen.create_binary.create_atomic_multipole`
- load: `multipiegen.util.util_atomic_multipole.load_atomic_multipole()`
- top-level keys: `lms`, `jml`, `lm`
- size: 79 KB

Description of the dictionary keys (header of the binary):

```text
* Matrix of atomic multipole (rank l <= 7, s = 0, 1) in the orbitals of L = 0..3.
- "lms/jml/lm" (str): (dict) matrix in LMs/JML/LM basis (info["harmonics"]["atomic_basis"]).
  - (L1,L2) (int,int): (dict) atomic multipole data for each bra-ket block.
    - (X,l,s,k,x) (SphericalMultipoleType): (ndarray(2l+1,2(2L1+1),2(2L2+1),sympy)) matrix for m = l, l-1, ..., -l.
NOTE:
  - half size of matrix for spinless "lm" (s = 0, up-spin part of "lms").
  - x = "q" for all.
  - multipoles are sorted in order of [Q/G/T/M, s, k, l].
```

## atomic_multipole_group.pkl

Atomic multipole SAMB of all PG in the JML/Lgs/Lg basis.

- create: `python -m multipiegen.create_binary.create_atomic_multipole_group`
- load: `multipiegen.util.util_atomic_multipole.load_atomic_multipole_group()`
- top-level keys: 47 keys: `PG:1`, `PG:2`, ..., `PG:47`
- size: 6,161 KB

Description of the dictionary keys (header of the binary):

```text
* Matrix of atomic multipole SAMB for all PG (rank l <= 7, s = 0, 1, orbitals of L = 0..3).
- PG_id (str): (dict) matrix data.
  - "jml/lgs/lg" (str): (dict) matrix of atomic multipole in (J,M;L)/(L,gamma,s)/(L,gamma) basis.
    - (L1,L2) (int,int): (dict) atomic multipole data for each bra-ket block.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,2(2L1+1),2(2L2+1),sympy)) matrix for each component.
NOTE:
  - the expression of the harmonics (X,l,Gamma,n) of a key is in harmonics.pkl.
  - (L,gamma) basis: orbitals of info["harmonics"]["atomic_basis"]["spinless"]["lg"].
  - half size of matrix for spinless "lg" (s = 0, up-spin part of "lgs" times sqrt(2)).
  - each component is normalized to the Frobenius norm 1 (diagonal block) or 1/sqrt(2) (off-diagonal block).
  - n is the multiplicity of the harmonics (X,l,Gamma) (its xyz polynomial, not renumbered); p = -1; x = "q".
  - multipoles are sorted in order of [Q/G/T/M, s, k, l, Gamma, n, p].
```

## cluster_scalar_basis.pkl

Orthonormal (scalar) cluster basis of the Oh and D6h clusters.

- create: `python -m multipiegen.create_binary.create_cluster_scalar_basis`
- load: `multipiegen.util.util_cluster_basis.load_cluster_basis()`
- top-level keys: `Oh`, `D6h`
- size: 5 KB

Description of the dictionary keys (header of the binary):

```text
Orthonormal (scalar) cluster basis of Oh and D6h clusters (binary: cluster_scalar_basis.pkl).

The cluster is made by substituting the representative site (site0) into the general points
of the Oh (48 points) or D6h (24 points) cluster. The value of each real top harmonics at
each cluster point gives a cluster vector. Different irreps (and components) are orthogonal
by symmetry. Within each (irrep, component), the vectors are orthogonalized by the
Gram-Schmidt process from lower rank (then lower multiplicity), and linearly dependent ones
are dropped. The result is an orthonormal basis of 48 (Oh) or 24 (D6h) vectors.

- cluster_scalar_basis (dict)::

    { "Oh/D6h": {
        "site0": [x0,y0,z0] (sympy),
        "cluster": cluster points, ndarray (N,3) (sympy),
        "basis": { PGMultipoleType: ndarray (dim,N) (sympy) },  # orthonormal, sum_i v_i^2 = 1 per component.
        "dropped": [ PGMultipoleType ],  # harmonics linearly dependent on lower ones (whole multiplets).
    } }

    PGMultipoleType = (X, l, Gamma, n, p, s, k, x) as the site SAMB of MultiPie Ver.2: X = "Q", rank l, irrep,
        multiplicity n (1-based, -1: unique) of the top harmonics, p = -1, s = 0, k = 0, x = "q".
        All components of a multiplet are in the value (rows, in the order of components).
```

## cluster_vector_basis.pkl

Orthonormal vector cluster basis (s=1) of the Oh and D6h clusters.

- create: `python -m multipiegen.create_binary.create_cluster_vector_basis`
- load: `multipiegen.util.util_cluster_vector_basis.load_cluster_vector_basis()`
- top-level keys: `Oh`, `D6h`
- size: 59 KB

Description of the dictionary keys (header of the binary):

```text
Orthonormal vector cluster basis (internal rank s=1) of Oh and D6h clusters (binary: cluster_vector_basis.pkl).

This is an intermediate basis to create the bond cluster basis (util_bond_cluster_basis, cluster_bond_basis.pkl); whether it is used
by itself is not decided yet.

Same procedure as the cluster basis (util_cluster_basis), for the vector harmonics. s=0 is the cluster basis itself.

- Vector harmonics. For a real top harmonics T (top key (l, irrep, mul, comp), u-matrix U) and k = -1, 0, 1,
  the vector harmonics is sum_m U[m] O_{lm}^{(s=1)}(k) (util_multipolar_harmonics, orbital rank l+k), a 3-vector
  whose components are in the order of the internal basis [px, py, pz] (data_basis.cubic_orbital_basis[1] and
  hexagonal_orbital_basis[1]).
- Type. For a polar internal basis (x="q"), it is polar (X="Q") if s+k is even and axial (X="G") if s+k is odd;
  for an axial internal basis (x="g"), Q and G are interchanged. A polar one has the irrep and component of T
  (coefficient 1); an axial one has the irrep, component and coefficient c of data_harmonics_mapping, G(key) = c T.
- Cluster vector. The values at the cluster points (N,3) (see util_cluster_basis.create_cluster).
- Orthonormalization. Within each (irrep, component), Gram-Schmidt in the order of rank l, then k = -1, 0, 1,
  then multiplicity; linearly dependent ones are dropped. The result is 3N = 144 (Oh) or 72 (D6h) vectors for each x.
  The linear dependence is first found numerically, then the kept vectors are orthonormalized exactly.

- cluster_vector_basis (dict)::

    { "Oh/D6h": {
        "site0": [x0,y0,z0] (sympy),
        "cluster": cluster points, ndarray (N,3) (sympy),
        "internal_basis": ["px","py","pz"],
        "basis": { PGMultipoleType: ndarray (dim,N,3) (sympy) },  # orthonormal, sum_{i,a} v_{ia}^2 = 1 per component.
        "dropped": [ PGMultipoleType ],  # vector harmonics linearly dependent on lower ones (whole multiplets).
    } }

    PGMultipoleType = (X, l, Gamma, n, p, s, k, x) as the vector SAMB of MultiPie Ver.2: type Q/G, rank, irrep,
        multiplicity (-1: unique), add_multiplicity p = -1 (multiplets are distinguished by k), s = 1, k = -1, 0, 1,
        internal type x = q/g. All components of a multiplet are in the value (rows, in the order of components).
```

## cluster_bond_basis.pkl

Orthonormal bond cluster basis of the Oh and D6h clusters.

- create: `python -m multipiegen.create_binary.create_cluster_bond_basis`
- load: `multipiegen.util.util_bond_cluster_basis.load_bond_cluster_basis()`
- top-level keys: `Oh`, `D6h`
- size: 7 KB

Description of the dictionary keys (header of the binary):

```text
Orthonormal bond cluster basis of Oh and D6h clusters (binary: cluster_bond_basis.pkl).

- Bond vectors. b_i is the general point i of the Oh or D6h cluster with x,y,z replaced by vector0 (unit vector,
  info["root_cluster"]["oh/d6h_vector0_cart"]).
- Bond cluster vectors. For each vector cluster basis v_key (util_cluster_vector_basis, (N,3)), the N-component
  vector f_key(i) = v_key(i) . b_i, with the same key. They are 3N vectors in the N-dimensional space.
- Orthonormalization. Within each (irrep, component), Gram-Schmidt in the order of rank, then k = -1, 0, 1, then
  multiplicity; linearly dependent ones (including vanishing ones) are dropped. The result is N = 48 (Oh) or
  24 (D6h) vectors for each internal basis type x = q, g. The linear dependence is first found numerically,
  then the kept vectors are orthonormalized exactly.
- Symmetry. f_key transforms as a scalar function on the cluster, (g f)(i) = s f(j) with R_g^{-1} p_i = p_j,
  s = 1 (x = q) or det R_g (x = g), with the representation matrices of its irrep.
- Key. As the anti-symmetric bond cluster SAMB (bond_a) of MultiPie Ver.2, a kept multiplet is a scalar function on the
  bonds keyed by PGMultipoleType(X, l, Gamma, n, p, s=0, k=0, x), with all components in the value. (X, l, Gamma, n)
  are those of the vector cluster basis. p (add_multiplicity) numbers the multiplets with the same
  (X, l, Gamma, n) [components] in the Gram-Schmidt order (k = -1, 0, 1) as 1, 2, ..., and is -1 if there is only
  one. x is the internal type of the vector cluster basis (MultiPie Ver.2 has x = q only).

- cluster_bond_basis (dict)::

    { "Oh/D6h": {
        "site0": [x0,y0,z0] (sympy),
        "vector0": [X0,Y0,Z0] (sympy),
        "cluster": cluster points, ndarray (N,3) (sympy),
        "bond": bond vectors b_i, ndarray (N,3) (sympy),
        "basis": { PGMultipoleType: ndarray (dim, N) (sympy) },  # orthonormal, rows = components.
        "source": { PGMultipoleType: PGMultipoleType },  # multiplet of the vector cluster basis (s=1).
        "dropped": [ PGMultipoleType ],  # multiplets of the vector cluster basis dropped.
    } }

    PGMultipoleType = (X, l, Gamma, n, p, s, k, x): type Q/G, rank, irrep, multiplicity (-1: unique),
        add_multiplicity (-1: unique), s = 0, k = 0, internal type of the vector cluster basis q/g.
    The vector cluster basis (cluster_vector_basis.pkl) is an intermediate basis to create this one.
```

## cluster_samb.pkl

Cluster SAMB of each Wyckoff site and bond of the real PG and SG.

- create: `python -m multipiegen.create_binary.create_cluster_samb`
- load: `multipiegen.util.util_cluster_samb.load_cluster_samb()`
- top-level keys: 267 keys: `PG:1`, `PG:2`, ..., `SG:230`
- size: 377 KB

Description of the dictionary keys (header of the binary):

```text
* Cluster SAMB for all real PG and SG.
- PG_id/SG_id (str): (dict) SAMB data.
  - "site" (str): (dict) site-cluster SAMB data.
    - sw_tag (str): (dict) SAMB data at site Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy)) SAMB for each component.
  - "bond_s" (str): (dict) symmetric-bond-cluster SAMB data.
    - bw_tag (str): (dict) SAMB data at bond Wyckoff, as "site".
  - "bond_a" (str): (dict) anti-symmetric-bond-cluster SAMB data.
    - bw_tag (str): (dict) SAMB data at bond Wyckoff, as "site".
NOTE:
  - the expression of the harmonics (X,l,Gamma,n) of a key is in harmonics.pkl.
  - created from the general-point SAMB of each group, i.e., cluster_scalar_basis ("site", "bond_s") and
    cluster_bond_basis (x = "q", "bond_a") of Oh or D6h restricted by the operation mapping (coset 0), and averaged
    by the mapping of each Wyckoff site and bond (with the sign of the bond reversal for "bond_a").
  - keys of the top group with the irrep and component of the group; n is the multiplicity of the harmonics (X,l,Gamma).
  - "-1" in multiplicity represents no multiplicity.
  - each component is normalized on the ns sites (bonds).
  - multipoles are sorted in order of [Gamma, l, n, p] ("site", "bond_s") or [Gamma, l, k, Q/G, n, p] ("bond_a").
```
