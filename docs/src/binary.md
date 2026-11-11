# Binary data

The binary data in `multipiegen/binary` are created from the data tables by the modules in `multipiegen/create_binary`.
They use the format of MultiPie Ver.2's `BinaryManager` (gzip-compressed pickle with a header), so they can be read without
MultiPieGen. `str(data)` shows the description below.

| file | create | load |
|---|---|---|
| `info.pkl` | `python -m multipiegen.create_binary.create_info` | `multipiegen.util.util_info.load_info()` |
| `group.pkl` | `python -m multipiegen.create_binary.create_group` | `multipiegen.util.util_group.load_group()` |
| `atomic_multipole.pkl` | `python -m multipiegen.create_binary.create_atomic_multipole` | `multipiegen.util.util_atomic_multipole.load_atomic_multipole()` |
| `atomic_multipole_group.pkl` | `python -m multipiegen.create_binary.create_atomic_multipole_group` | `multipiegen.util.util_atomic_multipole.load_atomic_multipole_group()` |
| `cluster_samb.pkl` | `python -m multipiegen.create_binary.create_cluster_samb` | `multipiegen.util.util_cluster_samb.load_cluster_samb()` |
| `multipolar_harmonics.pkl` | `python -m multipiegen.create_binary.create_multipolar_harmonics` | `multipiegen.util.util_harmonics_multipole.load_multipolar_harmonics()` |
| `response_tensor.pkl` | `python -m multipiegen.create_binary.create_response_tensor` | `multipiegen.util.util_response_tensor.load_response_tensor()` |
| `cluster_scalar_basis.pkl` | `python -m multipiegen.create_binary.create_cluster_scalar_basis` | `multipiegen.util.util_cluster_basis.load_cluster_basis()` |
| `cluster_vector_basis.pkl` | `python -m multipiegen.create_binary.create_cluster_vector_basis` | `multipiegen.util.util_cluster_vector_basis.load_cluster_vector_basis()` |
| `cluster_bond_basis.pkl` | `python -m multipiegen.create_binary.create_cluster_bond_basis` | `multipiegen.util.util_bond_cluster_basis.load_bond_cluster_basis()` |

## info.pkl

Overall information: group id/tag/name, character alias and compatibility relation, root cluster, orbital basis and
basis functions, response tensor.

```{literalinclude} ../../multipiegen/create_binary/create_info.py
:language: text
:start-after: h_info = """
:end-before: '"""'
```

## group.pkl

Data of all point groups (PG, 47), space groups (SG, 230), magnetic point groups (MPG, 157) and magnetic space groups
(MSG, 1651), keyed by group id.

```{literalinclude} ../../multipiegen/create_binary/create_group.py
:language: text
:start-after: h_group = """
:end-before: '"""'
```

## response_tensor.pkl

Response tensors (rank 0-4) of all magnetic point groups (157) in terms of active multipoles, keyed by MPG id.

```{literalinclude} ../../multipiegen/create_binary/create_response_tensor.py
:language: text
:start-after: h_response_tensor = """
:end-before: '"""'
```

## cluster_scalar_basis.pkl

Orthonormal cluster basis of the Oh (48 points) and D6h (24 points) root clusters.

```{literalinclude} ../../multipiegen/util/util_cluster_basis.py
:language: text
:start-after: '"""'
:end-before: '"""'
```

## cluster_vector_basis.pkl

Orthonormal vector cluster basis (internal rank s=1) of the Oh (144) and D6h (72) clusters, for polar (q) and axial (g)
internal basis. It is an intermediate basis to create `cluster_bond_basis.pkl`.

```{literalinclude} ../../multipiegen/util/util_cluster_vector_basis.py
:language: text
:start-after: '"""'
:end-before: '"""'
```

## cluster_bond_basis.pkl

Orthonormal bond cluster basis of the Oh (48) and D6h (24) clusters, from the vector cluster basis contracted with the
bond vectors (vector0) and re-orthonormalized.

```{literalinclude} ../../multipiegen/util/util_bond_cluster_basis.py
:language: text
:start-after: '"""'
:end-before: '"""'
```

## multipolar_harmonics.pkl

Multipolar harmonics (internal rank s = 1, 2, 3, rank l <= 3) of all 47 point groups with polar (q) and axial (g)
internal basis.

```{literalinclude} ../../multipiegen/util/util_harmonics_multipole.py
:language: text
:start-after: '"""'
:end-before: '"""'
```

## atomic_multipole.pkl, atomic_multipole_group.pkl

Matrices of atomic multipoles (rank l <= 7, s = 0, 1) in the s, p, d, f orbitals, and the atomic multipole SAMB of
all 47 point groups.

```{literalinclude} ../../multipiegen/util/util_atomic_multipole.py
:language: text
:start-after: '"""'
:end-before: '"""'
```

## cluster_samb.pkl

Cluster SAMB of each Wyckoff site and bond (site, symmetric bond, anti-symmetric bond) of the real point groups and
the space groups.

```{literalinclude} ../../multipiegen/util/util_cluster_samb.py
:language: text
:start-after: '"""'
:end-before: '"""'
```
