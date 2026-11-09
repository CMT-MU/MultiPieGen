# Examples

Sample code for each feature, in `docs/examples`. Run each file from anywhere after `pip install -e .`, e.g.
`python docs/examples/group.py`.

## cluster_basis.py

Load the orthonormal Oh/D6h cluster basis, inspect keys and values, expand a function on the cluster.

```{literalinclude} ../examples/cluster_basis.py
:language: python
:start-at: import numpy as np
```

## cluster_vector_basis.py

Load the orthonormal vector cluster basis (internal rank s=1, polar/axial internal basis) of the Oh/D6h clusters, expand a vector field.

```{literalinclude} ../examples/cluster_vector_basis.py
:language: python
:start-at: import numpy as np
```

## bond_cluster_basis.py

Load the orthonormal bond cluster basis of the Oh/D6h clusters (vector cluster basis contracted with the bond vectors and re-orthonormalized), keyed by PGMultipoleType as MultiPie's bond_a.

```{literalinclude} ../examples/bond_cluster_basis.py
:language: python
:start-at: from the vector cluster basis (cluster_vector_basis.pkl).
```

## group.py

Data of each group (group.pkl): symmetry operations, product, mapping to Oh/D6h, character, Wyckoff site and bond, active multipoles.

```{literalinclude} ../examples/group.py
:language: python
:start-at: from multipiegen.util.util_group import load_group
```

## representation_matrix.py

Representation matrices of a point group from the harmonics up to rank 11, and the check that all harmonics of an irrep give the same matrices.

```{literalinclude} ../examples/representation_matrix.py
:language: python
:start-at: from multipiegen.util.util_representation_matrix import (
```

## response_tensor.py

Response tensors (rank 0-4) of a magnetic point group in terms of active multipoles (response_tensor.pkl).

```{literalinclude} ../examples/response_tensor.py
:language: python
:start-at: from multipiegen.util.util_group import load_group
```

## harmonics_symbol.py

Polar and axial symbols of harmonics side by side for all 47 point groups (as in the harmonics PDF).

```{literalinclude} ../examples/harmonics_symbol.py
:language: python
:start-at: from multipiegen.util.util_harmonics_symbol import harmonics_symbol_table, polar_to_axial, multipole_symbol
```

## info.py

Overall information of groups (info.pkl): group id/tag/name, orbital basis and basis functions, root cluster.

```{literalinclude} ../examples/info.py
:language: python
:start-at: from multipiegen.util.util_info import load_info
```
