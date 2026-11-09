# Examples

Sample code for the features of MultiPieGen. Run each file from anywhere after `pip install -e .`.

| file | content |
|---|---|
| [cluster_basis.py](cluster_basis.py) | load the orthonormal Oh/D6h cluster basis, inspect keys and values, expand a function on the cluster |
| [cluster_vector_basis.py](cluster_vector_basis.py) | load the orthonormal vector cluster basis (internal rank s=1, polar/axial internal basis) of the Oh/D6h clusters, expand a vector field |
| [group.py](group.py) | data of each group (group.pkl): symmetry operations, product, mapping to Oh/D6h, character, Wyckoff site and bond, active multipoles |
| [representation_matrix.py](representation_matrix.py) | representation matrices of a point group from the harmonics up to rank 11, and the check that all harmonics of an irrep give the same matrices |
| [response_tensor.py](response_tensor.py) | response tensors (rank 0-4) of a magnetic point group in terms of active multipoles (response_tensor.pkl) |
| [harmonics_symbol.py](harmonics_symbol.py) | polar and axial symbols of harmonics side by side for all 47 point groups (as in the harmonics PDF) |
| [info.py](info.py) | overall information of groups (info.pkl): group id/tag/name, orbital basis and basis functions, root cluster |
