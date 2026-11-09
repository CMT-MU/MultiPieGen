# Examples

Sample code for the features of MultiPieGen. Run each file from anywhere after `pip install -e .`.

| file | content |
|---|---|
| [cluster_basis.py](cluster_basis.py) | load the orthonormal Oh/D6h cluster basis, inspect keys and values, expand a function on the cluster |
| [group.py](group.py) | data of each group (group.pkl): symmetry operations, product, mapping to Oh/D6h, character, Wyckoff site and bond, active multipoles |
| [representation_matrix.py](representation_matrix.py) | representation matrices of a point group from the harmonics up to rank 11, and the check that all harmonics of an irrep give the same matrices |
| [harmonics_symbol.py](harmonics_symbol.py) | polar and axial symbols of harmonics side by side for all 47 point groups (as in the harmonics PDF) |
| [info.py](info.py) | overall information of groups (info.pkl): group id/tag/name, orbital basis and basis functions, root cluster |
