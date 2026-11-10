# Examples

Sample code for the features of MultiPieGen. Run each file from anywhere after `pip install -e .`.

| file | content |
|---|---|
| [cluster_basis.py](cluster_basis.py) | load the orthonormal (scalar) Oh/D6h cluster basis (cluster_scalar_basis.pkl), keyed by PGMultipoleType, expand a function on the cluster |
| [cluster_vector_basis.py](cluster_vector_basis.py) | load the orthonormal vector cluster basis (internal rank s=1, polar/axial internal basis) of the Oh/D6h clusters, keyed by PGMultipoleType, expand a vector field |
| [bond_cluster_basis.py](bond_cluster_basis.py) | load the orthonormal bond cluster basis of the Oh/D6h clusters (vector cluster basis contracted with the bond vectors and re-orthonormalized), keyed by PGMultipoleType as MultiPie's bond_a |
| [group.py](group.py) | data of each group (group.pkl): symmetry operations, product, mapping to Oh/D6h, character, Wyckoff site and bond, active multipoles |
| [representation_matrix.py](representation_matrix.py) | representation matrices of a point group from the harmonics up to rank 11, and the check that all harmonics of an irrep give the same matrices |
| [response_tensor.py](response_tensor.py) | response tensors (rank 0-4) of a magnetic point group in terms of active multipoles (response_tensor.pkl) |
| [multipolar_harmonics.py](multipolar_harmonics.py) | multipolar harmonics (internal rank s=1,2,3, rank l<=3) with polar/axial internal basis of the data_basis orbitals (multipolar_harmonics.pkl) |
| [atomic_multipole.py](atomic_multipole.py) | atomic multipole matrices (atomic_multipole.pkl) and atomic multipole SAMB of a point group in the jml/lgs/lg basis, keyed by PGMultipoleType (atomic_multipole_group.pkl) |
| [harmonics_symbol.py](harmonics_symbol.py) | polar and axial symbols of harmonics side by side for all 47 point groups (as in the harmonics PDF) |
| [info.py](info.py) | overall information of groups (info.pkl): group id/tag/name, orbital basis and basis functions, root cluster |
