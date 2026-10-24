Pipeline (run inside this folder; needs spglib, numpy; step4 needs pyxtal's wyckoff_list.csv, path set at its top):
 step1_enum_t_subgroups.py  -> subs.json     (t-subgroup types of each space group)
 step2_enum_k_subgroups.py  -> res_m4.json   (k-subgroups up to lattice index 4)
 step3_select.py            -> chosen.json, parent_sg_candidates.csv (parent choice)
 step4_ita_order.py         -> ita_ops.json  (ITA-ordered operations in spglib settings)
 step5_embeddings.py        -> emb_all.json  (embedding x_S = M x_G + d, canonical choice)
 step6_emit_tables.py       -> parent_op_correspondence.{csv,py}, parent_transformations.csv
Intermediate json files (subs.json, chosen.json, pgmap.json, ita_ops.json, emb_all.json) are included.

Note: operation strings and (M,d) are in ITA origin choice 2 for the 24 groups with two origin choices (inversion at the origin).
spglib itself uses origin choice 1 internally; step4 converts.
 step7_check_general_positions.py -> prints the general-position consistency check (needs sympy)
 step8_coset_decomposition.py -> parent_coset_decomposition.{py,csv} (left cosets of the parent ops by the child image, ITA op numbers)
