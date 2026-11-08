# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
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
"""

# ==================================================
# representative site in cartesian coordinate (normalized when necessary), radius of scaled one is given by 1/scale.
#   d6h with principal axis [111] is equivanlent to oh.
root_cluster_cart_site0 = "[1,5,6]"
oh_cluster_cart_site0 = "[1,5,6]"
d6h_cluster_cart_site0 = "[2sqrt(2),sqrt(6),4sqrt(3)]"

# ==================================================
# representative vector in cartesian coordinate (normalized when necessary).
#   it should not be parallel nor perpendicular to site0 vector.
#   d6h with principal axis [111] is equivanlent to oh.
root_cluster_cart_vector0 = "[5,1,6]"
oh_cluster_cart_vector0 = "[5,1,6]"
d6h_cluster_cart_vector0 = "[-2*sqrt(2),sqrt(6),4*sqrt(3)]"

# ==================================================
# scale factor for space group purpose.
#   site0 is scaled as site0/scale.
cluster_scale = 5

# ==================================================
# positions of oh0 cluster in cartesian coordinate (48 points).
#   principal axis [001]. (same as general Wyckoff points)
oh_cluster_cart = "[[x,y,z],[-x,-y,z],[-x,y,-z],[x,-y,-z],[z,x,y],[z,-x,-y],[-z,-x,y],[-z,x,-y],[y,z,x],[-y,z,-x],[y,-z,-x],[-y,-z,x],[y,x,-z],[-y,-x,-z],[y,-x,z],[-y,x,z],[x,z,-y],[-x,z,y],[-x,-z,-y],[x,-z,y],[z,y,-x],[z,-y,x],[-z,y,x],[-z,-y,-x],[-x,-y,-z],[x,y,-z],[x,-y,z],[-x,y,z],[-z,-x,-y],[-z,x,y],[z,x,-y],[z,-x,y],[-y,-z,-x],[y,-z,x],[-y,z,x],[y,z,-x],[-y,-x,z],[y,x,z],[-y,x,-z],[y,-x,-z],[-x,-z,y],[x,-z,-y],[x,z,y],[-x,z,-y],[-z,-y,x],[-z,y,-x],[z,-y,-x],[z,y,x]]"

# ==================================================
# positions of D6h cluster in cartesian coordinate (24 points).
#   principal axis [001]. (same as general Wyckoff points)
d6h_cluster_cart = "[[x,y,z],[-x/2-sqrt(3)y/2,sqrt(3)x/2-y/2,z],[-x/2+sqrt(3)y/2,-sqrt(3)x/2-y/2,z],[-x,-y,z],[x/2+sqrt(3)y/2,-sqrt(3)x/2+y/2,z],[x/2-sqrt(3)y/2,sqrt(3)x/2+y/2,z],[-x/2+sqrt(3)y/2,sqrt(3)x/2+y/2,-z],[x,-y,-z],[-x/2-sqrt(3)y/2,-sqrt(3)x/2+y/2,-z],[x/2-sqrt(3)y/2,-sqrt(3)x/2-y/2,-z],[-x,y,-z],[x/2+sqrt(3)y/2,sqrt(3)x/2-y/2,-z],[-x,-y,-z],[x/2+sqrt(3)y/2,-sqrt(3)x/2+y/2,-z],[x/2-sqrt(3)y/2,sqrt(3)x/2+y/2,-z],[x,y,-z],[-x/2-sqrt(3)y/2,sqrt(3)x/2-y/2,-z],[-x/2+sqrt(3)y/2,-sqrt(3)x/2-y/2,-z],[x/2-sqrt(3)y/2,-sqrt(3)x/2-y/2,z],[-x,y,z],[x/2+sqrt(3)y/2,sqrt(3)x/2-y/2,z],[-x/2+sqrt(3)y/2,sqrt(3)x/2+y/2,z],[x,-y,z],[-x/2-sqrt(3)y/2,-sqrt(3)x/2+y/2,z]]"

# ==================================================
# positions of root cluster in cartesian coordinate (24x4=96 points).
#   contain general points of D6h (d6h0, d6h1, d6h2, d6h3) and Oh (oh0, oh1).
root_cluster_cart = "[[x,y,z],[z,x,y],[y,z,x],[-x/3+2y/3+2z/3,2x/3-y/3+2z/3,2x/3+2y/3-z/3],[2x/3+2y/3-z/3,-x/3+2y/3+2z/3,2x/3-y/3+2z/3],[2x/3-y/3+2z/3,2x/3+2y/3-z/3,-x/3+2y/3+2z/3],[-z,-y,-x],[-y,-x,-z],[-x,-z,-y],[-2x/3-2y/3+z/3,-2x/3+y/3-2z/3,x/3-2y/3-2z/3],[-2x/3+y/3-2z/3,x/3-2y/3-2z/3,-2x/3-2y/3+z/3],[x/3-2y/3-2z/3,-2x/3-2y/3+z/3,-2x/3+y/3-2z/3],[-x,-y,-z],[-z,-x,-y],[-y,-z,-x],[x/3-2y/3-2z/3,-2x/3+y/3-2z/3,-2x/3-2y/3+z/3],[-2x/3-2y/3+z/3,x/3-2y/3-2z/3,-2x/3+y/3-2z/3],[-2x/3+y/3-2z/3,-2x/3-2y/3+z/3,x/3-2y/3-2z/3],[z,y,x],[y,x,z],[x,z,y],[2x/3+2y/3-z/3,2x/3-y/3+2z/3,-x/3+2y/3+2z/3],[2x/3-y/3+2z/3,-x/3+2y/3+2z/3,2x/3+2y/3-z/3],[-x/3+2y/3+2z/3,2x/3+2y/3-z/3,2x/3-y/3+2z/3],[x,-z,-y],[z,-y,-x],[y,-x,-z],[-x/3+2y/3+2z/3,-2x/3-2y/3+z/3,-2x/3+y/3-2z/3],[2x/3+2y/3-z/3,-2x/3+y/3-2z/3,x/3-2y/3-2z/3],[2x/3-y/3+2z/3,x/3-2y/3-2z/3,-2x/3-2y/3+z/3],[-z,x,y],[-y,z,x],[-x,y,z],[-2x/3-2y/3+z/3,-x/3+2y/3+2z/3,2x/3-y/3+2z/3],[-2x/3+y/3-2z/3,2x/3+2y/3-z/3,-x/3+2y/3+2z/3],[x/3-2y/3-2z/3,2x/3-y/3+2z/3,2x/3+2y/3-z/3],[-x,z,y],[-z,y,x],[-y,x,z],[x/3-2y/3-2z/3,2x/3+2y/3-z/3,2x/3-y/3+2z/3],[-2x/3-2y/3+z/3,2x/3-y/3+2z/3,-x/3+2y/3+2z/3],[-2x/3+y/3-2z/3,-x/3+2y/3+2z/3,2x/3+2y/3-z/3],[z,-x,-y],[y,-z,-x],[x,-y,-z],[2x/3+2y/3-z/3,x/3-2y/3-2z/3,-2x/3+y/3-2z/3],[2x/3-y/3+2z/3,-2x/3-2y/3+z/3,x/3-2y/3-2z/3],[-x/3+2y/3+2z/3,-2x/3+y/3-2z/3,-2x/3-2y/3+z/3],[-z,y,-x],[-y,x,-z],[-x,z,-y],[-2x/3-2y/3+z/3,2x/3-y/3+2z/3,x/3-2y/3-2z/3],[-2x/3+y/3-2z/3,-x/3+2y/3+2z/3,-2x/3-2y/3+z/3],[x/3-2y/3-2z/3,2x/3+2y/3-z/3,-2x/3+y/3-2z/3],[x,-y,z],[z,-x,y],[y,-z,x],[-x/3+2y/3+2z/3,-2x/3+y/3-2z/3,2x/3+2y/3-z/3],[2x/3+2y/3-z/3,x/3-2y/3-2z/3,2x/3-y/3+2z/3],[2x/3-y/3+2z/3,-2x/3-2y/3+z/3,-x/3+2y/3+2z/3],[z,-y,x],[y,-x,z],[x,-z,y],[2x/3+2y/3-z/3,-2x/3+y/3-2z/3,-x/3+2y/3+2z/3],[2x/3-y/3+2z/3,x/3-2y/3-2z/3,2x/3+2y/3-z/3],[-x/3+2y/3+2z/3,-2x/3-2y/3+z/3,2x/3-y/3+2z/3],[-x,y,-z],[-z,x,-y],[-y,z,-x],[x/3-2y/3-2z/3,2x/3-y/3+2z/3,-2x/3-2y/3+z/3],[-2x/3-2y/3+z/3,-x/3+2y/3+2z/3,-2x/3+y/3-2z/3],[-2x/3+y/3-2z/3,2x/3+2y/3-z/3,x/3-2y/3-2z/3],[-y,-x,z],[-x,-z,y],[-z,-y,x],[-2x/3+y/3-2z/3,x/3-2y/3-2z/3,2x/3+2y/3-z/3],[x/3-2y/3-2z/3,-2x/3-2y/3+z/3,2x/3-y/3+2z/3],[-2x/3-2y/3+z/3,-2x/3+y/3-2z/3,-x/3+2y/3+2z/3],[y,z,-x],[x,y,-z],[z,x,-y],[2x/3-y/3+2z/3,2x/3+2y/3-z/3,x/3-2y/3-2z/3],[-x/3+2y/3+2z/3,2x/3-y/3+2z/3,-2x/3-2y/3+z/3],[2x/3+2y/3-z/3,-x/3+2y/3+2z/3,-2x/3+y/3-2z/3],[y,x,-z],[x,z,-y],[z,y,-x],[2x/3-y/3+2z/3,-x/3+2y/3+2z/3,-2x/3-2y/3+z/3],[-x/3+2y/3+2z/3,2x/3+2y/3-z/3,-2x/3+y/3-2z/3],[2x/3+2y/3-z/3,2x/3-y/3+2z/3,x/3-2y/3-2z/3],[-y,-z,x],[-x,-y,z],[-z,-x,y],[-2x/3+y/3-2z/3,-2x/3-2y/3+z/3,-x/3+2y/3+2z/3],[x/3-2y/3-2z/3,-2x/3+y/3-2z/3,2x/3+2y/3-z/3],[-2x/3-2y/3+z/3,x/3-2y/3-2z/3,2x/3-y/3+2z/3]]"
