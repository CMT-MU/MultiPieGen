# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
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
"""

# ==================================================
# representative site in cartesian coordinate (normalized when necessary), radius of scaled one is given by 1/scale.
#   d6h with principal axis [111] is equivanlent to oh.
oh_cluster_cart_site0 = "[1,5,6]"
d6h_cluster_cart_site0 = "[2sqrt(2),sqrt(6),4sqrt(3)]"

# ==================================================
# representative vector in cartesian coordinate (normalized when necessary).
#   it should not be parallel nor perpendicular to site0 vector.
#   d6h with principal axis [111] is equivanlent to oh.
oh_cluster_cart_vector0 = "[6,1,5]"
d6h_cluster_cart_vector0 = "[-5*sqrt(2)/2,sqrt(6)/2,4*sqrt(3)]"

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
