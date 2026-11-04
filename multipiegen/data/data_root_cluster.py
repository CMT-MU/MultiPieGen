"""
Definition of root cluseter.

- root_cluster_cart_site0 (str): representative site of root cluster.
- root_cluster_cart_vector0 (str): representative vector of root cluster.
- root_cluster_scale (int): scale factor (site -> site/scale).
- root_cluster_cart (str): sites of root cluster (first=[x,y,z]).
- oh_cluster_frac (str): sites of Oh cluster (first=[x,y,z]).
- d6h_cluster_frac (str): sites of D6h cluster (first=[x,y,z]).

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
# representative site in cartesian coordinate (to be normalized as r=1), radius of scaled one is given by 1/scale.
root_cluster_cart_site0 = "[5,sqrt(3),2]"

# ==================================================
# representative vector in cartesian coordinate (to be normalized as v=1).
#   it should not be parallel nor perpendicular to site0 vector.
root_cluster_cart_vector0 = "[sqrt(3),5,2]"

# ==================================================
# scale factor for space group purpose.
#   site0 is scaled as site0/scale.
root_cluster_scale = 5

# ==================================================
# positions of root cluster in cartesian coordinate (24x4=96 points).
#   contain general points of D6h (d6h0, d6h1, d6h2, d6h3) and Oh (oh0, oh1).
root_cluster_cart = "[[x,y,z],[z,x,y],[y,z,x],[-x/3+2y/3+2z/3,2x/3-y/3+2z/3,2x/3+2y/3-z/3],[2x/3+2y/3-z/3,-x/3+2y/3+2z/3,2x/3-y/3+2z/3],[2x/3-y/3+2z/3,2x/3+2y/3-z/3,-x/3+2y/3+2z/3],[-z,-y,-x],[-y,-x,-z],[-x,-z,-y],[-2x/3-2y/3+z/3,-2x/3+y/3-2z/3,x/3-2y/3-2z/3],[-2x/3+y/3-2z/3,x/3-2y/3-2z/3,-2x/3-2y/3+z/3],[x/3-2y/3-2z/3,-2x/3-2y/3+z/3,-2x/3+y/3-2z/3],[-x,-y,-z],[-z,-x,-y],[-y,-z,-x],[x/3-2y/3-2z/3,-2x/3+y/3-2z/3,-2x/3-2y/3+z/3],[-2x/3-2y/3+z/3,x/3-2y/3-2z/3,-2x/3+y/3-2z/3],[-2x/3+y/3-2z/3,-2x/3-2y/3+z/3,x/3-2y/3-2z/3],[z,y,x],[y,x,z],[x,z,y],[2x/3+2y/3-z/3,2x/3-y/3+2z/3,-x/3+2y/3+2z/3],[2x/3-y/3+2z/3,-x/3+2y/3+2z/3,2x/3+2y/3-z/3],[-x/3+2y/3+2z/3,2x/3+2y/3-z/3,2x/3-y/3+2z/3],[x,-z,-y],[z,-y,-x],[y,-x,-z],[-x/3+2y/3+2z/3,-2x/3-2y/3+z/3,-2x/3+y/3-2z/3],[2x/3+2y/3-z/3,-2x/3+y/3-2z/3,x/3-2y/3-2z/3],[2x/3-y/3+2z/3,x/3-2y/3-2z/3,-2x/3-2y/3+z/3],[-z,x,y],[-y,z,x],[-x,y,z],[-2x/3-2y/3+z/3,-x/3+2y/3+2z/3,2x/3-y/3+2z/3],[-2x/3+y/3-2z/3,2x/3+2y/3-z/3,-x/3+2y/3+2z/3],[x/3-2y/3-2z/3,2x/3-y/3+2z/3,2x/3+2y/3-z/3],[-x,z,y],[-z,y,x],[-y,x,z],[x/3-2y/3-2z/3,2x/3+2y/3-z/3,2x/3-y/3+2z/3],[-2x/3-2y/3+z/3,2x/3-y/3+2z/3,-x/3+2y/3+2z/3],[-2x/3+y/3-2z/3,-x/3+2y/3+2z/3,2x/3+2y/3-z/3],[z,-x,-y],[y,-z,-x],[x,-y,-z],[2x/3+2y/3-z/3,x/3-2y/3-2z/3,-2x/3+y/3-2z/3],[2x/3-y/3+2z/3,-2x/3-2y/3+z/3,x/3-2y/3-2z/3],[-x/3+2y/3+2z/3,-2x/3+y/3-2z/3,-2x/3-2y/3+z/3],[-z,y,-x],[-y,x,-z],[-x,z,-y],[-2x/3-2y/3+z/3,2x/3-y/3+2z/3,x/3-2y/3-2z/3],[-2x/3+y/3-2z/3,-x/3+2y/3+2z/3,-2x/3-2y/3+z/3],[x/3-2y/3-2z/3,2x/3+2y/3-z/3,-2x/3+y/3-2z/3],[x,-y,z],[z,-x,y],[y,-z,x],[-x/3+2y/3+2z/3,-2x/3+y/3-2z/3,2x/3+2y/3-z/3],[2x/3+2y/3-z/3,x/3-2y/3-2z/3,2x/3-y/3+2z/3],[2x/3-y/3+2z/3,-2x/3-2y/3+z/3,-x/3+2y/3+2z/3],[z,-y,x],[y,-x,z],[x,-z,y],[2x/3+2y/3-z/3,-2x/3+y/3-2z/3,-x/3+2y/3+2z/3],[2x/3-y/3+2z/3,x/3-2y/3-2z/3,2x/3+2y/3-z/3],[-x/3+2y/3+2z/3,-2x/3-2y/3+z/3,2x/3-y/3+2z/3],[-x,y,-z],[-z,x,-y],[-y,z,-x],[x/3-2y/3-2z/3,2x/3-y/3+2z/3,-2x/3-2y/3+z/3],[-2x/3-2y/3+z/3,-x/3+2y/3+2z/3,-2x/3+y/3-2z/3],[-2x/3+y/3-2z/3,2x/3+2y/3-z/3,x/3-2y/3-2z/3],[-y,-x,z],[-x,-z,y],[-z,-y,x],[-2x/3+y/3-2z/3,x/3-2y/3-2z/3,2x/3+2y/3-z/3],[x/3-2y/3-2z/3,-2x/3-2y/3+z/3,2x/3-y/3+2z/3],[-2x/3-2y/3+z/3,-2x/3+y/3-2z/3,-x/3+2y/3+2z/3],[y,z,-x],[x,y,-z],[z,x,-y],[2x/3-y/3+2z/3,2x/3+2y/3-z/3,x/3-2y/3-2z/3],[-x/3+2y/3+2z/3,2x/3-y/3+2z/3,-2x/3-2y/3+z/3],[2x/3+2y/3-z/3,-x/3+2y/3+2z/3,-2x/3+y/3-2z/3],[y,x,-z],[x,z,-y],[z,y,-x],[2x/3-y/3+2z/3,-x/3+2y/3+2z/3,-2x/3-2y/3+z/3],[-x/3+2y/3+2z/3,2x/3+2y/3-z/3,-2x/3+y/3-2z/3],[2x/3+2y/3-z/3,2x/3-y/3+2z/3,x/3-2y/3-2z/3],[-y,-z,x],[-x,-y,z],[-z,-x,y],[-2x/3+y/3-2z/3,-2x/3-2y/3+z/3,-x/3+2y/3+2z/3],[x/3-2y/3-2z/3,-2x/3+y/3-2z/3,2x/3+2y/3-z/3],[-2x/3-2y/3+z/3,x/3-2y/3-2z/3,2x/3-y/3+2z/3]]"

# ==================================================
# positions of oh0 cluster in fractional coordinate (48 points).
#   principal axis [001].
oh_cluster_frac = "[[x,y,z],[-x,-y,z],[-x,y,-z],[x,-y,-z],[z,x,y],[z,-x,-y],[-z,-x,y],[-z,x,-y],[y,z,x],[-y,z,-x],[y,-z,-x],[-y,-z,x],[y,x,-z],[-y,-x,-z],[y,-x,z],[-y,x,z],[x,z,-y],[-x,z,y],[-x,-z,-y],[x,-z,y],[z,y,-x],[z,-y,x],[-z,y,x],[-z,-y,-x],[-x,-y,-z],[x,y,-z],[x,-y,z],[-x,y,z],[-z,-x,-y],[-z,x,y],[z,x,-y],[z,-x,y],[-y,-z,-x],[y,-z,x],[-y,z,x],[y,z,-x],[-y,-x,z],[y,x,z],[-y,x,-z],[y,-x,-z],[-x,-z,y],[x,-z,-y],[x,z,y],[-x,z,-y],[-z,-y,x],[-z,y,-x],[z,-y,-x],[z,y,x]]"

# ==================================================
# positions of D6h cluster in fractional coordinate (24 points).
#   principal axis [001].
#   these points coincide with those in d6h0 cluster, when principal axis is taken as [111].
d6h_cluster_frac = "[[-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,-sqrt(3)x/3-sqrt(3)y/3-sqrt(3)z/3],[sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[-sqrt(2)x/3-sqrt(2)y/3+2sqrt(2)z/3,-2sqrt(2)x/3+sqrt(2)y/3+sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[-sqrt(2)x/3+2sqrt(2)y/3-sqrt(2)z/3,sqrt(2)x/3+sqrt(2)y/3-2sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3],[2sqrt(2)x/3-sqrt(2)y/3-sqrt(2)z/3,sqrt(2)x/3-2sqrt(2)y/3+sqrt(2)z/3,sqrt(3)x/3+sqrt(3)y/3+sqrt(3)z/3]]"
