"""
Correspondence between zj and atomic variable.
- correspondence for each bond cluster, dict[bond_name, dict[zj, expression in terms of atomic variables] ].
- only for SAMB with identity irrep.
- atomic variable of (m,n) component at bond 1 is given by g(m,n) +i h(m,n).
"""

N198L_var = {
    "Co;Co_001_1": {
        "z5": "2*sqrt(6)*g_{s_Co2,s_Co1}",
        "z6": "5*sqrt(6)*g_{px_Co2,s_Co1}/4+3*sqrt(2)*g_{py_Co2,s_Co1}/4+sqrt(6)*g_{pz_Co2,s_Co1}/2+5*sqrt(6)*g_{s_Co2,px_Co1}/4+3*sqrt(2)*g_{s_Co2,py_Co1}/4+sqrt(6)*g_{s_Co2,pz_Co1}/2",
        "z7": "5*sqrt(6)*g_{px_Co2,s_Co1}/4+3*sqrt(2)*g_{py_Co2,s_Co1}/4+sqrt(6)*g_{pz_Co2,s_Co1}/2-5*sqrt(6)*g_{s_Co2,px_Co1}/4-3*sqrt(2)*g_{s_Co2,py_Co1}/4-sqrt(6)*g_{s_Co2,pz_Co1}/2",
        "z8": "5*sqrt(2463)*g_{px_Co2,s_Co1}/821+63*sqrt(821)*g_{py_Co2,s_Co1}/821-44*sqrt(2463)*g_{pz_Co2,s_Co1}/821-5*sqrt(2463)*g_{s_Co2,px_Co1}/821-63*sqrt(821)*g_{s_Co2,py_Co1}/821+44*sqrt(2463)*g_{s_Co2,pz_Co1}/821",
        "z9": "-129*sqrt(1642)*g_{px_Co2,s_Co1}/3284+115*sqrt(4926)*g_{py_Co2,s_Co1}/3284+75*sqrt(1642)*g_{pz_Co2,s_Co1}/1642-129*sqrt(1642)*g_{s_Co2,px_Co1}/3284+115*sqrt(4926)*g_{s_Co2,py_Co1}/3284+75*sqrt(1642)*g_{s_Co2,pz_Co1}/1642",
        "z10": "5*sqrt(2463)*g_{px_Co2,s_Co1}/821+63*sqrt(821)*g_{py_Co2,s_Co1}/821-44*sqrt(2463)*g_{pz_Co2,s_Co1}/821+5*sqrt(2463)*g_{s_Co2,px_Co1}/821+63*sqrt(821)*g_{s_Co2,py_Co1}/821-44*sqrt(2463)*g_{s_Co2,pz_Co1}/821",
        "z11": "-129*sqrt(1642)*g_{px_Co2,s_Co1}/3284+115*sqrt(4926)*g_{py_Co2,s_Co1}/3284+75*sqrt(1642)*g_{pz_Co2,s_Co1}/1642+129*sqrt(1642)*g_{s_Co2,px_Co1}/3284-115*sqrt(4926)*g_{s_Co2,py_Co1}/3284-75*sqrt(1642)*g_{s_Co2,pz_Co1}/1642",
        "z12": "2*sqrt(2)*(g_{px_Co2,px_Co1}+g_{py_Co2,py_Co1}+g_{pz_Co2,pz_Co1})",
        "z13": "86*sqrt(4630)*g_{px_Co2,px_Co1}/2315+15*sqrt(24630)*g_{px_Co2,py_Co1}/1642+69*sqrt(8210)*g_{px_Co2,pz_Co1}/3284+15*sqrt(24630)*g_{py_Co2,px_Co1}/1642-46*sqrt(4630)*g_{py_Co2,py_Co1}/2315-129*sqrt(24630)*g_{py_Co2,pz_Co1}/16420+69*sqrt(8210)*g_{pz_Co2,px_Co1}/3284-129*sqrt(24630)*g_{pz_Co2,py_Co1}/16420-8*sqrt(4630)*g_{pz_Co2,pz_Co1}/463",
        "z14": "75*sqrt(1642)*g_{px_Co2,py_Co1}/1642-115*sqrt(4926)*g_{px_Co2,pz_Co1}/3284-75*sqrt(1642)*g_{py_Co2,px_Co1}/1642-129*sqrt(1642)*g_{py_Co2,pz_Co1}/3284+115*sqrt(4926)*g_{pz_Co2,px_Co1}/3284+129*sqrt(1642)*g_{pz_Co2,py_Co1}/3284",
        "z15": "sqrt(6)*g_{px_Co2,py_Co1}/2+3*sqrt(2)*g_{px_Co2,pz_Co1}/4+sqrt(6)*g_{py_Co2,px_Co1}/2+5*sqrt(6)*g_{py_Co2,pz_Co1}/4+3*sqrt(2)*g_{pz_Co2,px_Co1}/4+5*sqrt(6)*g_{pz_Co2,py_Co1}/4",
        "z16": "44*sqrt(2463)*g_{px_Co2,py_Co1}/821-63*sqrt(821)*g_{px_Co2,pz_Co1}/821+44*sqrt(2463)*g_{py_Co2,px_Co1}/821-5*sqrt(2463)*g_{py_Co2,pz_Co1}/821-63*sqrt(821)*g_{pz_Co2,px_Co1}/821-5*sqrt(2463)*g_{pz_Co2,py_Co1}/821",
        "z17": "86*sqrt(6945)*g_{px_Co2,px_Co1}/2315-15*sqrt(4105)*g_{px_Co2,py_Co1}/821-23*sqrt(12315)*g_{px_Co2,pz_Co1}/1642-15*sqrt(4105)*g_{py_Co2,px_Co1}/821-46*sqrt(6945)*g_{py_Co2,py_Co1}/2315+129*sqrt(4105)*g_{py_Co2,pz_Co1}/8210-23*sqrt(12315)*g_{pz_Co2,px_Co1}/1642+129*sqrt(4105)*g_{pz_Co2,py_Co1}/8210-8*sqrt(6945)*g_{pz_Co2,pz_Co1}/463",
        "z18": "sqrt(6)*g_{px_Co2,py_Co1}/2-3*sqrt(2)*g_{px_Co2,pz_Co1}/4-sqrt(6)*g_{py_Co2,px_Co1}/2+5*sqrt(6)*g_{py_Co2,pz_Co1}/4+3*sqrt(2)*g_{pz_Co2,px_Co1}/4-5*sqrt(6)*g_{pz_Co2,py_Co1}/4",
        "z19": "-44*sqrt(2463)*g_{px_Co2,py_Co1}/821-63*sqrt(821)*g_{px_Co2,pz_Co1}/821+44*sqrt(2463)*g_{py_Co2,px_Co1}/821+5*sqrt(2463)*g_{py_Co2,pz_Co1}/821+63*sqrt(821)*g_{pz_Co2,px_Co1}/821-5*sqrt(2463)*g_{pz_Co2,py_Co1}/821",
        "z20": "2*sqrt(1389)*(-g_{px_Co2,px_Co1}-21*g_{py_Co2,py_Co1}+22*g_{pz_Co2,pz_Co1})/463",
    }
}
