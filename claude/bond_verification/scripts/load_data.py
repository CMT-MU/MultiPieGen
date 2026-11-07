import ast, sys, pickle
src = open(sys.argv[1], encoding="utf-8").read()
tree = ast.parse(src)                       # parse only: the data file is never executed
for node in tree.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "wyckoff_bond_data":
        data = ast.literal_eval(node.value)
        break
pickle.dump(data, open(sys.argv[2], "wb"))
print("groups:", len(data), "entries:", sum(len(v) for v in data.values()))
