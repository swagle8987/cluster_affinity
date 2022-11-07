import dendropy


def readTree(path,schema):
    tree = dendropy.Tree.get(path=path,schema=schema)
    return tree

