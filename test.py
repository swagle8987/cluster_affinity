from src.extendedtree import ExtendedTree
from src.affinity import tree_affinity

if __name__ == "__main__":
    t1 = ExtendedTree.get(path="sample1.tre",schema="newick")
    t2 = ExtendedTree.get(path="sample2.tre",schema="newick")
    tree_affinity(t1,t2)
    string = ""
    string += t1.convertTreeToString()
    string += t2.convertTreeToString()
    print(string)
