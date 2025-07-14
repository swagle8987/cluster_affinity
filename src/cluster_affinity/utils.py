import matplotlib.pyplot as plt
from pathlib import Path
import ete4


def peek_line(path: "str") -> "str":
    with open(path, "r") as file:
        pos = file.tell()
        line = file.readline()
        file.seek(pos)
    return line.strip()


def convert_dict_to_2d_array(d):
    arr = []
    for i in d.values():
        l = []
        for j in i.values():
            l.append(j)
        arr.append(l)
    return arr




def check_input_trees(tlist: list[ete4.Tree]) -> bool:
    listtaxa = set(tlist[0].leaf_names())
    for i in tlist:
        for l in i.leaf_names():
            if l not in listtaxa:
                return False
    return True
