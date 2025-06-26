import matplotlib.pyplot as plt
from pathlib import Path
import ete4

def peek_line(path:'str')->'str':
    with open(path,"r") as file:
        pos=file.tell()
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

def make_matrix_image(matrix,output_path,xlabels=[],ylabels=[],title="",cmap=None):
    if not xlabels:
        xlabels = range(len(matrix))
    # we assume it is a square matrix
    if not ylabels:
        ylabels = range(len(matrix[0]))
    fig,ax = plt.subplots()
    im = ax.imshow(matrix,cmap,vmin=0,vmax=1)
    ax.set_xlabel("Target tree")
    ax.set_ylabel("Source tree")
    ax.set_xticks(range(len(xlabels)),labels=xlabels)
    ax.set_yticks(range(len(ylabels)),labels=ylabels)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            text = ax.text(j,i,round(matrix[i][j],3), ha="center", va="center", color="w")
    ax.set_title(title)
    fig.set_size_inches(len(matrix),len(matrix[0]))
    fig.tight_layout()
    plt.savefig(output_path,dpi=100)

def check_input_trees(tlist:list[ete4.Tree])->bool:
    listtaxa = set(tlist[0].leaf_names())
    for i in tlist:
        for l in i.leaf_names():
            if l not in listtaxa:
                raise RuntimeError("The input trees do not have the same taxa")
    return True
