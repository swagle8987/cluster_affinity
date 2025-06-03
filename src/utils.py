import matplotlib.pyplot as plt
from pathlib import Path

def peek_line(path):
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

def make_matrix_image(matrix,output_path,xlabels=[],ylabels=[],title=""):
    if not xlabels:
        xlabels = range(len(matrix))
    # we assume it is a square matrix
    if not ylabels:
        ylabels = range(len(matrix[0]))
    fig,ax = plt.subplots()
    im = ax.imshow(matrix)
    ax.set_xticks(range(len(xlabels)),labels=xlabels)
    ax.set_yticks(range(len(ylabels)),labels=ylabels)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            text = ax.text(j,i,round(matrix[i][j],2), ha="center", va="center", color="w")
    ax.set_title(title)
    fig.tight_layout()
    plt.savefig(output_path)

def check_input_trees(tlist):
    listtaxa = len(tlist[0])
    for i in tlist:
        if listtaxa != len(i):
            raise RuntimeWarning("The input trees do not have the same taxa")
    return True
