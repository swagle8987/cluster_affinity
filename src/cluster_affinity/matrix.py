import argparse
from .cluster_computation import (
    calculate_rooted_tau,
    compute_transfer_index,
    calculate_rooted_phi,
)
from ete4 import Tree, nexus


from .reader import FileFormatError
from pathlib import Path
import numpy as np
from .cli import get_dist

import os

from .config import COLOR_RANGE,cmap
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt

def get_matrix_args():
    parser = argparse.ArgumentParser(
        prog="Cluster Matrix",
        description="Computes the cluster matrix for pairwise costs",
    )
    parser.add_argument("t", help="The treefile containing the multiple trees",nargs="+",type=Path)
    parser.add_argument("outfile", help="The output image path",type=Path)
    parser.add_argument("--title", help="The title for the matrix",default="")
    parser.add_argument("--cost", help="The cost to use for the pairwise comparison")
    parser.add_argument("--filetype", help="The input filetype")
    parser.add_argument("--unrooted",help="Consider all trees to be unrooted",action="store_true")
    parser.add_argument(
        "--average",
        help="Adds a column displaying average cost to the matrix",
        action="store_true",
    )
    parser.add_argument("--csv_output",type=Path)
    parser.add_argument("--autoscale",help="Changes color scaling from 0-1 to min to max",action="store_true")
    return parser

def _get_trees(paths,ftype=None):
    trees = dict()
    for path in paths:
        if isinstance(path,Path):
            with open(path,"r") as file:
                data = file.read().strip()
        elif isinstance(path, str):
            data= path
        else:
            raise RuntimeError(f"Invalid input {path} to _get_tree")

        if not ftype:
            ftype = "nexus" if data[:6] == "#NEXUS" else "newick"

        if ftype == "newick":
            for ind, t in enumerate(data.split("\n")):
                print(f"Reading tree {ind}")
                trees[f"tree_{ind+1}"] = Tree(t.strip(), parser=1)
        elif ftype == "nexus":
            trees.update(nexus.loads(data))
        else:
            raise RuntimeWarning("Unsupported file format")
    return trees

def cluster_matrix():
    parser = get_matrix_args()
    args = parser.parse_args()
    if args.cost =="cluster_support":
        cost = "Cluster Support"
    else:
        cost = "Cluster Affinity"
    trees = _get_trees(args.t,args.filetype)
    if args.average:
        matrix = np.zeros([len(trees)+1,len(trees)+1])
    else:
        matrix = np.zeros([len(trees),len(trees)])
    for i,t1 in enumerate(trees.values()):
        for j,t2 in enumerate(trees.values()):
            matrix[i][j] = get_dist(cost,t1,t2,(not args.unrooted),(not args.unrooted))
    if args.average:
        xlabels = list(trees.keys()) + ["average"]
        matrix[-1] = np.average(matrix,axis=0)
        matrix[:,-1] = np.average(matrix,axis=1)
    else:
        xlabels = list(trees.keys())
    if args.autoscale:
        vmin = np.min(matrix)
        vmax = np.max(matrix)
    else:
        vmin=0
        vmax=1
    make_matrix_image(
        matrix, args.outfile, xlabels=xlabels, ylabels=xlabels, cmap=cmap, vmin=0,vmax=1
    )

def make_matrix_image(matrix, output_path, xlabels=[], ylabels=[], title="", cmap=None,vmin=0,vmax=1):
    if not xlabels:
        xlabels = range(len(matrix))
    # we assume it is a square matrix
    if not ylabels:
        ylabels = range(len(matrix[0]))
    fig, ax = plt.subplots()
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlabel("Target tree")
    ax.set_ylabel("Source tree")
    ax.set_xticks(range(len(xlabels)), labels=xlabels)
    ax.set_yticks(range(len(ylabels)), labels=ylabels)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            text = ax.text(
                j, i, round(matrix[i][j], 3), ha="center", va="center", color="w"
            )
    ax.set_title(title)
    fig.set_size_inches(len(matrix), len(matrix[0]))
    fig.tight_layout()
    plt.savefig(output_path, dpi=100)
