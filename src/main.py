import argparse

from cluster_computation import rooted_cluster_affinity,calculate_rooted_tau,rooted_cluster_support,calculate_rooted_phi
from ete4 import Tree,nexus
from ete4.smartview.explorer import add_tree,explore

from ete4.smartview import CircleFace, TextFace,Layout,PropFace,BASIC_LAYOUT,HeatmapFace,LegendFace
from utils import peek_line,convert_dict_to_2d_array, make_matrix_image, check_input_trees

import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap,rgb2hex,to_rgb


import os


COLOR_RANGE=[rgb2hex(to_rgb("xkcd:chartreuse")),rgb2hex(to_rgb("xkcd:orange red"))]
cmap = LinearSegmentedColormap.from_list("custom_tree_map",colors=COLOR_RANGE)

def draw_tree(tree):
    yield LegendFace("Range of the cost:\n 0 indicates a perfect match while 1 indicates the maximum cost.\n The relative width of the edges is also an indicator of cost", "continuous",value_range=(0,1),color_range=COLOR_RANGE[::-1],position="aligned",anchor=(0,0))
    yield { 'node-height-min':0,
            'collapsed':{"shape":"outline"}}

def draw_node(node):
    if node.is_leaf:
        yield PropFace("name",fs_min=11,position="right")
    else:
        yield PropFace("c_dist",fmt="%.3f",fs_min=11,position="right")
    color = rgb2hex(cmap(node.props["c_dist"])[:3])
    yield {
            "hz-line": {
                    "stroke":color,
                    "stroke-width":max(1,node.props["c_dist"]*5)
                },
            "vt-line":{
                "stroke":color,
                "stroke-width":max(1,node.props["c_dist"]*5)
                }
            }

def run_script(cost,args):

    ftype = args.filetype
    if not ftype:
        ftype = "nexus" if peek_line(args.t1) == "#NEXUS" else "newick"
    
    if ftype=="newick":
        t1 = Tree(open(args.t1),parser=1)
        t2 = Tree(open(args.t2),parser=1)
    else:
        t1 = nexus.get_trees(open(args.t1).read())["tree_1"]
        t2 = nexus.get_trees(open(args.t2).read())["tree_1"]

    if check_input_trees([t1,t2]):
        if cost=="cluster_affinity":
            dist = rooted_cluster_affinity(t1,t2)/calculate_rooted_tau(t1)
            name="Cluster Affinity Layout"
        else:
            dist = rooted_cluster_support(t1,t2)/calculate_rooted_phi(t1)
            name="Cluster Support Layout"
        if not args.cli:
            layout=Layout(name=name,draw_tree=draw_tree,draw_node=draw_node)
            add_tree(t1,name="Source tree",layouts=[layout])
            explore(t2,name="Target tree",layouts=[BASIC_LAYOUT])
            print("Press any key to stop the server and finish")
            input()
        else:
            print(dist)



def cluster_affinity_script():

    parser = argparse.ArgumentParser(
            prog='Cluster Affinity',
            description='Calculates the Asymmetric Cluster Affinity cost from t1 to t2',
    )

    parser.add_argument('t1', help='The source tree from which the cost is to be calculated')
    parser.add_argument('t2', help='The target tree to which the cost is to be calculated')
    parser.add_argument('-t','--filetype', help="The input file format")
    parser.add_argument('--cli',help="Disables interactive browser", action='store_true')


    args = parser.parse_args()

    run_script("cluster_affinity",args)


def cluster_support_script():

    parser = argparse.ArgumentParser(
            prog='Cluster Support',
            description='Calculates the Asymmetric Cluster Support cost from t1 to t2',
    )

    parser.add_argument('t1', help='The source tree from which the cost is to be calculated')
    parser.add_argument('t2', help='The target tree to which the cost is to be calculated')
    parser.add_argument('-t','--filetype', help="The input file format")
    parser.add_argument('--cli',help="Disables interactive browser", action='store_true')

    args = parser.parse_args()

    run_script("cluster_affinity",args)

if __name__=="__main__":
    cluster_affinity_script()
