#!/usr/bin/env python3
#
#
import importlib.metadata
import argparse

from reader import get_tree
from cluster_computation import (
    compute_transfer_index,
    calculate_rooted_phi,
    calculate_rooted_tau,
    calculate_unrooted_phi,
    calculate_unrooted_tau,
)

from visualization import start_web_server
import pathlib


def check_input_trees(tlist) -> bool:
    listtaxa = set(tlist[0].leaf_names())
    for i in tlist:
        for l in i.leaf_names():
            if l not in listtaxa:
                return False
    return True


def get_default_args(name, description):
    parser = argparse.ArgumentParser(
        prog=name,
        description=description,
    )

    parser.add_argument(
        "t1",
        help="The source tree from which the cost is to be calculated",
        type=pathlib.Path,
    )
    parser.add_argument(
        "t2",
        help="The target tree to which the cost is to be calculated",
        type=pathlib.Path,
    )
    parser.add_argument("-t", "--filetype", help="The input file format")
    parser.add_argument(
        "--cli", help="Disables interactive browser", action="store_true"
    )
    parser.add_argument(
        "--color_only",
        help="Disables node annotations(useful when visualizing dense trees)",
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Version number",
        action="version",
        version=importlib.metadata.version("cluster_affinity"),
    )
    return parser


def cluster_affinity():
    cluster_cost_script(support=False)


def cluster_support():
    cluster_cost_script(support=True)


def cluster_cost_script(support=False):
    if support:
        cost = "Cluster Support"
    else:
        cost = "Cluster Affinity"
    parser = get_default_args(
        name=cost,
        description="Calculates the Asymmetric {} cost from t1 to t2".format(cost),
    )
    args = parser.parse_args()
    t1, t1_is_rooted = get_tree(args.t1, args.filetype)
    t2, t2_is_rooted = get_tree(args.t2, args.filetype)
    if check_input_trees([t1, t2]):
        if t1_is_rooted:
            cost = "Rooted " + cost
            dist = compute_transfer_index(
                t1, t2, cost=cost, annotate_cost=(not args.cli)
            ) / calculate_rooted_tau(t1)
        else:
            cost = "Unrooted " + cost
            dist = compute_transfer_index(
                t1, t2, cost=cost, annotate_cost=(not args.cli)
            ) / calculate_unrooted_tau(t1)
        if not args.cli:
            print(t1.id)

            start_web_server(t1, t2, cost, args.color_only, args.t1.name, args.t2.name)
        else:
            print(dist)
    else:
        raise RuntimeError(
            "Tree {} and {} have different taxa set".format(args.t1, args.t2)
        )
