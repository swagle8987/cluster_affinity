#!/usr/bin/env python3

from ete4 import Tree, nexus
import os
import re
from pathlib import Path
import warnings


class FileFormatError(Exception):
    def __init__(self, path, error):
        self.path = path
        self.error = error
        super().__init__(self.error)

    def __str__(self):
        return "The following error occured when handling file {}: \n {}".format(
            self.path, self.error
        )


def get_tree(path, ftype=None):
    rooting = None
    tree = None
    if isinstance(path, Path):
        with open(path, "r") as file:
            data = file.read().strip()
    elif isinstance(path, str):
        data = path
    else:
        raise RuntimeError("Invalid input {} to get_tree".format(path))
    if not ftype:
        ftype = "nexus" if data[:6] == "#NEXUS" else "newick"
    if ftype == "newick":
        string_match = re.fullmatch(
            r"(\[&[RU]\])? ?([(\w\d:\.\-,\/)]+;)", data, flags=re.MULTILINE
        )
        if string_match:
            if string_match[1]:
                if string_match[1] == "[&U]":
                    rooting = False
                elif string_match[1] == "[&R]":
                    rooting = True
                else:
                    raise FileFormatError(path, "Invalid rooting state in newick file")
            else:
                rooting = True
            tree = Tree(string_match[2], parser=1)
        else:
            raise FileFormatError(
                path,
                "Could not match regex to data found in newick file;\n Is the data malformed?",
            )
    elif ftype == "nexus":
        tdict = nexus.loads(data)
        tree = list(tdict.values())[0]
        tree.add_prop("tree_name", list(tdict.keys())[0])
        if "[&U]" in data:
            rooting = False
        else:
            rooting = True
        
    else:
        raise FileFormatError(path, "Invalid file type")
    return tree, rooting
