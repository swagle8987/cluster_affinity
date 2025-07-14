## Helper file for configuration

import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, rgb2hex, to_rgb

COLOR_RANGE = [
    rgb2hex(to_rgb("xkcd:evergreen")),
    rgb2hex("#ffad01"),
    rgb2hex("#ff0000"),
]
cmap = LinearSegmentedColormap.from_list("custom_tree_map", colors=COLOR_RANGE)

