from ete4.smartview.explorer import add_tree, explore
from ete4.smartview import (
    CircleFace,
    TextFace,
    Layout,
    PropFace,
    BASIC_LAYOUT,
    HeatmapFace,
    LegendFace,
)
from .config import cmap,COLOR_RANGE
from matplotlib.colors import rgb2hex, to_rgb


def generate_layout(cost, color_only):

    def draw_tree(tree):
        yield LegendFace(
            "Range of the cost:\n 0 indicates a perfect match while 1 indicates the maximum cost.\n The relative width of the edges is also an indicator of cost",
            "continuous",
            value_range=(0, 1),
            color_range=COLOR_RANGE[::-1],
            position="aligned",
            anchor=(0, 0),
        )
        yield {"node-height-min": 0, "collapsed": {"shape": "outline"}}

    def draw_node(node):
        if not color_only:
            if node.is_leaf or node.is_root:
                yield PropFace("name", fs_min=11, position="right")
            else:
                yield PropFace("c_dist", fmt="%.3f", fs_min=11, position="right")
        if "c_dist" in node.props:
            color = rgb2hex(cmap(node.props["c_dist"])[:3])
        else:
            color = rgb2hex(cmap(0)[:3])  ## The root and leaves have zero cost always
        yield {
            "hz-line": {
                "stroke": color,
                "stroke-width": max(
                    1, node.props[cost] * 5 if cost in node.props else 1
                ),
            },
            "vt-line": {
                "stroke": color,
                "stroke-width": max(
                    1, node.props[cost] * 5 if cost in node.props else 1
                ),
            },
        }

    return Layout(name=cost, draw_node=draw_node, draw_tree=draw_tree)


def start_web_server(t1, t2, cost, color_only, t1_name, t2_name):
    layout = generate_layout(cost, color_only)
    add_tree(
        t1,
        name="Source tree ({})".format(t1_name),
        layouts=[layout],
    )
    explore(
        t2,
        name="Target tree ({})".format(t2_name),
        layouts=[BASIC_LAYOUT],
    )
    print("Press any key to stop the server and finish")
    input()
