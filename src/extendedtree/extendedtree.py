import colorsys
from dendropy import Tree,Node
from .counter import Counter
import warnings

class ExtendedTree(Tree):

    """Extended Tree class for dendropy """

    def __init__(self,*args,**kwargs):
        """TODO: to be defined. """
        Tree.__init__(self,*args,**kwargs)
        self.clusters = {}
        self.cluster_lookup = {}
        self.node_lookup = {}
        self.file_path = ""
        self.enrich_tree()

    def reset_tree(self):
        self.clusters = {}
        self.cluster_lookup  = {}
        self.node_lookup = {}
        for node in self.postorder_node_iter():
            node.label = None


    def get(**kwargs):
        newTree = Tree.get(**kwargs)
        newTree = ExtendedTree(newTree)
        return newTree
    
    def setFilePath(self,filepath):
       self.file_path = filepath


    def color_nodes(self,criteria):
        nodes = self.find_nodes(criteria)
        n = len(nodes)
        hsv_colors = [(i / n, 0.7, 0.7) for i in range(n)]
        rgb_colors = [colorsys.hsv_to_rgb(*hsv) for hsv in hsv_colors]
        rgb_colors = [[round(255 * color[i]) for i in range(len(color))] for color in rgb_colors]
        hex_colors = ['#%02x%02x%02x' % (color[0], color[1], color[2]) for color in rgb_colors]
        for i,color in zip(nodes,hex_colors):
            i.annotations.add_new("!color",color)

    # SPR stuff
    def spr(self,src_node,dest_node):
        if src_node == dest_node:
            raise ValueError("The source node is the dest_node")
        if src_node.parent_node == None:
            raise ValueError("The source node is the root node")
  #      elif self.cluster_lookup[dest_node.label] < self.cluster_lookup[src_node.label]:
  #          raise ValueError("The dest_node is a descendant of the src_node")
        else:
            if src_node.parent_node == dest_node:
                warnings.warn("The destination node is the parent of the source node")
            # If dest_node is full then add new node
            if len(dest_node.child_nodes()) >= 2:
                n = Node()
                if dest_node.parent_node:
                    p = dest_node.parent_node
                    p.remove_child(dest_node)
                    p.add_child(n)
                else:
                    self.seed_node = n
                n.add_child(dest_node)
                dest_node = n
            l = src_node.parent_node
            l.remove_child(src_node)
            dest_node.add_child(src_node)
            self.encode_bipartitions(True, False)
            self.reset_tree()
            self.enrich_tree()


    def enrich_tree(self):
        c = Counter()
        for node in self.postorder_node_iter():
            cluster = set()
            if node.is_leaf():
                cluster.add(str(node.taxon))
            else:
                for child in node.child_node_iter():
                    cluster = cluster | self.cluster_lookup[child.label]
            if not node.label:
                node.label = str(c.get_count())
                node.annotations.add_new("node_number",c.get_count())
                c.count_up()
            self.cluster_lookup[node.label] = frozenset(cluster)
            self.node_lookup[node.label] = node
    
    def generate_clusters(self):
        for node in self.postorder_node_iter():
            cluster = set()
            if node.is_leaf():
                cluster.add(str(node.taxon))
                self.cluster_lookup[node.label] = frozenset(cluster)
            else:
                for child in node.child_node_iter():
                    cluster = cluster | self.cluster_lookup[child.label]
                self.cluster_lookup[node.label] = frozenset(cluster)
        self.clusters = list(self.cluster_lookup.values())

    def get_cluster(self,label):
        return self.cluster_lookup[label]
        
    def writeTreeToFile(self,path,schema="newick"):
        self.write(path=path,schema=schema,suppress_internal_node_labels=False,suppress_annotations=False,suppress_rooting=True)
        
    def convertTreeToString(self,schema="newick"):
        return self.as_string(schema="newick",suppress_internal_node_labels=False,suppress_annotations=False,suppress_rooting=True)


    def __str__(self):
        return self.convertTreeToString(schema="newick")
