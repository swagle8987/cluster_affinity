import colorsys
from dendropy import Tree
from .counter import Counter

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

    def label_tree(self):
        for node in self.postorder_node_iter():
            if not node.label:
                node.label = get_annotation()

    def __str__(self):
        return self.convertTreeToString(schema="newick")
