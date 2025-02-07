import dendropy

class ExtendedTree(dendropy.Tree):
    def __init__(self,*args,**kwargs) -> None:
        super().__init__(*args,**kwargs)
        self.taxon_mapping = dict()
        self.compute_sizes()

    @classmethod
    def node_factory(cls, **kwargs):
        return NewNode(**kwargs)

    def compute_sizes(self):
        ind = 0
        def recursive_compute(node):
            nonlocal ind
            if node.is_leaf() and node.parent_node: ## Because the seed node is also a leaf??
                self.taxon_mapping[node.taxon.label] = node
                node.size = 1
            else:
                node.size = sum([recursive_compute(i) for i in node.child_nodes()])
            node.index = ind
            ind += 1
            return node.size
        recursive_compute(self.seed_node)


    def get_leaf_from_taxonlabel(self,l):
        return self.taxon_mapping[l]



class NewNode(dendropy.Node):
    def __init__(self,*args,**kwargs) -> None:
        super().__init__(*args,**kwargs)
        self.size = -1
        self.index = -1
        self.heavy = False

    def is_heavy(self):
        if not self._parent_node:
            return True
        elif self.parent_node.size <= 2*self.size:
            return True
        else:
            return False
