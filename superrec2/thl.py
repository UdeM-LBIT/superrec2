from collections import defaultdict
from collections.abc import Sequence
from functional import product, total_ordering
from typing import Optional, List
from enum import Enum, auto
from dataclasses import dataclass
from ete3 import PhyloTree, PhyloTreeNode


class Event(Enum):
    Leaf = auto()
    Speciation = auto()
    Duplication = auto()
    Transfer = auto()


@dataclass
@total_ordering
class THLEventInfo:
    cost: int
    event: Event
    species: PhyloTreeNode
    species_left: Optional[PhyloTreeNode] = None
    species_right: Optional[PhyloTreeNode] = None

    def __lt__(self, other):
        return self.cost < other.cost

    def __eq__(self, other):
        return self.cost == other.cost


class THLMinEvents(Sequence):
    def __init__(self):
        self._items = []

    def add(self, *infos: THLEventInfo):
        for info in infos:
            if info.cost < float('inf'):
                if not self._items or self._items[0] > info:
                    self._items = [info]
                elif self._items[0] == info:
                    self._items.append(info)

    def min(self):
        if self._items:
            return self._items[0].cost
        else:
            return float('inf')

    def __getitem__(self, idx):
        return self._items[idx]

    def __len__(self):
        return len(self._items)


def thl_mapping(gene_tree, species_tree):
    min_costs = defaultdict(THLMinEvents)

    for v in gene_tree.traverse("postorder"):
        if v.is_leaf():
            s = species_tree.get_leaves_by_name(v.species)[0]
            min_costs[(v, s)].add(
                THLEventInfo(
                    cost=0,
                    event=Event.Leaf,
                    species=s,
                )
            )
        else:
            vl, vr = v.children

            for u in species_tree.traverse("postorder"):
                options = THLMinEvents()

                if not u.is_leaf():
                    ul, ur = u.children

                    for ul_ch, ur_ch in product(ul.traverse(), ur.traverse()):
                        options.add(
                            THLEventInfo(
                                cost=(
                                    min_costs[(vl, ul_ch)]
                                    + min_costs[(vr, ur_ch)]
                                ),
                                event=Event.Speciation,
                                species=u,
                                species_left=ul_ch,
                                species_right=ur_ch,
                            ),
                            THLEventInfo(
                                cost=(
                                    min_costs[(vl, ur_ch)]
                                    + min_costs[(vr, ul_ch)]
                                ),
                                event=Event.Speciation,
                                species=u,
                                species_left=ur_ch,
                                species_right=ul_ch,
                            ),
                        )

                    for w in species_tree.traverse("postorder"):












if __name__ == "__main__":
    def get_species_name(node_name_string):
        return node_name_string.split("_")[0]

    # ga = PhyloTree("((x1,(x2,(y1,z1))),(y2,z2));")
    # gb = PhyloTree("(x1,(x2,(y1,z1)));")
    # gc = PhyloTree("(x1,(y2,z2));")
    st = PhyloTree("((x_1,(x_2,(y_1,z_1))),(y_2,z_2));", sp_naming_function=get_species_name)
    sp = PhyloTree("(x,(y,z));", sp_naming_function=get_species_name)
