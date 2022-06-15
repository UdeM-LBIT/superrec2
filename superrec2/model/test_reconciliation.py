import unittest
from ete3 import Tree
from infinity import inf
from ..utils.trees import LowestCommonAncestor
from .reconciliation import (
    ReconciliationInput,
    SuperReconciliationInput,
    ReconciliationOutput,
    SuperReconciliationOutput,
    EdgeEvent,
    NodeEvent,
    EventCosts,
    EventsAtNode,
)


class TestModelReconciliation(unittest.TestCase):
    def test_rec_events_at_node(self):
        gene_tree = Tree("(x_1,z_1)O1;", format=1)
        species_tree = Tree("((X,Y)S2,Z)S1;", format=1)
        species_lca = LowestCommonAncestor(species_tree)

        costs = {
            NodeEvent.DUPLICATION: 1,
            NodeEvent.HORIZONTAL_TRANSFER: 4,
            EdgeEvent.SEGMENTAL_LOSS: 1,
            EdgeEvent.FULL_LOSS: 1,
            EdgeEvent.UNSAMPLED_TRANSFER: inf,
            EdgeEvent.TRANSFER_LOSS: inf,
        }

        leaf_object_species = {
            gene_tree & "x_1": species_tree & "X",
            gene_tree & "z_1": species_tree & "Z",
        }

        object_species = {
            **leaf_object_species,
            gene_tree & "O1": species_tree & "S1",
        }

        rec_input = ReconciliationInput(
            gene_tree,
            species_lca,
            leaf_object_species,
            EventCosts(costs),
        )

        rec_output = ReconciliationOutput(
            rec_input,
            object_species,
        )

        self.assertEqual(
            rec_output.events_at_node(gene_tree & "O1"),
            EventsAtNode(
                node_event=NodeEvent.SPECIATION,
                target_child=None,
                edge_events=[
                    {EdgeEvent.FULL_LOSS: 1},
                    {},
                ]
            )
        )


    # @classmethod
    # def setUpClass(cls):
    #     cls.gene_tree = Tree(
    #         """
    #         (
    #             ((x_1,z_1)3,(w_1,w_2)4)2,
    #             (
    #                 (
    #                     (x_2,y_4)7,
    #                     ((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8
    #                 )6,
    #                 (w_3,(z_3,(t_1,t_2)14)13)12
    #             )5
    #         )1;
    #     """,
    #         format=1,
    #     )

    #     cls.species_tree = Tree("(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;", format=1)
    #     cls.species_lca = LowestCommonAncestor(cls.species_tree)

    #     cls.leaf_object_species = {
    #         cls.gene_tree & "x_1": cls.species_tree & "X",
    #         cls.gene_tree & "x_2": cls.species_tree & "X",
    #         cls.gene_tree & "x_3": cls.species_tree & "X",
    #         cls.gene_tree & "y_1": cls.species_tree & "Y",
    #         cls.gene_tree & "y_2": cls.species_tree & "Y",
    #         cls.gene_tree & "y_3": cls.species_tree & "Y",
    #         cls.gene_tree & "y_4": cls.species_tree & "Y",
    #         cls.gene_tree & "z_1": cls.species_tree & "Z",
    #         cls.gene_tree & "z_2": cls.species_tree & "Z",
    #         cls.gene_tree & "z_3": cls.species_tree & "Z",
    #         cls.gene_tree & "w_1": cls.species_tree & "W",
    #         cls.gene_tree & "w_2": cls.species_tree & "W",
    #         cls.gene_tree & "w_3": cls.species_tree & "W",
    #         cls.gene_tree & "t_1": cls.species_tree & "T",
    #         cls.gene_tree & "t_2": cls.species_tree & "T",
    #     }

    #     cls.object_species = {
    #         **cls.leaf_object_species,
    #         cls.gene_tree & "1": cls.species_tree & "XYZWT",
    #         cls.gene_tree & "2": cls.species_tree & "XYZ",
    #         cls.gene_tree & "3": cls.species_tree & "XYZ",
    #         cls.gene_tree & "4": cls.species_tree & "W",
    #         cls.gene_tree & "5": cls.species_tree & "XYZWT",
    #         cls.gene_tree & "6": cls.species_tree & "XYZ",
    #         cls.gene_tree & "7": cls.species_tree & "XY",
    #         cls.gene_tree & "8": cls.species_tree & "XYZ",
    #         cls.gene_tree & "9": cls.species_tree & "XY",
    #         cls.gene_tree & "10": cls.species_tree & "Y",
    #         cls.gene_tree & "11": cls.species_tree & "Y",
    #         cls.gene_tree & "12": cls.species_tree & "WT",
    #         cls.gene_tree & "13": cls.species_tree & "T",
    #         cls.gene_tree & "14": cls.species_tree & "T",
    #     }

    #     cls.leaf_syntenies = {
    #         cls.gene_tree & "x_1": "abcd",
    #         cls.gene_tree & "z_1": "abcd",
    #         cls.gene_tree & "w_1": "ab",
    #         cls.gene_tree & "w_2": "abc",
    #         cls.gene_tree & "x_2": "defg",
    #         cls.gene_tree & "y_4": "def",
    #         cls.gene_tree & "x_3": "cdef",
    #         cls.gene_tree & "y_1": "ce",
    #         cls.gene_tree & "y_2": "cde",
    #         cls.gene_tree & "y_3": "cde",
    #         cls.gene_tree & "z_2": "cef",
    #         cls.gene_tree & "w_3": "defg",
    #         cls.gene_tree & "z_3": "defg",
    #         cls.gene_tree & "t_1": "def",
    #         cls.gene_tree & "t_2": "defg",
    #     }

    #     cls.syntenies = {
    #         **cls.leaf_syntenies,
    #         cls.gene_tree & "1": "abcdefg",
    #         cls.gene_tree & "2": "abcd",
    #         cls.gene_tree & "3": "abcd",
    #         cls.gene_tree & "4": "abc",
    #         cls.gene_tree & "5": "abcdefg",
    #         cls.gene_tree & "6": "abcdefg",
    #         cls.gene_tree & "7": "defg",
    #         cls.gene_tree & "8": "cdef",
    #         cls.gene_tree & "9": "cdef",
    #         cls.gene_tree & "10": "cdef",
    #         cls.gene_tree & "11": "cde",
    #         cls.gene_tree & "12": "defg",
    #         cls.gene_tree & "13": "defg",
    #         cls.gene_tree & "14": "defg",
    #     }

    #     cls.direct_costs = EventCosts({
    #         EdgeEvent.TRANSFER_LOSS: inf,
    #         EdgeEvent.UNSAMPLED_TRANSFER: inf,
    #     })

    #     cls.skip_costs = EventCosts()

    #     cls.rec_input = ReconciliationInput(
    #         cls.gene_tree,
    #         cls.species_lca,
    #         cls.leaf_object_species,
    #         cls.direct_costs,
    #     )

    #     cls.srec_input = SuperReconciliationInput(
    #         cls.gene_tree,
    #         cls.species_lca,
    #         cls.leaf_object_species,
    #         cls.leaf_syntenies,
    #         cls.direct_costs,
    #     )

    #     cls.rec_output = ReconciliationOutput(
    #         cls.rec_input,
    #         cls.object_species,
    #     )

    #     cls.srec_output = SuperReconciliationOutput(
    #         cls.rec_input, cls.object_species, cls.syntenies, ordered=True
    #     )

    #     cls.usrec_output = SuperReconciliationOutput(
    #         cls.rec_input, cls.object_species, cls.syntenies, ordered=False
    #     )

    #     print(cls.rec_output.events_at_node(cls.gene_tree))

    # def test_node_event(self):
    #     expected_events = {
    #         "1": NodeEvent.DUPLICATION,
    #         "2": NodeEvent.HORIZONTAL_TRANSFER,
    #         "3": NodeEvent.SPECIATION,
    #         "4": NodeEvent.DUPLICATION,
    #         "5": NodeEvent.SPECIATION,
    #         "6": NodeEvent.DUPLICATION,
    #         "7": NodeEvent.SPECIATION,
    #         "8": NodeEvent.SPECIATION,
    #         "9": NodeEvent.SPECIATION,
    #         "10": NodeEvent.DUPLICATION,
    #         "11": NodeEvent.DUPLICATION,
    #         "12": NodeEvent.SPECIATION,
    #         "13": NodeEvent.HORIZONTAL_TRANSFER,
    #         "14": NodeEvent.DUPLICATION,
    #     }

    #     for name, event in expected_events.items():
    #         self.assertEqual(
    #             self.rec_output.node_event(self.gene_tree & name),
    #             event,
    #         )

    # def test_labeling_cost(self):
    #     self.assertEqual(self.srec_output.labeling_cost(), 6)
    #     self.assertEqual(self.usrec_output.labeling_cost(), 5)

    # def test_reconciliation_cost(self):
    #     expected_costs = {
    #         (1, 0, 0): 6,
    #         (0, 1, 0): 2,
    #         (0, 0, 1): 3,
    #         (1, 1, 1): 11,
    #     }

    #     for (dup, hgt, loss), value in expected_costs.items():
    #         self.rec_input.costs.values[NodeEvent.DUPLICATION] = dup
    #         self.rec_input.costs.values[NodeEvent.HORIZONTAL_TRANSFER] = hgt
    #         self.rec_input.costs.values[EdgeEvent.FULL_LOSS] = loss
    #         self.assertEqual(self.rec_output.cost(), value)
