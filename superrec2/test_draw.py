import unittest
import textwrap
from ete3 import Tree
from .model.synteny import parse_synteny_mapping
from .model.tree_mapping import get_species_mapping, parse_tree_mapping
from .model.reconciliation import (
    SuperReconciliationInput,
    SuperReconciliationOutput,
)
from .utils.trees import LowestCommonAncestor
from .draw import compute_layout, render_to_tikz


class TestReconciliationDraw(unittest.TestCase):
    def assertRender(
        self, gene_text, species_text, rec_text, labeling_text, expect
    ):
        gene_tree = Tree(gene_text, format=1)
        species_tree = Tree(species_text, format=1)
        srec_input = SuperReconciliationInput(
            object_tree=gene_tree,
            species_lca=LowestCommonAncestor(species_tree),
            leaf_object_species={},
            costs={},
            leaf_syntenies={},
        )
        srec_output = SuperReconciliationOutput(
            input=srec_input,
            object_species={
                **get_species_mapping(gene_tree, species_tree),
                **parse_tree_mapping(gene_tree, species_tree, rec_text),
            },
            syntenies=parse_synteny_mapping(gene_tree, labeling_text),
        )
        layout_info = compute_layout(srec_output)
        out = render_to_tikz(srec_output, layout_info)
        self.assertEqual(
            out,
            textwrap.dedent(
                r"""
                \tikzset{
                    x={1pt},
                    y={-1pt},
                    species border/.style={
                        line width={1pt},
                        shorten <={-1pt / 2 + 0.05pt},
                        shorten >={-1pt / 2 + 0.05pt},
                    },
                    species label/.style={
                        font=\bfseries,
                        midway,
                        yshift=-10,
                    },
                    branch/.style={
                        line width={0.5pt},
                        preaction={
                            draw, white, -{},
                            line width={4pt},
                            shorten <={0.5pt},
                            shorten >={0.5pt},
                        },
                    },
                    transfer branch/.style={
                        branch,
                        -Stealth,
                    },
                    loss/.style={
                        branch, dashed,
                    },
                    extant gene/.style={
                        circle, fill,
                        outer sep=0pt, inner sep=0pt,
                        minimum size={3},
                        label={
                            [font={\strut},
                             inner xsep=0pt,
                             inner ysep=2pt,
                             outer ysep=0pt,
                             fill=white]
                            below:#1
                        },
                    },
                    branch node/.style={
                        draw, fill=white,
                        outer sep=0pt, inner sep=0pt,
                        line width={0.5pt},
                        font={\strut},
                    },
                    speciation/.style={
                        branch node, rounded rectangle,
                        minimum width={10},
                        minimum height={10},
                    },
                    duplication/.style={
                        branch node, rectangle,
                        inner xsep=4pt,
                        minimum width={10},
                        minimum height={10},
                    },
                    horizontal gene transfer/.style={
                        branch node, signal, signal to=east and west,
                        minimum width={10},
                        minimum height={10},
                    },
                }
                \begin{tikzpicture}
                """
            ).lstrip()
            + expect
            + "\\end{tikzpicture}\n",
        )

    def test_speciations(self):
        self.assertRender(
            gene_text="((x_1,y_1)2,z_1)1;",
            species_text="((X,Y)XY,Z)XYZ;",
            rec_text="1:XYZ,2:XY",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (29.762999999999998,42.5) |- (92.026,10) -- (92.026,0);
                \draw[species border] (153.449,85.0) |- (124.526,10) -- (124.526,0);
                \draw[species border] (62.263,42.5) |- (92.026,42.5) -| (124.526,85.0);
                \draw[species border] (0,85.0) |- (29.762999999999998,52.5) -- (29.762999999999998,42.5);
                \draw[species border] (92.026,85.0) |- (62.263,52.5) -- (62.263,42.5);
                \draw[species border] (29.762999999999998,85.0) |- (29.762999999999998,85.0) -| (62.263,85.0);
                \draw[species border] (0,85.0) -- (0,115.0) -- node[species label] {X} (29.762999999999998,115.0) -- (29.762999999999998,85.0);
                \draw[species border] (62.263,85.0) -- (62.263,115.0) -- node[species label] {Y} (92.026,115.0) -- (92.026,85.0);
                \draw[species border] (124.526,85.0) -- (124.526,115.0) -- node[species label] {Z} (153.449,115.0) -- (153.449,85.0);
                % gene branches
                \draw[branch] (108.276,26.25) -- (108.276,0);
                \draw[branch] (46.013,42.5) |- (108.276,26.25) -| (138.9875,85.0);
                \draw[branch] (46.013,68.75) -- (46.013,42.5);
                \draw[branch] (14.881499999999999,85.0) |- (46.013,68.75) -| (77.1445,85.0);
                \draw[branch] (14.881499999999999,104.5) -- (14.881499999999999,85.0);
                \draw[branch] (77.1445,104.5) -- (77.1445,85.0);
                \draw[branch] (138.9875,104.5) -- (138.9875,85.0);
                % gene transfers
                % events
                \node[speciation] at (108.276,26.25) {};
                \node[speciation] at (46.013,68.75) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,96.5) {};
                \node[extant gene={y\textsubscript{1}}] at (77.1445,96.5) {};
                \node[extant gene={z\textsubscript{1}}] at (138.9875,96.5) {};
                """
            ).lstrip(),
        )

    def test_duplications(self):
        self.assertRender(
            gene_text="""(
                ((x_3,y_3)8,z_3)7,
                (((x_1,y_2)4,z_2)3,((x_2,y_1)6,z_1)5)2
            )1;""",
            species_text="((X,Y)XY,Z)XYZ;",
            rec_text="1:XYZ,2:XYZ,3:XYZ,4:XY,5:XYZ,6:XY,7:XYZ,8:XY",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (59.288999999999994,112.5) |- (186.07799999999997,45.0) -- (186.07799999999997,0);
                \draw[species border] (310.347,190.0) |- (253.57799999999997,45.0) -- (253.57799999999997,0);
                \draw[species border] (126.78899999999999,112.5) |- (186.07799999999997,112.5) -| (253.57799999999997,190.0);
                \draw[species border] (0,190.0) |- (59.288999999999994,122.5) -- (59.288999999999994,112.5);
                \draw[species border] (186.07799999999997,190.0) |- (126.78899999999999,122.5) -- (126.78899999999999,112.5);
                \draw[species border] (59.288999999999994,190.0) |- (59.288999999999994,190.0) -| (126.78899999999999,190.0);
                \draw[species border] (0,190.0) -- (0,220.0) -- node[species label] {X} (59.288999999999994,220.0) -- (59.288999999999994,190.0);
                \draw[species border] (126.78899999999999,190.0) -- (126.78899999999999,220.0) -- node[species label] {Y} (186.07799999999997,220.0) -- (186.07799999999997,190.0);
                \draw[species border] (253.57799999999997,190.0) -- (253.57799999999997,220.0) -- node[species label] {Z} (310.347,220.0) -- (310.347,190.0);
                % gene branches
                \draw[branch] (75.53899999999999,112.5) |- (202.32799999999997,61.25) -| (268.0395,190.0);
                \draw[branch] (93.03899999999999,112.5) |- (219.82799999999997,78.75) -| (281.9625,190.0);
                \draw[branch] (110.53899999999999,112.5) |- (237.32799999999997,96.25) -| (295.8855,190.0);
                \draw[branch] (219.82799999999997,78.75) |- (228.57799999999997,38.75) -| (237.32799999999997,96.25);
                \draw[branch] (215.45299999999997,16.25) -- (215.45299999999997,0);
                \draw[branch] (202.32799999999997,61.25) |- (215.45299999999997,16.25) -| (228.57799999999997,38.75);
                \draw[branch] (75.53899999999999,138.75) -- (75.53899999999999,112.5);
                \draw[branch] (14.881499999999999,190.0) |- (75.53899999999999,138.75) -| (141.67049999999998,190.0);
                \draw[branch] (93.03899999999999,156.25) -- (93.03899999999999,112.5);
                \draw[branch] (29.6445,190.0) |- (93.03899999999999,156.25) -| (156.43349999999998,190.0);
                \draw[branch] (110.53899999999999,173.75) -- (110.53899999999999,112.5);
                \draw[branch] (44.4075,190.0) |- (110.53899999999999,173.75) -| (171.1965,190.0);
                \draw[branch] (14.881499999999999,209.5) -- (14.881499999999999,190.0);
                \draw[branch] (29.644499999999997,209.5) -- (29.6445,190.0);
                \draw[branch] (44.4075,209.5) -- (44.4075,190.0);
                \draw[branch] (141.67049999999998,209.5) -- (141.67049999999998,190.0);
                \draw[branch] (156.43349999999998,209.5) -- (156.43349999999998,190.0);
                \draw[branch] (171.1965,209.5) -- (171.1965,190.0);
                \draw[branch] (268.0395,209.5) -- (268.0395,190.0);
                \draw[branch] (281.9625,209.5) -- (281.9625,190.0);
                \draw[branch] (295.8855,209.5) -- (295.8855,190.0);
                % gene transfers
                % events
                \node[speciation] at (202.32799999999997,61.25) {};
                \node[speciation] at (219.82799999999997,78.75) {};
                \node[speciation] at (237.32799999999997,96.25) {};
                \node[duplication] at (228.57799999999997,38.75) {};
                \node[duplication] at (215.45299999999997,16.25) {};
                \node[speciation] at (75.53899999999999,138.75) {};
                \node[speciation] at (93.03899999999999,156.25) {};
                \node[speciation] at (110.53899999999999,173.75) {};
                \node[extant gene={x\textsubscript{3}}] at (14.881499999999999,201.5) {};
                \node[extant gene={x\textsubscript{1}}] at (29.644499999999997,201.5) {};
                \node[extant gene={x\textsubscript{2}}] at (44.4075,201.5) {};
                \node[extant gene={y\textsubscript{3}}] at (141.67049999999998,201.5) {};
                \node[extant gene={y\textsubscript{2}}] at (156.43349999999998,201.5) {};
                \node[extant gene={y\textsubscript{1}}] at (171.1965,201.5) {};
                \node[extant gene={z\textsubscript{3}}] at (268.0395,201.5) {};
                \node[extant gene={z\textsubscript{2}}] at (281.9625,201.5) {};
                \node[extant gene={z\textsubscript{1}}] at (295.8855,201.5) {};
                """
            ).lstrip(),
        )

    def test_speciations_losses(self):
        self.assertRender(
            gene_text="(x_1,z_1)1;",
            species_text="(X,(Y,(Z,W)ZW)YZW)XYZW;",
            rec_text="1:XYZW",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (0,102.5) |- (29.762999999999998,10) -- (29.762999999999998,0);
                \draw[species border] (102.263,42.5) |- (62.263,10) -- (62.263,0);
                \draw[species border] (29.762999999999998,102.5) |- (29.762999999999998,42.5) -| (82.263,42.5);
                \draw[species border] (0,102.5) -- (0,132.5) -- node[species label] {X} (29.762999999999998,132.5) -- (29.762999999999998,102.5);
                \draw[species border] (62.263,121.5) |- (82.263,52.5) -- (82.263,42.5);
                \draw[species border] (151.186,72.5) |- (102.263,52.5) -- (102.263,42.5);
                \draw[species border] (82.263,121.5) |- (82.263,72.5) -| (131.186,72.5);
                \draw[species border] (62.263,121.5) -- (62.263,132.5) -- node[species label] {Y} (82.263,132.5) -- (82.263,121.5);
                \draw[species border] (102.263,102.5) |- (131.186,82.5) -- (131.186,72.5);
                \draw[species border] (171.186,121.5) |- (151.186,82.5) -- (151.186,72.5);
                \draw[species border] (131.186,102.5) |- (131.186,102.5) -| (151.186,121.5);
                \draw[species border] (102.263,102.5) -- (102.263,132.5) -- node[species label] {Z} (131.186,132.5) -- (131.186,102.5);
                \draw[species border] (151.186,121.5) -- (151.186,132.5) -- node[species label] {W} (171.186,132.5) -- (171.186,121.5);
                % gene branches
                \draw[branch] (46.013,26.25) -- (46.013,0);
                \draw[branch] (14.881499999999999,102.5) |- (46.013,26.25) -| (92.263,42.5);
                \draw[branch] (14.881499999999999,122.0) -- (14.881499999999999,102.5);
                \draw[branch] (92.263,62.5) -- (92.263,42.5);
                \draw[loss] (92.263,62.5) -- ++(-20, 0);
                \draw[branch] (92.263,62.5) -| (141.186,72.5);
                \draw[branch] (141.186,92.5) -- (141.186,72.5);
                \draw[loss] (141.186,92.5) -- ++(20, 0);
                \draw[branch] (141.186,92.5) -| (116.7245,102.5);
                \draw[branch] (116.7245,122.0) -- (116.7245,102.5);
                % gene transfers
                % events
                \node[speciation] at (46.013,26.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,114.0) {};
                \node[extant gene={z\textsubscript{1}}] at (116.7245,114.0) {};
                """
            ).lstrip(),
        )

    def test_speciation_swapped(self):
        self.maxDiff = None
        self.assertRender(
            gene_text="(w_1,x_1)1;",
            species_text="((X,Y)XY,(Z,W)ZW)XYZW;",
            rec_text="1:XYZW",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (29.762999999999998,42.5) |- (69.763,10) -- (69.763,0);
                \draw[species border] (142.263,42.5) |- (102.263,10) -- (102.263,0);
                \draw[species border] (49.763,42.5) |- (69.763,42.5) -| (122.263,42.5);
                \draw[species border] (0,72.5) |- (29.762999999999998,52.5) -- (29.762999999999998,42.5);
                \draw[species border] (69.763,91.5) |- (49.763,52.5) -- (49.763,42.5);
                \draw[species border] (29.762999999999998,72.5) |- (29.762999999999998,72.5) -| (49.763,91.5);
                \draw[species border] (0,72.5) -- (0,102.5) -- node[species label] {X} (29.762999999999998,102.5) -- (29.762999999999998,72.5);
                \draw[species border] (49.763,91.5) -- (49.763,102.5) -- node[species label] {Y} (69.763,102.5) -- (69.763,91.5);
                \draw[species border] (102.263,91.5) |- (122.263,52.5) -- (122.263,42.5);
                \draw[species border] (173.966,72.5) |- (142.263,52.5) -- (142.263,42.5);
                \draw[species border] (122.263,91.5) |- (122.263,72.5) -| (142.263,72.5);
                \draw[species border] (102.263,91.5) -- (102.263,102.5) -- node[species label] {Z} (122.263,102.5) -- (122.263,91.5);
                \draw[species border] (142.263,72.5) -- (142.263,102.5) -- node[species label] {W} (173.966,102.5) -- (173.966,72.5);
                % gene branches
                \draw[branch] (86.013,26.25) -- (86.013,0);
                \draw[branch] (39.763,42.5) |- (86.013,26.25) -| (132.263,42.5);
                \draw[branch] (39.763,62.5) -- (39.763,42.5);
                \draw[loss] (39.763,62.5) -- ++(20, 0);
                \draw[branch] (39.763,62.5) -| (14.881499999999999,72.5);
                \draw[branch] (14.881499999999999,92.0) -- (14.881499999999999,72.5);
                \draw[branch] (132.263,62.5) -- (132.263,42.5);
                \draw[loss] (132.263,62.5) -- ++(-20, 0);
                \draw[branch] (132.263,62.5) -| (158.1145,72.5);
                \draw[branch] (158.1145,92.0) -- (158.1145,72.5);
                % gene transfers
                % events
                \node[speciation] at (86.013,26.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,84.0) {};
                \node[extant gene={w\textsubscript{1}}] at (158.1145,84.0) {};
                """
            ).lstrip(),
        )

    def test_duplications_losses(self):
        self.assertRender(
            gene_text="(z_3,(((x_1,y_2)4,z_2)3,((x_2,y_1)6,z_1)5)2)1;",
            species_text="((X,Y)XY,Z)XYZ;",
            rec_text="1:XYZ,2:XYZ,3:XYZ,4:XY,5:XYZ,6:XY",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (44.525999999999996,100.0) |- (139.052,45.0) -- (139.052,0);
                \draw[species border] (250.821,160.0) |- (194.052,45.0) -- (194.052,0);
                \draw[species border] (94.526,100.0) |- (139.052,100.0) -| (194.052,160.0);
                \draw[species border] (0,160.0) |- (44.525999999999996,110.0) -- (44.525999999999996,100.0);
                \draw[species border] (139.052,160.0) |- (94.526,110.0) -- (94.526,100.0);
                \draw[species border] (44.525999999999996,160.0) |- (44.525999999999996,160.0) -| (94.526,160.0);
                \draw[species border] (0,160.0) -- (0,190.0) -- node[species label] {X} (44.525999999999996,190.0) -- (44.525999999999996,160.0);
                \draw[species border] (94.526,160.0) -- (94.526,190.0) -- node[species label] {Y} (139.052,190.0) -- (139.052,160.0);
                \draw[species border] (194.052,160.0) -- (194.052,190.0) -- node[species label] {Z} (250.821,190.0) -- (250.821,160.0);
                % gene branches
                \draw[branch] (60.775999999999996,100.0) |- (155.302,61.25) -| (222.4365,160.0);
                \draw[branch] (78.276,100.0) |- (172.802,78.75) -| (236.3595,160.0);
                \draw[branch] (155.302,61.25) |- (164.052,38.75) -| (172.802,78.75);
                \draw[loss] (184.052,90.0) -- ++(-20, 0);
                \draw[branch] (184.052,90.0) -| (208.5135,160.0);
                \draw[branch] (174.052,16.25) -- (174.052,0);
                \draw[branch] (184.052,90.0) |- (174.052,16.25) -| (164.052,38.75);
                \draw[branch] (60.775999999999996,126.25) -- (60.775999999999996,100.0);
                \draw[branch] (14.881499999999999,160.0) |- (60.775999999999996,126.25) -| (109.4075,160.0);
                \draw[branch] (78.276,143.75) -- (78.276,100.0);
                \draw[branch] (29.6445,160.0) |- (78.276,143.75) -| (124.1705,160.0);
                \draw[branch] (14.881499999999999,179.5) -- (14.881499999999999,160.0);
                \draw[branch] (29.644499999999997,179.5) -- (29.6445,160.0);
                \draw[branch] (109.4075,179.5) -- (109.4075,160.0);
                \draw[branch] (124.17049999999999,179.5) -- (124.1705,160.0);
                \draw[branch] (208.5135,179.5) -- (208.5135,160.0);
                \draw[branch] (222.4365,179.5) -- (222.4365,160.0);
                \draw[branch] (236.3595,179.5) -- (236.3595,160.0);
                % gene transfers
                % events
                \node[speciation] at (155.302,61.25) {};
                \node[speciation] at (172.802,78.75) {};
                \node[duplication] at (164.052,38.75) {};
                \node[duplication] at (174.052,16.25) {};
                \node[speciation] at (60.775999999999996,126.25) {};
                \node[speciation] at (78.276,143.75) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,171.5) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,171.5) {};
                \node[extant gene={y\textsubscript{2}}] at (109.4075,171.5) {};
                \node[extant gene={y\textsubscript{1}}] at (124.17049999999999,171.5) {};
                \node[extant gene={z\textsubscript{3}}] at (208.5135,171.5) {};
                \node[extant gene={z\textsubscript{2}}] at (222.4365,171.5) {};
                \node[extant gene={z\textsubscript{1}}] at (236.3595,171.5) {};
                """
            ).lstrip(),
        )

    def test_nested_duplications(self):
        self.assertRender(
            gene_text="(((x_1,y_1)4,(x_2,y_2)5)2,((x_3,y_3)6,(x_4,y_4)7)3)1;",
            species_text="(X,Y)XY;",
            rec_text="1:XY,2:XY,3:XY,4:XY,5:XY,6:XY,7:XY",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (0,130.0) |- (74.052,45.0) -- (74.052,0);
                \draw[species border] (233.10400000000004,130.0) |- (159.05200000000002,45.0) -- (159.05200000000002,0);
                \draw[species border] (74.052,130.0) |- (74.052,130.0) -| (159.05200000000002,130.0);
                \draw[species border] (0,130.0) -- (0,160.0) -- node[species label] {X} (74.052,160.0) -- (74.052,130.0);
                \draw[species border] (159.05200000000002,130.0) -- (159.05200000000002,160.0) -- node[species label] {Y} (233.10400000000004,160.0) -- (233.10400000000004,130.0);
                % gene branches
                \draw[branch] (14.881499999999999,130.0) |- (90.302,61.25) -| (173.9335,130.0);
                \draw[branch] (29.6445,130.0) |- (107.802,78.75) -| (188.69650000000001,130.0);
                \draw[branch] (90.302,61.25) |- (99.052,38.75) -| (107.802,78.75);
                \draw[branch] (44.4075,130.0) |- (125.302,96.25) -| (203.45950000000002,130.0);
                \draw[branch] (59.170500000000004,130.0) |- (142.80200000000002,113.75) -| (218.22250000000003,130.0);
                \draw[branch] (125.302,96.25) |- (134.05200000000002,38.75) -| (142.80200000000002,113.75);
                \draw[branch] (116.552,16.25) -- (116.552,0);
                \draw[branch] (99.052,38.75) |- (116.552,16.25) -| (134.05200000000002,38.75);
                \draw[branch] (14.881499999999999,149.5) -- (14.881499999999999,130.0);
                \draw[branch] (29.644499999999997,149.5) -- (29.6445,130.0);
                \draw[branch] (44.4075,149.5) -- (44.4075,130.0);
                \draw[branch] (59.170500000000004,149.5) -- (59.170500000000004,130.0);
                \draw[branch] (173.9335,149.5) -- (173.9335,130.0);
                \draw[branch] (188.69650000000001,149.5) -- (188.69650000000001,130.0);
                \draw[branch] (203.45950000000002,149.5) -- (203.45950000000002,130.0);
                \draw[branch] (218.22250000000003,149.5) -- (218.22250000000003,130.0);
                % gene transfers
                % events
                \node[speciation] at (90.302,61.25) {};
                \node[speciation] at (107.802,78.75) {};
                \node[duplication] at (99.052,38.75) {};
                \node[speciation] at (125.302,96.25) {};
                \node[speciation] at (142.80200000000002,113.75) {};
                \node[duplication] at (134.05200000000002,38.75) {};
                \node[duplication] at (116.552,16.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,141.5) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,141.5) {};
                \node[extant gene={x\textsubscript{3}}] at (44.4075,141.5) {};
                \node[extant gene={x\textsubscript{4}}] at (59.170500000000004,141.5) {};
                \node[extant gene={y\textsubscript{1}}] at (173.9335,141.5) {};
                \node[extant gene={y\textsubscript{2}}] at (188.69650000000001,141.5) {};
                \node[extant gene={y\textsubscript{3}}] at (203.45950000000002,141.5) {};
                \node[extant gene={y\textsubscript{4}}] at (218.22250000000003,141.5) {};
                """
            ).lstrip(),
        )

    def test_duplication_swapped(self):
        self.assertRender(
            gene_text="(y_1,x_1)1;",
            species_text="((X,Y)XY,Z)XYZ;",
            rec_text="1:XYZ",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (29.762999999999998,47.5) |- (84.526,22.5) -- (84.526,0);
                \draw[species border] (135.026,101.5) |- (115.026,22.5) -- (115.026,0);
                \draw[species border] (54.763,47.5) |- (84.526,47.5) -| (115.026,101.5);
                \draw[species border] (0,82.5) |- (29.762999999999998,57.5) -- (29.762999999999998,47.5);
                \draw[species border] (84.526,82.5) |- (54.763,57.5) -- (54.763,47.5);
                \draw[species border] (29.762999999999998,82.5) |- (29.762999999999998,82.5) -| (54.763,82.5);
                \draw[species border] (0,82.5) -- (0,112.5) -- node[species label] {X} (29.762999999999998,112.5) -- (29.762999999999998,82.5);
                \draw[species border] (54.763,82.5) -- (54.763,112.5) -- node[species label] {Y} (84.526,112.5) -- (84.526,82.5);
                \draw[species border] (115.026,101.5) -- (115.026,112.5) -- node[species label] {Z} (135.026,112.5) -- (135.026,101.5);
                % gene branches
                \draw[loss] (97.276,32.5) -- ++(20, 0);
                \draw[branch] (97.276,32.5) -| (39.763,47.5);
                \draw[loss] (102.276,37.5) -- ++(20, 0);
                \draw[branch] (102.276,37.5) -| (44.763,47.5);
                \draw[branch] (99.776,16.25) -- (99.776,0);
                \draw[branch] (97.276,32.5) |- (99.776,16.25) -| (102.276,37.5);
                \draw[branch] (39.763,67.5) -- (39.763,47.5);
                \draw[loss] (39.763,67.5) -- ++(-20, 0);
                \draw[branch] (39.763,67.5) -| (69.6445,82.5);
                \draw[branch] (44.763,72.5) -- (44.763,47.5);
                \draw[loss] (44.763,72.5) -- ++(20, 0);
                \draw[branch] (44.763,72.5) -| (14.881499999999999,82.5);
                \draw[branch] (14.881499999999999,102.0) -- (14.881499999999999,82.5);
                \draw[branch] (69.6445,102.0) -- (69.6445,82.5);
                % gene transfers
                % events
                \node[duplication] at (99.776,16.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,94.0) {};
                \node[extant gene={y\textsubscript{1}}] at (69.6445,94.0) {};
                """
            ).lstrip(),
        )

    def test_transfers(self):
        self.assertRender(
            gene_text="((x_1,y_1)2,(((x_2,y_2)5,z_1)4,z_2)3)1;",
            species_text="((X,Y)XY,Z)XYZ;",
            rec_text="1:XYZ,2:XY,3:XYZ,4:XY,5:XY",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (44.525999999999996,60.0) |- (139.052275,22.5) -- (139.052275,0);
                \draw[species border] (219.398275,132.5) |- (176.552275,22.5) -- (176.552275,0);
                \draw[species border] (94.526275,60.0) |- (139.052275,60.0) -| (176.552275,132.5);
                \draw[species border] (0,132.5) |- (44.525999999999996,82.5) -- (44.525999999999996,60.0);
                \draw[species border] (139.052275,132.5) |- (94.526275,82.5) -- (94.526275,60.0);
                \draw[species border] (44.525999999999996,132.5) |- (44.525999999999996,132.5) -| (94.526275,132.5);
                \draw[species border] (0,132.5) -- (0,162.5) -- node[species label] {X} (44.525999999999996,162.5) -- (44.525999999999996,132.5);
                \draw[species border] (94.526275,132.5) -- (94.526275,162.5) -- node[species label] {Y} (139.052275,162.5) -- (139.052275,132.5);
                \draw[species border] (176.552275,132.5) -- (176.552275,162.5) -- node[species label] {Z} (219.398275,162.5) -- (219.398275,132.5);
                % gene branches
                \draw[branch] (78.276,60.0) |- (155.302275,38.75) -| (204.936775,132.5);
                \draw[loss] (166.552275,50.0) -- ++(20, 0);
                \draw[branch] (166.552275,50.0) -| (60.775999999999996,60.0);
                \draw[branch] (160.927275,16.25) -- (160.927275,0);
                \draw[branch] (166.552275,50.0) |- (160.927275,16.25) -| (155.302275,38.75);
                \draw[branch] (60.775999999999996,98.75) -- (60.775999999999996,60.0);
                \draw[branch] (14.881499999999999,132.5) |- (60.775999999999996,98.75) -| (109.407775,132.5);
                \draw[branch] (29.6445,132.5) |- (78.276,116.25) -| (124.17077499999999,132.5);
                \draw[branch] (78.276,76.25) -- (78.276,60.0);
                \draw[branch] (78.276,116.25) |- (78.276,76.25);
                \draw[branch] (14.881499999999999,152.0) -- (14.881499999999999,132.5);
                \draw[branch] (29.644499999999997,152.0) -- (29.6445,132.5);
                \draw[branch] (109.407775,152.0) -- (109.407775,132.5);
                \draw[branch] (124.17077499999999,152.0) -- (124.17077499999999,132.5);
                \draw[branch] (191.013775,152.0) -- (191.013775,132.5);
                \draw[branch] (204.936775,152.0) -- (204.936775,132.5);
                % gene transfers
                \draw[transfer branch] (78.276,76.25) to[bend left=35] (191.013775,132.5);
                % events
                \node[speciation] at (155.302275,38.75) {};
                \node[duplication] at (160.927275,16.25) {};
                \node[speciation] at (60.775999999999996,98.75) {};
                \node[speciation] at (78.276,116.25) {};
                \node[horizontal gene transfer] at (78.276,76.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,144.0) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,144.0) {};
                \node[extant gene={y\textsubscript{1}}] at (109.407775,144.0) {};
                \node[extant gene={y\textsubscript{2}}] at (124.17077499999999,144.0) {};
                \node[extant gene={z\textsubscript{1}}] at (191.013775,144.0) {};
                \node[extant gene={z\textsubscript{2}}] at (204.936775,144.0) {};
                """
            ).lstrip(),
        )

    def test_empty_root(self):
        self.assertRender(
            gene_text="(x_1,(y_1,z_1)2)1;",
            species_text="(X,(Y,Z)YZ)XYZ;",
            rec_text="1:YZ,2:YZ",
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (0,75.0) |- (29.762999999999998,10) -- (29.762999999999998,0);
                \draw[species border] (112.02655,20.0) |- (49.763,10) -- (49.763,0);
                \draw[species border] (29.762999999999998,75.0) |- (29.762999999999998,20) -| (79.526,20.0);
                \draw[species border] (0,75.0) -- (0,105.0) -- node[species label] {X} (29.762999999999998,105.0) -- (29.762999999999998,75.0);
                \draw[species border] (49.763,75.0) |- (79.526,42.5) -- (79.526,20.0);
                \draw[species border] (140.94955,75.0) |- (112.02655,42.5) -- (112.02655,20.0);
                \draw[species border] (79.526,75.0) |- (79.526,75.0) -| (112.02655,75.0);
                \draw[species border] (49.763,75.0) -- (49.763,105.0) -- node[species label] {Y} (79.526,105.0) -- (79.526,75.0);
                \draw[species border] (112.02655,75.0) -- (112.02655,105.0) -- node[species label] {Z} (140.94955,105.0) -- (140.94955,75.0);
                % gene branches
                \draw[branch] (14.881499999999999,94.5) -- (14.881499999999999,75.0);
                \draw[branch] (64.6445,75.0) |- (95.776275,58.75) -| (126.48805,75.0);
                \draw[branch] (95.776275,36.25) -- (95.776275,20.0);
                \draw[branch] (95.776275,58.75) |- (95.776275,36.25);
                \draw[branch] (64.6445,94.5) -- (64.6445,75.0);
                \draw[branch] (126.48805,94.5) -- (126.48805,75.0);
                % gene transfers
                \draw[transfer branch] (95.776275,36.25) to[bend right=35] (14.881499999999999,75.0);
                % events
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,86.5) {};
                \node[speciation] at (95.776275,58.75) {};
                \node[horizontal gene transfer] at (95.776275,36.25) {};
                \node[extant gene={y\textsubscript{1}}] at (64.6445,86.5) {};
                \node[extant gene={z\textsubscript{1}}] at (126.48805,86.5) {};
                """
            ).lstrip(),
        )

    def test_all(self):
        self.assertRender(
            gene_text="""(
                ((x_1,z_1)3,(w_1,w_2)4)2,
                (
                    ((x_2,y_4)7,((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8)6,
                    (w_3,(z_3,(t_1,t_2)14)13)12
                )5
            )1;""",
            species_text="(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;",
            rec_text=(
                "1:XYZWT,2:XYZ,3:XYZ,4:W,5:XYZWT,6:XYZ,7:XY,"
                "8:XYZ,9:XY,10:Y,11:Y,12:WT,13:T,14:T"
            ),
            labeling_text="",
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (188.341,60.0) |- (300.110275,22.5) -- (300.110275,0);
                \draw[species border] (435.219275,160.0) |- (337.610275,22.5) -- (337.610275,0);
                \draw[species border] (243.341275,60.0) |- (300.110275,60.0) -| (402.719275,160.0);
                \draw[species border] (59.288999999999994,137.5) |- (188.341,82.5) -- (188.341,60.0);
                \draw[species border] (300.110275,247.5) |- (243.341275,82.5) -- (243.341275,60.0);
                \draw[species border] (114.28899999999999,137.5) |- (188.341,137.5) -| (243.341275,247.5);
                \draw[species border] (0,247.5) |- (59.288999999999994,147.5) -- (59.288999999999994,137.5);
                \draw[species border] (188.341,202.5) |- (114.28899999999999,147.5) -- (114.28899999999999,137.5);
                \draw[species border] (59.288999999999994,247.5) |- (59.288999999999994,202.5) -| (114.28899999999999,202.5);
                \draw[species border] (0,247.5) -- (0,277.5) -- node[species label] {X} (59.288999999999994,277.5) -- (59.288999999999994,247.5);
                \draw[species border] (114.28899999999999,202.5) -- (114.28899999999999,277.5) -- node[species label] {Y} (188.341,277.5) -- (188.341,202.5);
                \draw[species border] (243.341275,247.5) -- (243.341275,277.5) -- node[species label] {Z} (300.110275,277.5) -- (300.110275,247.5);
                \draw[species border] (337.610275,225.0) |- (402.719275,170.0) -- (402.719275,160.0);
                \draw[species border] (476.96525499999996,202.5) |- (435.219275,170.0) -- (435.219275,160.0);
                \draw[species border] (402.719275,225.0) |- (402.719275,202.5) -| (435.219275,202.5);
                \draw[species border] (337.610275,225.0) -- (337.610275,277.5) -- node[species label] {W} (402.719275,277.5) -- (402.719275,225.0);
                \draw[species border] (435.219275,202.5) -- (435.219275,277.5) -- node[species label] {T} (476.96525499999996,277.5) -- (476.96525499999996,202.5);
                % gene branches
                \draw[branch] (227.716275,60.0) |- (316.360275,38.75) -| (418.969275,160.0);
                \draw[loss] (327.610275,50.0) -- ++(20, 0);
                \draw[branch] (327.610275,50.0) -| (204.591275,60.0);
                \draw[branch] (321.985275,16.25) -- (321.985275,0);
                \draw[branch] (327.610275,50.0) |- (321.985275,16.25) -| (316.360275,38.75);
                \draw[branch] (104.28899999999999,137.5) |- (204.591275,98.75) -| (257.802775,247.5);
                \draw[branch] (204.591275,76.25) -- (204.591275,60.0);
                \draw[branch] (204.591275,98.75) |- (204.591275,76.25);
                \draw[branch] (93.03899999999999,137.5) |- (222.091275,116.25) -| (271.725775,247.5);
                \draw[loss] (233.341275,127.5) -- ++(20, 0);
                \draw[branch] (233.341275,127.5) -| (75.53899999999999,137.5);
                \draw[branch] (227.716275,76.25) -- (227.716275,60.0);
                \draw[branch] (233.341275,127.5) |- (227.716275,76.25) -| (222.091275,116.25);
                \draw[branch] (75.53899999999999,163.75) -- (75.53899999999999,137.5);
                \draw[branch] (29.6445,247.5) |- (75.53899999999999,163.75) -| (129.17049999999998,202.5);
                \draw[branch] (93.03899999999999,181.25) -- (93.03899999999999,137.5);
                \draw[branch] (44.4075,247.5) |- (93.03899999999999,181.25) -| (155.00574999999998,202.5);
                \draw[branch] (104.28899999999999,192.5) -- (104.28899999999999,137.5);
                \draw[loss] (104.28899999999999,192.5) -- ++(20, 0);
                \draw[branch] (104.28899999999999,192.5) -| (14.881499999999999,247.5);
                \draw[branch] (14.881499999999999,267.0) -- (14.881499999999999,247.5);
                \draw[branch] (29.644499999999997,267.0) -- (29.6445,247.5);
                \draw[branch] (44.4075,267.0) -- (44.4075,247.5);
                \draw[branch] (129.17049999999998,267.0) -- (129.17049999999998,202.5);
                \draw[branch] (158.6965,267.0) |- (166.07799999999997,241.25) -| (173.4595,267.0);
                \draw[branch] (155.00574999999998,218.75) -- (155.00574999999998,202.5);
                \draw[branch] (143.93349999999998,267.0) |- (155.00574999999998,218.75) -| (166.07799999999997,241.25);
                \draw[branch] (257.802775,267.0) -- (257.802775,247.5);
                \draw[branch] (271.725775,267.0) -- (271.725775,247.5);
                \draw[branch] (285.648775,267.0) -- (285.648775,247.5);
                \draw[branch] (418.969275,186.25) -- (418.969275,160.0);
                \draw[branch] (386.867775,225.0) |- (418.969275,186.25) -| (456.092265,202.5);
                \draw[branch] (361.813275,241.25) -- (361.813275,225.0);
                \draw[branch] (353.461775,267.0) |- (361.813275,241.25) -| (370.164775,267.0);
                \draw[branch] (386.867775,267.0) -- (386.867775,225.0);
                \draw[branch] (449.40576999999996,267.0) |- (456.092265,241.25) -| (462.77876,267.0);
                \draw[branch] (456.092265,218.75) -- (456.092265,202.5);
                \draw[branch] (456.092265,241.25) |- (456.092265,218.75);
                % gene transfers
                \draw[transfer branch] (204.591275,76.25) to[bend left=35] (361.813275,225.0);
                \draw[transfer branch] (456.092265,218.75) to[bend right=35] (285.648775,247.5);
                % events
                \node[speciation] at (316.360275,38.75) {};
                \node[duplication] at (321.985275,16.25) {};
                \node[speciation] at (204.591275,98.75) {};
                \node[horizontal gene transfer] at (204.591275,76.25) {};
                \node[speciation] at (222.091275,116.25) {};
                \node[duplication] at (227.716275,76.25) {};
                \node[speciation] at (75.53899999999999,163.75) {};
                \node[speciation] at (93.03899999999999,181.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,259.0) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,259.0) {};
                \node[extant gene={x\textsubscript{3}}] at (44.4075,259.0) {};
                \node[extant gene={y\textsubscript{4}}] at (129.17049999999998,259.0) {};
                \node[extant gene={y\textsubscript{1}}] at (143.93349999999998,259.0) {};
                \node[extant gene={y\textsubscript{2}}] at (158.6965,259.0) {};
                \node[extant gene={y\textsubscript{3}}] at (173.4595,259.0) {};
                \node[duplication] at (166.07799999999997,241.25) {};
                \node[duplication] at (155.00574999999998,218.75) {};
                \node[extant gene={z\textsubscript{1}}] at (257.802775,259.0) {};
                \node[extant gene={z\textsubscript{2}}] at (271.725775,259.0) {};
                \node[extant gene={z\textsubscript{3}}] at (285.648775,259.0) {};
                \node[speciation] at (418.969275,186.25) {};
                \node[extant gene={w\textsubscript{1}}] at (353.461775,259.0) {};
                \node[extant gene={w\textsubscript{2}}] at (370.164775,259.0) {};
                \node[duplication] at (361.813275,241.25) {};
                \node[extant gene={w\textsubscript{3}}] at (386.867775,259.0) {};
                \node[extant gene={t\textsubscript{1}}] at (449.40576999999996,259.0) {};
                \node[extant gene={t\textsubscript{2}}] at (462.77876,259.0) {};
                \node[duplication] at (456.092265,241.25) {};
                \node[horizontal gene transfer] at (456.092265,218.75) {};
                """
            ).lstrip(),
        )

    def test_syntenies(self):
        self.assertRender(
            gene_text="""(
                ((x_1,z_1)3,(w_1,w_2)4)2,
                (
                    ((x_2,y_4)7,((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8)6,
                    (w_3,(z_3,(t_1,t_2)14)13)12
                )5
            )1;""",
            species_text="(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;",
            rec_text=(
                "1:XYZWT,2:XYZ,3:XYZ,4:W,5:XYZWT,6:XYZ,7:XY,"
                "8:XYZ,9:XY,10:Y,11:Y,12:WT,13:T,14:T"
            ),
            labeling_text=(
                "1:abcdefg,2:abcd,3:abcd,x_1:abcd,z_1:abcd,4:abc,w_1:ab,"
                "w_2:abc,5:abcdefg,6:abcdefg,7:defg,x_2:defg,y_4:def,8:cdef,"
                "9:cdef,x_3:cdef,10:cdef,y_1:cef,11:cde,y_2:cde,y_3:cde,"
                "z_2:cef,12:defg,w_3:defg,13:defg,z_3:defg,14:defg,t_1:def,"
                "t_2:defg"
            ),
            expect=textwrap.dedent(
                r"""
                % species
                \draw[species border] (265.84000000000003,60.0) |- (450.94027500000004,22.5) -- (450.94027500000004,0);
                \draw[species border] (653.2002750000001,160.0) |- (528.7402750000001,22.5) -- (528.7402750000001,0);
                \draw[species border] (370.100275,60.0) |- (450.94027500000004,60.0) -| (602.6402750000001,160.0);
                \draw[species border] (86.4,137.5) |- (265.84000000000003,82.5) -- (265.84000000000003,60.0);
                \draw[species border] (450.94027500000004,247.5) |- (370.100275,82.5) -- (370.100275,60.0);
                \draw[species border] (176.96,137.5) |- (265.84000000000003,137.5) -| (370.100275,247.5);
                \draw[species border] (0,247.5) |- (86.4,147.5) -- (86.4,137.5);
                \draw[species border] (265.84000000000003,202.5) |- (176.96,147.5) -- (176.96,137.5);
                \draw[species border] (86.4,247.5) |- (86.4,202.5) -| (176.96,202.5);
                \draw[species border] (0,247.5) -- (0,277.5) -- node[species label] {X} (86.4,277.5) -- (86.4,247.5);
                \draw[species border] (176.96,202.5) -- (176.96,277.5) -- node[species label] {Y} (265.84000000000003,277.5) -- (265.84000000000003,202.5);
                \draw[species border] (370.100275,247.5) -- (370.100275,277.5) -- node[species label] {Z} (450.94027500000004,277.5) -- (450.94027500000004,247.5);
                \draw[species border] (528.7402750000001,225.0) |- (602.6402750000001,170.0) -- (602.6402750000001,160.0);
                \draw[species border] (709.3202750000002,202.5) |- (653.2002750000001,170.0) -- (653.2002750000001,160.0);
                \draw[species border] (602.6402750000001,225.0) |- (602.6402750000001,202.5) -| (653.2002750000001,202.5);
                \draw[species border] (528.7402750000001,225.0) -- (528.7402750000001,277.5) -- node[species label] {W} (602.6402750000001,277.5) -- (602.6402750000001,225.0);
                \draw[species border] (653.2002750000001,202.5) -- (653.2002750000001,277.5) -- node[species label] {T} (709.3202750000002,277.5) -- (709.3202750000002,202.5);
                % gene branches
                \draw[branch] (339.18027500000005,60.0) |- (483.86027500000006,38.75) -| (627.9202750000001,160.0);
                \draw[loss] (511.7802750000001,50.0) -- ++(20, 0);
                \draw[branch] (511.7802750000001,50.0) -| (292.51027500000004,60.0);
                \draw[branch] (497.82027500000004,16.25) -- (497.82027500000004,0);
                \draw[branch] (511.7802750000001,50.0) |- (497.82027500000004,16.25) -| (483.86027500000006,38.75);
                \draw[branch] (166.96,137.5) |- (292.51027500000004,98.75) -| (390.520275,247.5);
                \draw[branch] (292.51027500000004,76.25) -- (292.51027500000004,60.0);
                \draw[branch] (292.51027500000004,98.75) |- (292.51027500000004,76.25);
                \draw[branch] (146.96,137.5) |- (329.18027500000005,116.25) -| (411.910275,247.5);
                \draw[loss] (349.18027500000005,127.5) -- ++(20, 0);
                \draw[branch] (349.18027500000005,127.5) -| (111.68,137.5);
                \draw[branch] (339.18027500000005,76.25) -- (339.18027500000005,60.0);
                \draw[branch] (349.18027500000005,127.5) |- (339.18027500000005,76.25) -| (329.18027500000005,116.25);
                \draw[branch] (111.68,163.75) -- (111.68,137.5);
                \draw[branch] (44.87,247.5) |- (111.68,163.75) -| (193.49,202.5);
                \draw[branch] (146.96,181.25) -- (146.96,137.5);
                \draw[branch] (67.65,247.5) |- (146.96,181.25) -| (224.945,202.5);
                \draw[branch] (166.96,192.5) -- (166.96,137.5);
                \draw[loss] (166.96,192.5) -- ++(20, 0);
                \draw[branch] (166.96,192.5) -| (20.42,247.5);
                \draw[branch] (20.42,267.0) -- (20.42,247.5);
                \draw[branch] (44.870000000000005,267.0) -- (44.87,247.5);
                \draw[branch] (67.65,267.0) -- (67.65,247.5);
                \draw[branch] (193.49,267.0) -- (193.49,202.5);
                \draw[branch] (229.18,267.0) |- (238.9,241.25) -| (248.62,267.0);
                \draw[branch] (224.945,218.75) -- (224.945,202.5);
                \draw[branch] (210.99,267.0) |- (224.945,218.75) -| (238.9,241.25);
                \draw[branch] (390.520275,267.0) -- (390.520275,247.5);
                \draw[branch] (411.910275,267.0) -- (411.910275,247.5);
                \draw[branch] (431.910275,267.0) -- (431.910275,247.5);
                \draw[branch] (627.9202750000001,186.25) -- (627.9202750000001,160.0);
                \draw[branch] (583.6102750000001,225.0) |- (627.9202750000001,186.25) -| (680.0102750000001,202.5);
                \draw[branch] (552.9802750000001,241.25) -- (552.9802750000001,225.0);
                \draw[branch] (544.0202750000001,267.0) |- (552.9802750000001,241.25) -| (561.9402750000002,267.0);
                \draw[branch] (583.6102750000001,267.0) -- (583.6102750000001,225.0);
                \draw[branch] (669.7302750000001,267.0) |- (680.0102750000001,241.25) -| (690.2902750000002,267.0);
                \draw[branch] (680.0102750000001,218.75) -- (680.0102750000001,202.5);
                \draw[branch] (680.0102750000001,241.25) |- (680.0102750000001,218.75);
                % gene transfers
                \draw[transfer branch] (292.51027500000004,76.25) to[bend left=35] (552.9802750000001,225.0);
                \draw[transfer branch] (680.0102750000001,218.75) to[bend right=35] (431.910275,247.5);
                % events
                \node[speciation] at (483.86027500000006,38.75) {abcdefg};
                \node[duplication] at (497.82027500000004,16.25) {abcdefg};
                \node[speciation] at (292.51027500000004,98.75) {abcd};
                \node[horizontal gene transfer] at (292.51027500000004,76.25) {abcd};
                \node[speciation] at (329.18027500000005,116.25) {cdef};
                \node[duplication] at (339.18027500000005,76.25) {abcdefg};
                \node[speciation] at (111.68,163.75) {defg};
                \node[speciation] at (146.96,181.25) {cdef};
                \node[extant gene={abcd}] at (20.42,259.0) {};
                \node[extant gene={defg}] at (44.870000000000005,259.0) {};
                \node[extant gene={cdef}] at (67.65,259.0) {};
                \node[extant gene={def}] at (193.49,259.0) {};
                \node[extant gene={cef}] at (210.99,259.0) {};
                \node[extant gene={cde}] at (229.18,259.0) {};
                \node[extant gene={cde}] at (248.62,259.0) {};
                \node[duplication] at (238.9,241.25) {cde};
                \node[duplication] at (224.945,218.75) {cdef};
                \node[extant gene={abcd}] at (390.520275,259.0) {};
                \node[extant gene={cef}] at (411.910275,259.0) {};
                \node[extant gene={defg}] at (431.910275,259.0) {};
                \node[speciation] at (627.9202750000001,186.25) {defg};
                \node[extant gene={ab}] at (544.0202750000001,259.0) {};
                \node[extant gene={abc}] at (561.9402750000002,259.0) {};
                \node[duplication] at (552.9802750000001,241.25) {abc};
                \node[extant gene={defg}] at (583.6102750000001,259.0) {};
                \node[extant gene={def}] at (669.7302750000001,259.0) {};
                \node[extant gene={defg}] at (690.2902750000002,259.0) {};
                \node[duplication] at (680.0102750000001,241.25) {defg};
                \node[horizontal gene transfer] at (680.0102750000001,218.75) {defg};
                """
            ).lstrip(),
        )
