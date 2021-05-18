import unittest
import textwrap
from ete3 import PhyloTree
from .tools import (
    get_species_name,
    reconcile_leaves,
    parse_labeling,
    parse_reconciliation,
)
from .draw import compute_layout, render_to_tikz


class TestReconciliationDraw(unittest.TestCase):
    def assertRender(
        self, gene_text, species_text, rec_text, labeling_text, expect
    ):
        gene_tree = PhyloTree(
            gene_text, sp_naming_function=get_species_name, format=1
        )
        species_tree = PhyloTree(species_text, format=1)
        labeling = parse_labeling(gene_tree, labeling_text)
        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            **parse_reconciliation(gene_tree, species_tree, rec_text),
        }
        layout_info = compute_layout(gene_tree, species_tree, rec, labeling)
        out = render_to_tikz(species_tree, rec, layout_info)
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
                        minimum width={3pt},
                        minimum height={3pt},
                        label={[font={\strut}, inner xsep=0pt]below:#1},
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
        self.maxDiff = None
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
                \draw[species border] (0,85.0) -- (0,126.86595) -- node[species label] {X} (29.762999999999998,126.86595) -- (29.762999999999998,85.0);
                \draw[species border] (62.263,85.0) -- (62.263,126.86595) -- node[species label] {Y} (92.026,126.86595) -- (92.026,85.0);
                \draw[species border] (124.526,85.0) -- (124.526,126.86595) -- node[species label] {Z} (153.449,126.86595) -- (153.449,85.0);
                % gene branches
                \draw[branch] (108.276,26.25) -- (108.276,0);
                \draw[branch] (46.013,42.5) |- (108.276,26.25) -| (138.9875,85.0);
                \draw[branch] (46.013,68.75) -- (46.013,42.5);
                \draw[branch] (14.881499999999999,85.0) |- (46.013,68.75) -| (77.1445,85.0);
                \draw[branch] (14.881499999999999,105.932975) -- (14.881499999999999,85.0);
                \draw[branch] (77.1445,105.932975) -- (77.1445,85.0);
                \draw[branch] (138.9875,105.932975) -- (138.9875,85.0);
                % gene transfers
                % events
                \node[speciation] at (108.276,26.25) {};
                \node[speciation] at (46.013,68.75) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,105.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (77.1445,105.932975) {};
                \node[extant gene={z\textsubscript{1}}] at (138.9875,105.932975) {};
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
                \draw[species border] (59.288999999999994,107.5) |- (186.07799999999997,40.0) -- (186.07799999999997,0);
                \draw[species border] (310.347,185.0) |- (253.57799999999997,40.0) -- (253.57799999999997,0);
                \draw[species border] (126.78899999999999,107.5) |- (186.07799999999997,107.5) -| (253.57799999999997,185.0);
                \draw[species border] (0,185.0) |- (59.288999999999994,117.5) -- (59.288999999999994,107.5);
                \draw[species border] (186.07799999999997,185.0) |- (126.78899999999999,117.5) -- (126.78899999999999,107.5);
                \draw[species border] (59.288999999999994,185.0) |- (59.288999999999994,185.0) -| (126.78899999999999,185.0);
                \draw[species border] (0,185.0) -- (0,226.86595) -- node[species label] {X} (59.288999999999994,226.86595) -- (59.288999999999994,185.0);
                \draw[species border] (126.78899999999999,185.0) -- (126.78899999999999,226.86595) -- node[species label] {Y} (186.07799999999997,226.86595) -- (186.07799999999997,185.0);
                \draw[species border] (253.57799999999997,185.0) -- (253.57799999999997,226.86595) -- node[species label] {Z} (310.347,226.86595) -- (310.347,185.0);
                % gene branches
                \draw[branch] (75.53899999999999,107.5) |- (202.32799999999997,56.25) -| (268.0395,185.0);
                \draw[branch] (93.03899999999999,107.5) |- (219.82799999999997,73.75) -| (281.9625,185.0);
                \draw[branch] (110.53899999999999,107.5) |- (237.32799999999997,91.25) -| (295.8855,185.0);
                \draw[branch] (219.82799999999997,73.75) |- (228.57799999999997,33.75) -| (237.32799999999997,91.25);
                \draw[branch] (215.45299999999997,16.25) -- (215.45299999999997,0);
                \draw[branch] (202.32799999999997,56.25) |- (215.45299999999997,16.25) -| (228.57799999999997,33.75);
                \draw[branch] (75.53899999999999,133.75) -- (75.53899999999999,107.5);
                \draw[branch] (14.881499999999999,185.0) |- (75.53899999999999,133.75) -| (141.67049999999998,185.0);
                \draw[branch] (93.03899999999999,151.25) -- (93.03899999999999,107.5);
                \draw[branch] (29.6445,185.0) |- (93.03899999999999,151.25) -| (156.43349999999998,185.0);
                \draw[branch] (110.53899999999999,168.75) -- (110.53899999999999,107.5);
                \draw[branch] (44.4075,185.0) |- (110.53899999999999,168.75) -| (171.1965,185.0);
                \draw[branch] (14.881499999999999,205.932975) -- (14.881499999999999,185.0);
                \draw[branch] (29.644499999999997,205.932975) -- (29.6445,185.0);
                \draw[branch] (44.4075,205.932975) -- (44.4075,185.0);
                \draw[branch] (141.67049999999998,205.932975) -- (141.67049999999998,185.0);
                \draw[branch] (156.43349999999998,205.932975) -- (156.43349999999998,185.0);
                \draw[branch] (171.1965,205.932975) -- (171.1965,185.0);
                \draw[branch] (268.0395,205.932975) -- (268.0395,185.0);
                \draw[branch] (281.9625,205.932975) -- (281.9625,185.0);
                \draw[branch] (295.8855,205.932975) -- (295.8855,185.0);
                % gene transfers
                % events
                \node[speciation] at (202.32799999999997,56.25) {};
                \node[speciation] at (219.82799999999997,73.75) {};
                \node[speciation] at (237.32799999999997,91.25) {};
                \node[duplication] at (228.57799999999997,33.75) {};
                \node[duplication] at (215.45299999999997,16.25) {};
                \node[speciation] at (75.53899999999999,133.75) {};
                \node[speciation] at (93.03899999999999,151.25) {};
                \node[speciation] at (110.53899999999999,168.75) {};
                \node[extant gene={x\textsubscript{3}}] at (14.881499999999999,205.932975) {};
                \node[extant gene={x\textsubscript{1}}] at (29.644499999999997,205.932975) {};
                \node[extant gene={x\textsubscript{2}}] at (44.4075,205.932975) {};
                \node[extant gene={y\textsubscript{3}}] at (141.67049999999998,205.932975) {};
                \node[extant gene={y\textsubscript{2}}] at (156.43349999999998,205.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (171.1965,205.932975) {};
                \node[extant gene={z\textsubscript{3}}] at (268.0395,205.932975) {};
                \node[extant gene={z\textsubscript{2}}] at (281.9625,205.932975) {};
                \node[extant gene={z\textsubscript{1}}] at (295.8855,205.932975) {};
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
                \draw[species border] (0,102.5) -- (0,144.36595) -- node[species label] {X} (29.762999999999998,144.36595) -- (29.762999999999998,102.5);
                \draw[species border] (62.263,124.36595) |- (82.263,52.5) -- (82.263,42.5);
                \draw[species border] (151.186,72.5) |- (102.263,52.5) -- (102.263,42.5);
                \draw[species border] (82.263,124.36595) |- (82.263,72.5) -| (131.186,72.5);
                \draw[species border] (62.263,124.36595) -- (62.263,144.36595) -- node[species label] {Y} (82.263,144.36595) -- (82.263,124.36595);
                \draw[species border] (102.263,102.5) |- (131.186,82.5) -- (131.186,72.5);
                \draw[species border] (171.186,124.36595) |- (151.186,82.5) -- (151.186,72.5);
                \draw[species border] (131.186,102.5) |- (131.186,102.5) -| (151.186,124.36595);
                \draw[species border] (102.263,102.5) -- (102.263,144.36595) -- node[species label] {Z} (131.186,144.36595) -- (131.186,102.5);
                \draw[species border] (151.186,124.36595) -- (151.186,144.36595) -- node[species label] {W} (171.186,144.36595) -- (171.186,124.36595);
                % gene branches
                \draw[branch] (46.013,26.25) -- (46.013,0);
                \draw[branch] (14.881499999999999,102.5) |- (46.013,26.25) -| (92.263,42.5);
                \draw[branch] (14.881499999999999,123.432975) -- (14.881499999999999,102.5);
                \draw[branch] (92.263,62.5) -- (92.263,42.5);
                \draw[loss] (92.263,62.5) -- ++(-20, 0);
                \draw[branch] (92.263,62.5) -| (141.186,72.5);
                \draw[branch] (141.186,92.5) -- (141.186,72.5);
                \draw[loss] (141.186,92.5) -- ++(20, 0);
                \draw[branch] (141.186,92.5) -| (116.7245,102.5);
                \draw[branch] (116.7245,123.432975) -- (116.7245,102.5);
                % gene transfers
                % events
                \node[speciation] at (46.013,26.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,123.432975) {};
                \node[extant gene={z\textsubscript{1}}] at (116.7245,123.432975) {};
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
                \draw[species border] (69.763,94.36595) |- (49.763,52.5) -- (49.763,42.5);
                \draw[species border] (29.762999999999998,72.5) |- (29.762999999999998,72.5) -| (49.763,94.36595);
                \draw[species border] (0,72.5) -- (0,114.36595) -- node[species label] {X} (29.762999999999998,114.36595) -- (29.762999999999998,72.5);
                \draw[species border] (49.763,94.36595) -- (49.763,114.36595) -- node[species label] {Y} (69.763,114.36595) -- (69.763,94.36595);
                \draw[species border] (102.263,94.36595) |- (122.263,52.5) -- (122.263,42.5);
                \draw[species border] (173.966,72.5) |- (142.263,52.5) -- (142.263,42.5);
                \draw[species border] (122.263,94.36595) |- (122.263,72.5) -| (142.263,72.5);
                \draw[species border] (102.263,94.36595) -- (102.263,114.36595) -- node[species label] {Z} (122.263,114.36595) -- (122.263,94.36595);
                \draw[species border] (142.263,72.5) -- (142.263,114.36595) -- node[species label] {W} (173.966,114.36595) -- (173.966,72.5);
                % gene branches
                \draw[branch] (86.013,26.25) -- (86.013,0);
                \draw[branch] (39.763,42.5) |- (86.013,26.25) -| (132.263,42.5);
                \draw[branch] (39.763,62.5) -- (39.763,42.5);
                \draw[loss] (39.763,62.5) -- ++(20, 0);
                \draw[branch] (39.763,62.5) -| (14.881499999999999,72.5);
                \draw[branch] (14.881499999999999,93.432975) -- (14.881499999999999,72.5);
                \draw[branch] (132.263,62.5) -- (132.263,42.5);
                \draw[loss] (132.263,62.5) -- ++(-20, 0);
                \draw[branch] (132.263,62.5) -| (158.1145,72.5);
                \draw[branch] (158.1145,93.432975) -- (158.1145,72.5);
                % gene transfers
                % events
                \node[speciation] at (86.013,26.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,93.432975) {};
                \node[extant gene={w\textsubscript{1}}] at (158.1145,93.432975) {};
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
                \draw[species border] (44.525999999999996,95.0) |- (139.052,40.0) -- (139.052,0);
                \draw[species border] (250.821,155.0) |- (194.052,40.0) -- (194.052,0);
                \draw[species border] (94.526,95.0) |- (139.052,95.0) -| (194.052,155.0);
                \draw[species border] (0,155.0) |- (44.525999999999996,105.0) -- (44.525999999999996,95.0);
                \draw[species border] (139.052,155.0) |- (94.526,105.0) -- (94.526,95.0);
                \draw[species border] (44.525999999999996,155.0) |- (44.525999999999996,155.0) -| (94.526,155.0);
                \draw[species border] (0,155.0) -- (0,196.86595) -- node[species label] {X} (44.525999999999996,196.86595) -- (44.525999999999996,155.0);
                \draw[species border] (94.526,155.0) -- (94.526,196.86595) -- node[species label] {Y} (139.052,196.86595) -- (139.052,155.0);
                \draw[species border] (194.052,155.0) -- (194.052,196.86595) -- node[species label] {Z} (250.821,196.86595) -- (250.821,155.0);
                % gene branches
                \draw[branch] (60.775999999999996,95.0) |- (155.302,56.25) -| (222.4365,155.0);
                \draw[branch] (78.276,95.0) |- (172.802,73.75) -| (236.3595,155.0);
                \draw[branch] (155.302,56.25) |- (164.052,33.75) -| (172.802,73.75);
                \draw[loss] (184.052,85.0) -- ++(-20, 0);
                \draw[branch] (184.052,85.0) -| (208.5135,155.0);
                \draw[branch] (174.052,16.25) -- (174.052,0);
                \draw[branch] (184.052,85.0) |- (174.052,16.25) -| (164.052,33.75);
                \draw[branch] (60.775999999999996,121.25) -- (60.775999999999996,95.0);
                \draw[branch] (14.881499999999999,155.0) |- (60.775999999999996,121.25) -| (109.4075,155.0);
                \draw[branch] (78.276,138.75) -- (78.276,95.0);
                \draw[branch] (29.6445,155.0) |- (78.276,138.75) -| (124.1705,155.0);
                \draw[branch] (14.881499999999999,175.932975) -- (14.881499999999999,155.0);
                \draw[branch] (29.644499999999997,175.932975) -- (29.6445,155.0);
                \draw[branch] (109.4075,175.932975) -- (109.4075,155.0);
                \draw[branch] (124.17049999999999,175.932975) -- (124.1705,155.0);
                \draw[branch] (208.5135,175.932975) -- (208.5135,155.0);
                \draw[branch] (222.4365,175.932975) -- (222.4365,155.0);
                \draw[branch] (236.3595,175.932975) -- (236.3595,155.0);
                % gene transfers
                % events
                \node[speciation] at (155.302,56.25) {};
                \node[speciation] at (172.802,73.75) {};
                \node[duplication] at (164.052,33.75) {};
                \node[duplication] at (174.052,16.25) {};
                \node[speciation] at (60.775999999999996,121.25) {};
                \node[speciation] at (78.276,138.75) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,175.932975) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,175.932975) {};
                \node[extant gene={y\textsubscript{2}}] at (109.4075,175.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (124.17049999999999,175.932975) {};
                \node[extant gene={z\textsubscript{3}}] at (208.5135,175.932975) {};
                \node[extant gene={z\textsubscript{2}}] at (222.4365,175.932975) {};
                \node[extant gene={z\textsubscript{1}}] at (236.3595,175.932975) {};
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
                \draw[species border] (0,125.0) |- (74.052,40.0) -- (74.052,0);
                \draw[species border] (233.10400000000004,125.0) |- (159.05200000000002,40.0) -- (159.05200000000002,0);
                \draw[species border] (74.052,125.0) |- (74.052,125.0) -| (159.05200000000002,125.0);
                \draw[species border] (0,125.0) -- (0,166.86595) -- node[species label] {X} (74.052,166.86595) -- (74.052,125.0);
                \draw[species border] (159.05200000000002,125.0) -- (159.05200000000002,166.86595) -- node[species label] {Y} (233.10400000000004,166.86595) -- (233.10400000000004,125.0);
                % gene branches
                \draw[branch] (14.881499999999999,125.0) |- (90.302,56.25) -| (173.9335,125.0);
                \draw[branch] (29.6445,125.0) |- (107.802,73.75) -| (188.69650000000001,125.0);
                \draw[branch] (90.302,56.25) |- (99.052,33.75) -| (107.802,73.75);
                \draw[branch] (44.4075,125.0) |- (125.302,91.25) -| (203.45950000000002,125.0);
                \draw[branch] (59.170500000000004,125.0) |- (142.80200000000002,108.75) -| (218.22250000000003,125.0);
                \draw[branch] (125.302,91.25) |- (134.05200000000002,33.75) -| (142.80200000000002,108.75);
                \draw[branch] (116.552,16.25) -- (116.552,0);
                \draw[branch] (99.052,33.75) |- (116.552,16.25) -| (134.05200000000002,33.75);
                \draw[branch] (14.881499999999999,145.932975) -- (14.881499999999999,125.0);
                \draw[branch] (29.644499999999997,145.932975) -- (29.6445,125.0);
                \draw[branch] (44.4075,145.932975) -- (44.4075,125.0);
                \draw[branch] (59.170500000000004,145.932975) -- (59.170500000000004,125.0);
                \draw[branch] (173.9335,145.932975) -- (173.9335,125.0);
                \draw[branch] (188.69650000000001,145.932975) -- (188.69650000000001,125.0);
                \draw[branch] (203.45950000000002,145.932975) -- (203.45950000000002,125.0);
                \draw[branch] (218.22250000000003,145.932975) -- (218.22250000000003,125.0);
                % gene transfers
                % events
                \node[speciation] at (90.302,56.25) {};
                \node[speciation] at (107.802,73.75) {};
                \node[duplication] at (99.052,33.75) {};
                \node[speciation] at (125.302,91.25) {};
                \node[speciation] at (142.80200000000002,108.75) {};
                \node[duplication] at (134.05200000000002,33.75) {};
                \node[duplication] at (116.552,16.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,145.932975) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,145.932975) {};
                \node[extant gene={x\textsubscript{3}}] at (44.4075,145.932975) {};
                \node[extant gene={x\textsubscript{4}}] at (59.170500000000004,145.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (173.9335,145.932975) {};
                \node[extant gene={y\textsubscript{2}}] at (188.69650000000001,145.932975) {};
                \node[extant gene={y\textsubscript{3}}] at (203.45950000000002,145.932975) {};
                \node[extant gene={y\textsubscript{4}}] at (218.22250000000003,145.932975) {};
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
                \draw[species border] (29.762999999999998,53.0) |- (84.526,22.5) -- (84.526,0);
                \draw[species border] (135.026,109.86595) |- (115.026,22.5) -- (115.026,0);
                \draw[species border] (54.763,53.0) |- (84.526,47.5) -| (115.026,109.86595);
                \draw[species border] (0,88.0) |- (29.762999999999998,63.0) -- (29.762999999999998,53.0);
                \draw[species border] (84.526,88.0) |- (54.763,63.0) -- (54.763,53.0);
                \draw[species border] (29.762999999999998,88.0) |- (29.762999999999998,88.0) -| (54.763,88.0);
                \draw[species border] (0,88.0) -- (0,129.86595) -- node[species label] {X} (29.762999999999998,129.86595) -- (29.762999999999998,88.0);
                \draw[species border] (54.763,88.0) -- (54.763,129.86595) -- node[species label] {Y} (84.526,129.86595) -- (84.526,88.0);
                \draw[species border] (115.026,109.86595) -- (115.026,129.86595) -- node[species label] {Z} (135.026,129.86595) -- (135.026,109.86595);
                % gene branches
                \draw[loss] (97.276,32.5) -- ++(20, 0);
                \draw[branch] (97.276,32.5) -| (39.763,53.0);
                \draw[loss] (102.276,37.5) -- ++(20, 0);
                \draw[branch] (102.276,37.5) -| (44.763,53.0);
                \draw[branch] (99.776,16.25) -- (99.776,0);
                \draw[branch] (97.276,32.5) |- (99.776,16.25) -| (102.276,37.5);
                \draw[branch] (39.763,73.0) -- (39.763,53.0);
                \draw[loss] (39.763,73.0) -- ++(-20, 0);
                \draw[branch] (39.763,73.0) -| (69.6445,88.0);
                \draw[branch] (44.763,78.0) -- (44.763,53.0);
                \draw[loss] (44.763,78.0) -- ++(20, 0);
                \draw[branch] (44.763,78.0) -| (14.881499999999999,88.0);
                \draw[branch] (14.881499999999999,108.932975) -- (14.881499999999999,88.0);
                \draw[branch] (69.6445,108.932975) -- (69.6445,88.0);
                % gene transfers
                % events
                \node[duplication] at (99.776,16.25) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,108.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (69.6445,108.932975) {};
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
                \draw[species border] (44.525999999999996,59.999999999999986) |- (139.052275,22.5) -- (139.052275,0);
                \draw[species border] (219.398275,132.500275) |- (176.552275,22.5) -- (176.552275,0);
                \draw[species border] (94.526275,59.999999999999986) |- (139.052275,60.0) -| (176.552275,132.500275);
                \draw[species border] (0,132.500275) |- (44.525999999999996,82.49999999999999) -- (44.525999999999996,59.999999999999986);
                \draw[species border] (139.052275,132.500275) |- (94.526275,82.49999999999999) -- (94.526275,59.999999999999986);
                \draw[species border] (44.525999999999996,132.500275) |- (44.525999999999996,132.5) -| (94.526275,132.500275);
                \draw[species border] (0,132.500275) -- (0,174.366225) -- node[species label] {X} (44.525999999999996,174.366225) -- (44.525999999999996,132.500275);
                \draw[species border] (94.526275,132.500275) -- (94.526275,174.366225) -- node[species label] {Y} (139.052275,174.366225) -- (139.052275,132.500275);
                \draw[species border] (176.552275,132.500275) -- (176.552275,174.366225) -- node[species label] {Z} (219.398275,174.366225) -- (219.398275,132.500275);
                % gene branches
                \draw[branch] (78.276,59.999999999999986) |- (155.302275,38.75) -| (204.936775,132.500275);
                \draw[loss] (166.552275,50.0) -- ++(20, 0);
                \draw[branch] (166.552275,50.0) -| (60.775999999999996,59.999999999999986);
                \draw[branch] (160.927275,16.25) -- (160.927275,0);
                \draw[branch] (166.552275,50.0) |- (160.927275,16.25) -| (155.302275,38.75);
                \draw[branch] (60.775999999999996,98.74999999999999) -- (60.775999999999996,59.999999999999986);
                \draw[branch] (14.881499999999999,132.500275) |- (60.775999999999996,98.74999999999999) -| (109.407775,132.500275);
                \draw[branch] (29.6445,132.500275) |- (78.276,116.24999999999999) -| (124.17077499999999,132.500275);
                \draw[branch] (78.276,76.24999999999999) -- (78.276,59.999999999999986);
                \draw[branch] (78.276,116.24999999999999) |- (78.276,76.24999999999999);
                \draw[branch] (14.881499999999999,153.43325) -- (14.881499999999999,132.500275);
                \draw[branch] (29.644499999999997,153.43325) -- (29.6445,132.500275);
                \draw[branch] (109.407775,153.43325) -- (109.407775,132.500275);
                \draw[branch] (124.17077499999999,153.43325) -- (124.17077499999999,132.500275);
                \draw[branch] (191.013775,153.43325) -- (191.013775,132.500275);
                \draw[branch] (204.936775,153.43325) -- (204.936775,132.500275);
                % gene transfers
                \draw[transfer branch] (78.276,76.24999999999999) to[bend left=35] (191.013775,132.500275);
                % events
                \node[speciation] at (155.302275,38.75) {};
                \node[duplication] at (160.927275,16.25) {};
                \node[speciation] at (60.775999999999996,98.74999999999999) {};
                \node[speciation] at (78.276,116.24999999999999) {};
                \node[horizontal gene transfer] at (78.276,76.24999999999999) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,153.43325) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,153.43325) {};
                \node[extant gene={y\textsubscript{1}}] at (109.407775,153.43325) {};
                \node[extant gene={y\textsubscript{2}}] at (124.17077499999999,153.43325) {};
                \node[extant gene={z\textsubscript{1}}] at (191.013775,153.43325) {};
                \node[extant gene={z\textsubscript{2}}] at (204.936775,153.43325) {};
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
                \draw[species border] (0,85.00055) |- (29.762999999999998,10) -- (29.762999999999998,0);
                \draw[species border] (112.02655,30.0) |- (49.763,10) -- (49.763,0);
                \draw[species border] (29.762999999999998,85.00055) |- (29.762999999999998,20) -| (79.526,30.0);
                \draw[species border] (0,85.00055) -- (0,126.8665) -- node[species label] {X} (29.762999999999998,126.8665) -- (29.762999999999998,85.00055);
                \draw[species border] (49.763,85.00055) |- (79.526,52.5) -- (79.526,30.0);
                \draw[species border] (140.94955,85.00055) |- (112.02655,52.5) -- (112.02655,30.0);
                \draw[species border] (79.526,85.00055) |- (79.526,85.0) -| (112.02655,85.00055);
                \draw[species border] (49.763,85.00055) -- (49.763,126.8665) -- node[species label] {Y} (79.526,126.8665) -- (79.526,85.00055);
                \draw[species border] (112.02655,85.00055) -- (112.02655,126.8665) -- node[species label] {Z} (140.94955,126.8665) -- (140.94955,85.00055);
                % gene branches
                \draw[branch] (14.881499999999999,105.933525) -- (14.881499999999999,85.00055);
                \draw[branch] (64.6445,85.00055) |- (95.776275,68.75) -| (126.48805,85.00055);
                \draw[branch] (95.776275,46.25) -- (95.776275,30.0);
                \draw[branch] (95.776275,68.75) |- (95.776275,46.25);
                \draw[branch] (64.6445,105.933525) -- (64.6445,85.00055);
                \draw[branch] (126.48805,105.933525) -- (126.48805,85.00055);
                % gene transfers
                \draw[transfer branch] (95.776275,46.25) to[bend right=35] (14.881499999999999,85.00055);
                % events
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,105.933525) {};
                \node[speciation] at (95.776275,68.75) {};
                \node[horizontal gene transfer] at (95.776275,46.25) {};
                \node[extant gene={y\textsubscript{1}}] at (64.6445,105.933525) {};
                \node[extant gene={z\textsubscript{1}}] at (126.48805,105.933525) {};
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
                \draw[species border] (435.219275,160.000275) |- (337.610275,22.5) -- (337.610275,0);
                \draw[species border] (243.341275,60.0) |- (300.110275,60.0) -| (402.719275,160.000275);
                \draw[species border] (59.288999999999994,137.500275) |- (188.341,82.5) -- (188.341,60.0);
                \draw[species border] (300.110275,237.500275) |- (243.341275,82.5) -- (243.341275,60.0);
                \draw[species border] (114.28899999999999,137.500275) |- (188.341,137.5) -| (243.341275,237.500275);
                \draw[species border] (0,237.500275) |- (59.288999999999994,147.500275) -- (59.288999999999994,137.500275);
                \draw[species border] (188.341,202.500275) |- (114.28899999999999,147.500275) -- (114.28899999999999,137.500275);
                \draw[species border] (59.288999999999994,237.500275) |- (59.288999999999994,202.500275) -| (114.28899999999999,202.500275);
                \draw[species border] (0,237.500275) -- (0,279.366225) -- node[species label] {X} (59.288999999999994,279.366225) -- (59.288999999999994,237.500275);
                \draw[species border] (114.28899999999999,202.500275) -- (114.28899999999999,279.366225) -- node[species label] {Y} (188.341,279.366225) -- (188.341,202.500275);
                \draw[species border] (243.341275,237.500275) -- (243.341275,279.366225) -- node[species label] {Z} (300.110275,279.366225) -- (300.110275,237.500275);
                \draw[species border] (337.610275,220.000275) |- (402.719275,170.000275) -- (402.719275,160.000275);
                \draw[species border] (476.96525499999996,202.500275) |- (435.219275,170.000275) -- (435.219275,160.000275);
                \draw[species border] (402.719275,220.000275) |- (402.719275,202.500275) -| (435.219275,202.500275);
                \draw[species border] (337.610275,220.000275) -- (337.610275,279.366225) -- node[species label] {W} (402.719275,279.366225) -- (402.719275,220.000275);
                \draw[species border] (435.219275,202.500275) -- (435.219275,279.366225) -- node[species label] {T} (476.96525499999996,279.366225) -- (476.96525499999996,202.500275);
                % gene branches
                \draw[branch] (227.716275,60.0) |- (316.360275,38.75) -| (418.969275,160.000275);
                \draw[loss] (327.610275,50.0) -- ++(20, 0);
                \draw[branch] (327.610275,50.0) -| (204.591275,60.0);
                \draw[branch] (321.985275,16.25) -- (321.985275,0);
                \draw[branch] (327.610275,50.0) |- (321.985275,16.25) -| (316.360275,38.75);
                \draw[branch] (104.28899999999999,137.500275) |- (204.591275,98.75) -| (257.802775,237.500275);
                \draw[branch] (204.591275,76.25) -- (204.591275,60.0);
                \draw[branch] (204.591275,98.75) |- (204.591275,76.25);
                \draw[branch] (93.03899999999999,137.500275) |- (222.091275,116.25) -| (271.725775,237.500275);
                \draw[loss] (233.341275,127.5) -- ++(20, 0);
                \draw[branch] (233.341275,127.5) -| (75.53899999999999,137.500275);
                \draw[branch] (227.716275,76.25) -- (227.716275,60.0);
                \draw[branch] (233.341275,127.5) |- (227.716275,76.25) -| (222.091275,116.25);
                \draw[branch] (75.53899999999999,163.750275) -- (75.53899999999999,137.500275);
                \draw[branch] (29.6445,237.500275) |- (75.53899999999999,163.750275) -| (129.17049999999998,202.500275);
                \draw[branch] (93.03899999999999,181.250275) -- (93.03899999999999,137.500275);
                \draw[branch] (44.4075,237.500275) |- (93.03899999999999,181.250275) -| (155.00574999999998,202.500275);
                \draw[branch] (104.28899999999999,192.500275) -- (104.28899999999999,137.500275);
                \draw[loss] (104.28899999999999,192.500275) -- ++(20, 0);
                \draw[branch] (104.28899999999999,192.500275) -| (14.881499999999999,237.500275);
                \draw[branch] (14.881499999999999,258.43325) -- (14.881499999999999,237.500275);
                \draw[branch] (29.644499999999997,258.43325) -- (29.6445,237.500275);
                \draw[branch] (44.4075,258.43325) -- (44.4075,237.500275);
                \draw[branch] (129.17049999999998,258.43325) -- (129.17049999999998,202.500275);
                \draw[branch] (158.6965,258.43325) |- (166.07799999999997,236.250275) -| (173.4595,258.43325);
                \draw[branch] (155.00574999999998,218.750275) -- (155.00574999999998,202.500275);
                \draw[branch] (143.93349999999998,258.43325) |- (155.00574999999998,218.750275) -| (166.07799999999997,236.250275);
                \draw[branch] (257.802775,258.43325) -- (257.802775,237.500275);
                \draw[branch] (271.725775,258.43325) -- (271.725775,237.500275);
                \draw[branch] (285.648775,258.43325) -- (285.648775,237.500275);
                \draw[branch] (418.969275,186.250275) -- (418.969275,160.000275);
                \draw[branch] (386.867775,220.000275) |- (418.969275,186.250275) -| (456.092265,202.500275);
                \draw[branch] (361.813275,236.250275) -- (361.813275,220.000275);
                \draw[branch] (353.461775,258.43325) |- (361.813275,236.250275) -| (370.164775,258.43325);
                \draw[branch] (386.867775,258.43325) -- (386.867775,220.000275);
                \draw[branch] (449.40576999999996,258.43325) |- (456.092265,236.250275) -| (462.77876,258.43325);
                \draw[branch] (456.092265,218.750275) -- (456.092265,202.500275);
                \draw[branch] (456.092265,236.250275) |- (456.092265,218.750275);
                % gene transfers
                \draw[transfer branch] (204.591275,76.25) to[bend left=35] (361.813275,220.000275);
                \draw[transfer branch] (456.092265,218.750275) to[bend right=35] (285.648775,237.500275);
                % events
                \node[speciation] at (316.360275,38.75) {};
                \node[duplication] at (321.985275,16.25) {};
                \node[speciation] at (204.591275,98.75) {};
                \node[horizontal gene transfer] at (204.591275,76.25) {};
                \node[speciation] at (222.091275,116.25) {};
                \node[duplication] at (227.716275,76.25) {};
                \node[speciation] at (75.53899999999999,163.750275) {};
                \node[speciation] at (93.03899999999999,181.250275) {};
                \node[extant gene={x\textsubscript{1}}] at (14.881499999999999,258.43325) {};
                \node[extant gene={x\textsubscript{2}}] at (29.644499999999997,258.43325) {};
                \node[extant gene={x\textsubscript{3}}] at (44.4075,258.43325) {};
                \node[extant gene={y\textsubscript{4}}] at (129.17049999999998,258.43325) {};
                \node[extant gene={y\textsubscript{1}}] at (143.93349999999998,258.43325) {};
                \node[extant gene={y\textsubscript{2}}] at (158.6965,258.43325) {};
                \node[extant gene={y\textsubscript{3}}] at (173.4595,258.43325) {};
                \node[duplication] at (166.07799999999997,236.250275) {};
                \node[duplication] at (155.00574999999998,218.750275) {};
                \node[extant gene={z\textsubscript{1}}] at (257.802775,258.43325) {};
                \node[extant gene={z\textsubscript{2}}] at (271.725775,258.43325) {};
                \node[extant gene={z\textsubscript{3}}] at (285.648775,258.43325) {};
                \node[speciation] at (418.969275,186.250275) {};
                \node[extant gene={w\textsubscript{1}}] at (353.461775,258.43325) {};
                \node[extant gene={w\textsubscript{2}}] at (370.164775,258.43325) {};
                \node[duplication] at (361.813275,236.250275) {};
                \node[extant gene={w\textsubscript{3}}] at (386.867775,258.43325) {};
                \node[extant gene={t\textsubscript{1}}] at (449.40576999999996,258.43325) {};
                \node[extant gene={t\textsubscript{2}}] at (462.77876,258.43325) {};
                \node[duplication] at (456.092265,236.250275) {};
                \node[horizontal gene transfer] at (456.092265,218.750275) {};
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
                \draw[species border] (265.84000000000003,100.30000000000001) |- (450.94027500000004,22.5) -- (450.94027500000004,0);
                \draw[species border] (653.2002750000001,267.06027500000005) |- (528.7402750000001,22.5) -- (528.7402750000001,0);
                \draw[species border] (370.100275,100.30000000000001) |- (450.94027500000004,60.0) -| (602.6402750000001,267.06027500000005);
                \draw[species border] (86.4,227.06027500000005) |- (265.84000000000003,122.80000000000001) -- (265.84000000000003,100.30000000000001);
                \draw[species border] (450.94027500000004,362.62027500000005) |- (370.100275,122.80000000000001) -- (370.100275,100.30000000000001);
                \draw[species border] (176.96,227.06027500000005) |- (265.84000000000003,177.8) -| (370.100275,362.62027500000005);
                \draw[species border] (0,362.62027500000005) |- (86.4,237.06027500000005) -- (86.4,227.06027500000005);
                \draw[species border] (265.84000000000003,327.62027500000005) |- (176.96,237.06027500000005) -- (176.96,227.06027500000005);
                \draw[species border] (86.4,362.62027500000005) |- (86.4,292.06027500000005) -| (176.96,327.62027500000005);
                \draw[species border] (0,362.62027500000005) -- (0,404.48622500000005) -- node[species label] {X} (86.4,404.48622500000005) -- (86.4,362.62027500000005);
                \draw[species border] (176.96,327.62027500000005) -- (176.96,404.48622500000005) -- node[species label] {Y} (265.84000000000003,404.48622500000005) -- (265.84000000000003,327.62027500000005);
                \draw[species border] (370.100275,362.62027500000005) -- (370.100275,404.48622500000005) -- node[species label] {Z} (450.94027500000004,404.48622500000005) -- (450.94027500000004,362.62027500000005);
                \draw[species border] (528.7402750000001,345.12027500000005) |- (602.6402750000001,277.06027500000005) -- (602.6402750000001,267.06027500000005);
                \draw[species border] (709.3202750000002,327.62027500000005) |- (653.2002750000001,277.06027500000005) -- (653.2002750000001,267.06027500000005);
                \draw[species border] (602.6402750000001,345.12027500000005) |- (602.6402750000001,309.56027500000005) -| (653.2002750000001,327.62027500000005);
                \draw[species border] (528.7402750000001,345.12027500000005) -- (528.7402750000001,404.48622500000005) -- node[species label] {W} (602.6402750000001,404.48622500000005) -- (602.6402750000001,345.12027500000005);
                \draw[species border] (653.2002750000001,327.62027500000005) -- (653.2002750000001,404.48622500000005) -- node[species label] {T} (709.3202750000002,404.48622500000005) -- (709.3202750000002,327.62027500000005);
                % gene branches
                \draw[branch] (339.18027500000005,100.30000000000001) |- (483.86027500000006,38.75) -| (627.9202750000001,267.06027500000005);
                \draw[loss] (511.7802750000001,50.0) -- ++(20, 0);
                \draw[branch] (511.7802750000001,50.0) -| (292.51027500000004,100.30000000000001);
                \draw[branch] (497.82027500000004,16.25) -- (497.82027500000004,0);
                \draw[branch] (511.7802750000001,50.0) |- (497.82027500000004,16.25) -| (483.86027500000006,38.75);
                \draw[branch] (166.96,227.06027500000005) |- (292.51027500000004,139.05) -| (390.520275,362.62027500000005);
                \draw[branch] (292.51027500000004,116.55000000000001) -- (292.51027500000004,100.30000000000001);
                \draw[branch] (292.51027500000004,139.05) |- (292.51027500000004,116.55000000000001);
                \draw[branch] (146.96,227.06027500000005) |- (329.18027500000005,156.55) -| (411.910275,362.62027500000005);
                \draw[loss] (349.18027500000005,167.8) -- ++(20, 0);
                \draw[branch] (349.18027500000005,167.8) -| (111.68,227.06027500000005);
                \draw[branch] (339.18027500000005,116.55000000000001) -- (339.18027500000005,100.30000000000001);
                \draw[branch] (349.18027500000005,167.8) |- (339.18027500000005,116.55000000000001) -| (329.18027500000005,156.55);
                \draw[branch] (111.68,253.31027500000005) -- (111.68,227.06027500000005);
                \draw[branch] (44.87,362.62027500000005) |- (111.68,253.31027500000005) -| (193.49,327.62027500000005);
                \draw[branch] (146.96,270.81027500000005) -- (146.96,227.06027500000005);
                \draw[branch] (67.65,362.62027500000005) |- (146.96,270.81027500000005) -| (224.945,327.62027500000005);
                \draw[branch] (166.96,282.06027500000005) -- (166.96,227.06027500000005);
                \draw[loss] (166.96,282.06027500000005) -- ++(20, 0);
                \draw[branch] (166.96,282.06027500000005) -| (20.42,362.62027500000005);
                \draw[branch] (20.42,383.55325000000005) -- (20.42,362.62027500000005);
                \draw[branch] (44.870000000000005,383.55325000000005) -- (44.87,362.62027500000005);
                \draw[branch] (67.65,383.55325000000005) -- (67.65,362.62027500000005);
                \draw[branch] (193.49,383.55325000000005) -- (193.49,327.62027500000005);
                \draw[branch] (229.18,383.55325000000005) |- (238.9,361.37027500000005) -| (248.62,383.55325000000005);
                \draw[branch] (224.945,343.87027500000005) -- (224.945,327.62027500000005);
                \draw[branch] (210.99,383.55325000000005) |- (224.945,343.87027500000005) -| (238.9,361.37027500000005);
                \draw[branch] (390.520275,383.55325000000005) -- (390.520275,362.62027500000005);
                \draw[branch] (411.910275,383.55325000000005) -- (411.910275,362.62027500000005);
                \draw[branch] (431.910275,383.55325000000005) -- (431.910275,362.62027500000005);
                \draw[branch] (627.9202750000001,293.31027500000005) -- (627.9202750000001,267.06027500000005);
                \draw[branch] (583.6102750000001,345.12027500000005) |- (627.9202750000001,293.31027500000005) -| (680.0102750000001,327.62027500000005);
                \draw[branch] (552.9802750000001,361.37027500000005) -- (552.9802750000001,345.12027500000005);
                \draw[branch] (544.0202750000001,383.55325000000005) |- (552.9802750000001,361.37027500000005) -| (561.9402750000002,383.55325000000005);
                \draw[branch] (583.6102750000001,383.55325000000005) -- (583.6102750000001,345.12027500000005);
                \draw[branch] (669.7302750000001,383.55325000000005) |- (680.0102750000001,361.37027500000005) -| (690.2902750000002,383.55325000000005);
                \draw[branch] (680.0102750000001,343.87027500000005) -- (680.0102750000001,327.62027500000005);
                \draw[branch] (680.0102750000001,361.37027500000005) |- (680.0102750000001,343.87027500000005);
                % gene transfers
                \draw[transfer branch] (292.51027500000004,116.55000000000001) to[bend left=35] (552.9802750000001,345.12027500000005);
                \draw[transfer branch] (680.0102750000001,343.87027500000005) to[bend right=35] (431.910275,362.62027500000005);
                % events
                \node[speciation] at (483.86027500000006,38.75) {abcdefg};
                \node[duplication] at (497.82027500000004,16.25) {abcdefg};
                \node[speciation] at (292.51027500000004,139.05) {abcd};
                \node[horizontal gene transfer] at (292.51027500000004,116.55000000000001) {abcd};
                \node[speciation] at (329.18027500000005,156.55) {cdef};
                \node[duplication] at (339.18027500000005,116.55000000000001) {abcdefg};
                \node[speciation] at (111.68,253.31027500000005) {defg};
                \node[speciation] at (146.96,270.81027500000005) {cdef};
                \node[extant gene={abcd}] at (20.42,383.55325000000005) {};
                \node[extant gene={defg}] at (44.870000000000005,383.55325000000005) {};
                \node[extant gene={cdef}] at (67.65,383.55325000000005) {};
                \node[extant gene={def}] at (193.49,383.55325000000005) {};
                \node[extant gene={cef}] at (210.99,383.55325000000005) {};
                \node[extant gene={cde}] at (229.18,383.55325000000005) {};
                \node[extant gene={cde}] at (248.62,383.55325000000005) {};
                \node[duplication] at (238.9,361.37027500000005) {cde};
                \node[duplication] at (224.945,343.87027500000005) {cdef};
                \node[extant gene={abcd}] at (390.520275,383.55325000000005) {};
                \node[extant gene={cef}] at (411.910275,383.55325000000005) {};
                \node[extant gene={defg}] at (431.910275,383.55325000000005) {};
                \node[speciation] at (627.9202750000001,293.31027500000005) {defg};
                \node[extant gene={ab}] at (544.0202750000001,383.55325000000005) {};
                \node[extant gene={abc}] at (561.9402750000002,383.55325000000005) {};
                \node[duplication] at (552.9802750000001,361.37027500000005) {abc};
                \node[extant gene={defg}] at (583.6102750000001,383.55325000000005) {};
                \node[extant gene={def}] at (669.7302750000001,383.55325000000005) {};
                \node[extant gene={defg}] at (690.2902750000002,383.55325000000005) {};
                \node[duplication] at (680.0102750000001,361.37027500000005) {defg};
                \node[horizontal gene transfer] at (680.0102750000001,343.87027500000005) {defg};
                """
            ).lstrip(),
        )
