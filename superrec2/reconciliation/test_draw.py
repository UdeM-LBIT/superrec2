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
        layout_info = compute_layout(species_tree, rec, labeling)
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
                \draw[species border] (49.763,72.5) |- (152.026,20) -- (152.026,0);
                \draw[species border] (253.449,145.0) |- (204.526,20) -- (204.526,0);
                \draw[species border] (102.263,72.5) |- (152.026,72.5) -| (204.526,145.0);
                \draw[species border] (0,145.0) |- (49.763,92.5) -- (49.763,72.5);
                \draw[species border] (152.026,145.0) |- (102.263,92.5) -- (102.263,72.5);
                \draw[species border] (49.763,145.0) |- (49.763,145.0) -| (102.263,145.0);
                \draw[species border] (0,145.0) -- (0,196.86595) -- node[species label] {X} (49.763,196.86595) -- (49.763,145.0);
                \draw[species border] (102.263,145.0) -- (102.263,196.86595) -- node[species label] {Y} (152.026,196.86595) -- (152.026,145.0);
                \draw[species border] (204.526,145.0) -- (204.526,196.86595) -- node[species label] {Z} (253.449,196.86595) -- (253.449,145.0);
                % gene branches
                \draw[branch] (178.276,46.25) -- (178.276,0);
                \draw[branch] (76.013,72.5) |- (178.276,46.25) -| (228.9875,145.0);
                \draw[branch] (76.013,118.75) -- (76.013,72.5);
                \draw[branch] (24.8815,145.0) |- (76.013,118.75) -| (127.14450000000001,145.0);
                \draw[branch] (24.8815,175.932975) -- (24.8815,145.0);
                \draw[branch] (127.14450000000001,175.932975) -- (127.14450000000001,145.0);
                \draw[branch] (228.9875,175.932975) -- (228.9875,145.0);
                % gene transfers
                % events
                \node[speciation] at (178.276,46.25) {};
                \node[speciation] at (76.013,118.75) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,175.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (127.14450000000001,175.932975) {};
                \node[extant gene={z\textsubscript{1}}] at (228.9875,175.932975) {};
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
                \draw[species border] (89.289,152.5) |- (276.078,55.0) -- (276.078,0);
                \draw[species border] (460.347,270.0) |- (373.578,55.0) -- (373.578,0);
                \draw[species border] (186.789,152.5) |- (276.078,152.5) -| (373.578,270.0);
                \draw[species border] (0,270.0) |- (89.289,172.5) -- (89.289,152.5);
                \draw[species border] (276.078,270.0) |- (186.789,172.5) -- (186.789,152.5);
                \draw[species border] (89.289,270.0) |- (89.289,270.0) -| (186.789,270.0);
                \draw[species border] (0,270.0) -- (0,321.86595) -- node[species label] {X} (89.289,321.86595) -- (89.289,270.0);
                \draw[species border] (186.789,270.0) -- (186.789,321.86595) -- node[species label] {Y} (276.078,321.86595) -- (276.078,270.0);
                \draw[species border] (373.578,270.0) -- (373.578,321.86595) -- node[species label] {Z} (460.347,321.86595) -- (460.347,270.0);
                % gene branches
                \draw[branch] (115.539,152.5) |- (302.328,81.25) -| (416.9625,270.0);
                \draw[branch] (138.039,152.5) |- (324.828,103.75) -| (435.8855,270.0);
                \draw[branch] (160.539,152.5) |- (347.328,126.25) -| (398.0395,270.0);
                \draw[branch] (302.328,81.25) |- (313.578,48.75) -| (324.828,103.75);
                \draw[branch] (330.453,26.25) -- (330.453,0);
                \draw[branch] (347.328,126.25) |- (330.453,26.25) -| (313.578,48.75);
                \draw[branch] (115.539,198.75) -- (115.539,152.5);
                \draw[branch] (44.6445,270.0) |- (115.539,198.75) -| (231.43349999999998,270.0);
                \draw[branch] (138.039,221.25) -- (138.039,152.5);
                \draw[branch] (64.4075,270.0) |- (138.039,221.25) -| (251.1965,270.0);
                \draw[branch] (160.539,243.75) -- (160.539,152.5);
                \draw[branch] (24.8815,270.0) |- (160.539,243.75) -| (211.67049999999998,270.0);
                \draw[branch] (24.8815,300.932975) -- (24.8815,270.0);
                \draw[branch] (44.6445,300.932975) -- (44.6445,270.0);
                \draw[branch] (64.4075,300.932975) -- (64.4075,270.0);
                \draw[branch] (211.67049999999998,300.932975) -- (211.67049999999998,270.0);
                \draw[branch] (231.43349999999998,300.932975) -- (231.43349999999998,270.0);
                \draw[branch] (251.1965,300.932975) -- (251.1965,270.0);
                \draw[branch] (398.0395,300.932975) -- (398.0395,270.0);
                \draw[branch] (416.9625,300.932975) -- (416.9625,270.0);
                \draw[branch] (435.8855,300.932975) -- (435.8855,270.0);
                % gene transfers
                % events
                \node[speciation] at (302.328,81.25) {};
                \node[speciation] at (324.828,103.75) {};
                \node[speciation] at (347.328,126.25) {};
                \node[duplication] at (313.578,48.75) {};
                \node[duplication] at (330.453,26.25) {};
                \node[speciation] at (115.539,198.75) {};
                \node[speciation] at (138.039,221.25) {};
                \node[speciation] at (160.539,243.75) {};
                \node[extant gene={x\textsubscript{3}}] at (24.8815,300.932975) {};
                \node[extant gene={x\textsubscript{1}}] at (44.6445,300.932975) {};
                \node[extant gene={x\textsubscript{2}}] at (64.4075,300.932975) {};
                \node[extant gene={y\textsubscript{3}}] at (211.67049999999998,300.932975) {};
                \node[extant gene={y\textsubscript{2}}] at (231.43349999999998,300.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (251.1965,300.932975) {};
                \node[extant gene={z\textsubscript{3}}] at (398.0395,300.932975) {};
                \node[extant gene={z\textsubscript{2}}] at (416.9625,300.932975) {};
                \node[extant gene={z\textsubscript{1}}] at (435.8855,300.932975) {};
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
                \draw[species border] (0,192.5) |- (49.763,20) -- (49.763,0);
                \draw[species border] (182.263,72.5) |- (102.263,20) -- (102.263,0);
                \draw[species border] (49.763,192.5) |- (49.763,72.5) -| (142.263,72.5);
                \draw[species border] (0,192.5) -- (0,244.36595) -- node[species label] {X} (49.763,244.36595) -- (49.763,192.5);
                \draw[species border] (102.263,214.36595) |- (142.263,92.5) -- (142.263,72.5);
                \draw[species border] (271.18600000000004,132.5) |- (182.263,92.5) -- (182.263,72.5);
                \draw[species border] (142.263,214.36595) |- (142.263,132.5) -| (231.186,132.5);
                \draw[species border] (102.263,214.36595) -- (102.263,244.36595) -- node[species label] {Y} (142.263,244.36595) -- (142.263,214.36595);
                \draw[species border] (182.263,192.5) |- (231.186,152.5) -- (231.186,132.5);
                \draw[species border] (311.18600000000004,214.36595) |- (271.18600000000004,152.5) -- (271.18600000000004,132.5);
                \draw[species border] (231.186,192.5) |- (231.186,192.5) -| (271.18600000000004,214.36595);
                \draw[species border] (182.263,192.5) -- (182.263,244.36595) -- node[species label] {Z} (231.186,244.36595) -- (231.186,192.5);
                \draw[species border] (271.18600000000004,214.36595) -- (271.18600000000004,244.36595) -- node[species label] {W} (311.18600000000004,244.36595) -- (311.18600000000004,214.36595);
                % gene branches
                \draw[branch] (76.013,46.25) -- (76.013,0);
                \draw[branch] (24.8815,192.5) |- (76.013,46.25) -| (162.263,72.5);
                \draw[branch] (24.8815,223.432975) -- (24.8815,192.5);
                \draw[branch] (162.263,112.5) -- (162.263,72.5);
                \draw[loss] (162.263,112.5) -- ++(-20, 0);
                \draw[branch] (162.263,112.5) -| (251.186,132.5);
                \draw[branch] (251.186,172.5) -- (251.186,132.5);
                \draw[loss] (251.186,172.5) -- ++(20, 0);
                \draw[branch] (251.186,172.5) -| (206.7245,192.5);
                \draw[branch] (206.7245,223.432975) -- (206.7245,192.5);
                % gene transfers
                % events
                \node[speciation] at (76.013,46.25) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,223.432975) {};
                \node[extant gene={z\textsubscript{1}}] at (206.7245,223.432975) {};
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
                \draw[species border] (49.763,72.5) |- (129.763,20) -- (129.763,0);
                \draw[species border] (262.26300000000003,72.5) |- (182.263,20) -- (182.263,0);
                \draw[species border] (89.763,72.5) |- (129.763,72.5) -| (222.263,72.5);
                \draw[species border] (0,132.5) |- (49.763,92.5) -- (49.763,72.5);
                \draw[species border] (129.763,154.36595) |- (89.763,92.5) -- (89.763,72.5);
                \draw[species border] (49.763,132.5) |- (49.763,132.5) -| (89.763,154.36595);
                \draw[species border] (0,132.5) -- (0,184.36595) -- node[species label] {X} (49.763,184.36595) -- (49.763,132.5);
                \draw[species border] (89.763,154.36595) -- (89.763,184.36595) -- node[species label] {Y} (129.763,184.36595) -- (129.763,154.36595);
                \draw[species border] (182.263,154.36595) |- (222.263,92.5) -- (222.263,72.5);
                \draw[species border] (313.966,132.5) |- (262.26300000000003,92.5) -- (262.26300000000003,72.5);
                \draw[species border] (222.263,154.36595) |- (222.263,132.5) -| (262.26300000000003,132.5);
                \draw[species border] (182.263,154.36595) -- (182.263,184.36595) -- node[species label] {Z} (222.263,184.36595) -- (222.263,154.36595);
                \draw[species border] (262.26300000000003,132.5) -- (262.26300000000003,184.36595) -- node[species label] {W} (313.966,184.36595) -- (313.966,132.5);
                % gene branches
                \draw[branch] (156.013,46.25) -- (156.013,0);
                \draw[branch] (69.763,72.5) |- (156.013,46.25) -| (242.263,72.5);
                \draw[branch] (69.763,112.5) -- (69.763,72.5);
                \draw[loss] (69.763,112.5) -- ++(20, 0);
                \draw[branch] (69.763,112.5) -| (24.8815,132.5);
                \draw[branch] (24.8815,163.432975) -- (24.8815,132.5);
                \draw[branch] (242.263,112.5) -- (242.263,72.5);
                \draw[loss] (242.263,112.5) -- ++(-20, 0);
                \draw[branch] (242.263,112.5) -| (288.1145,132.5);
                \draw[branch] (288.1145,163.432975) -- (288.1145,132.5);
                % gene transfers
                % events
                \node[speciation] at (156.013,46.25) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,163.432975) {};
                \node[extant gene={w\textsubscript{1}}] at (288.1145,163.432975) {};
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
                \draw[species border] (69.526,140.0) |- (214.05200000000002,55.0) -- (214.05200000000002,0);
                \draw[species border] (385.821,235.0) |- (299.052,55.0) -- (299.052,0);
                \draw[species border] (144.526,140.0) |- (214.05200000000002,140.0) -| (299.052,235.0);
                \draw[species border] (0,235.0) |- (69.526,160.0) -- (69.526,140.0);
                \draw[species border] (214.05200000000002,235.0) |- (144.526,160.0) -- (144.526,140.0);
                \draw[species border] (69.526,235.0) |- (69.526,235.0) -| (144.526,235.0);
                \draw[species border] (0,235.0) -- (0,286.86595) -- node[species label] {X} (69.526,286.86595) -- (69.526,235.0);
                \draw[species border] (144.526,235.0) -- (144.526,286.86595) -- node[species label] {Y} (214.05200000000002,286.86595) -- (214.05200000000002,235.0);
                \draw[species border] (299.052,235.0) -- (299.052,286.86595) -- node[species label] {Z} (385.821,286.86595) -- (385.821,235.0);
                % gene branches
                \draw[branch] (95.776,140.0) |- (240.30200000000002,81.25) -| (342.4365,235.0);
                \draw[branch] (118.276,140.0) |- (262.802,103.75) -| (361.3595,235.0);
                \draw[branch] (240.30200000000002,81.25) |- (251.55200000000002,48.75) -| (262.802,103.75);
                \draw[loss] (279.052,120.0) -- ++(-20, 0);
                \draw[branch] (279.052,120.0) -| (323.5135,235.0);
                \draw[branch] (265.302,26.25) -- (265.302,0);
                \draw[branch] (279.052,120.0) |- (265.302,26.25) -| (251.55200000000002,48.75);
                \draw[branch] (95.776,186.25) -- (95.776,140.0);
                \draw[branch] (24.8815,235.0) |- (95.776,186.25) -| (169.4075,235.0);
                \draw[branch] (118.276,208.75) -- (118.276,140.0);
                \draw[branch] (44.6445,235.0) |- (118.276,208.75) -| (189.1705,235.0);
                \draw[branch] (24.8815,265.932975) -- (24.8815,235.0);
                \draw[branch] (44.6445,265.932975) -- (44.6445,235.0);
                \draw[branch] (169.4075,265.932975) -- (169.4075,235.0);
                \draw[branch] (189.1705,265.932975) -- (189.1705,235.0);
                \draw[branch] (323.5135,265.932975) -- (323.5135,235.0);
                \draw[branch] (342.4365,265.932975) -- (342.4365,235.0);
                \draw[branch] (361.3595,265.932975) -- (361.3595,235.0);
                % gene transfers
                % events
                \node[speciation] at (240.30200000000002,81.25) {};
                \node[speciation] at (262.802,103.75) {};
                \node[duplication] at (251.55200000000002,48.75) {};
                \node[duplication] at (265.302,26.25) {};
                \node[speciation] at (95.776,186.25) {};
                \node[speciation] at (118.276,208.75) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,265.932975) {};
                \node[extant gene={x\textsubscript{2}}] at (44.6445,265.932975) {};
                \node[extant gene={y\textsubscript{2}}] at (169.4075,265.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (189.1705,265.932975) {};
                \node[extant gene={z\textsubscript{3}}] at (323.5135,265.932975) {};
                \node[extant gene={z\textsubscript{2}}] at (342.4365,265.932975) {};
                \node[extant gene={z\textsubscript{1}}] at (361.3595,265.932975) {};
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
                \draw[species border] (0,175.0) |- (109.05199999999999,55.0) -- (109.05199999999999,0);
                \draw[species border] (338.104,175.0) |- (229.052,55.0) -- (229.052,0);
                \draw[species border] (109.05199999999999,175.0) |- (109.05199999999999,175.0) -| (229.052,175.0);
                \draw[species border] (0,175.0) -- (0,226.86595) -- node[species label] {X} (109.05199999999999,226.86595) -- (109.05199999999999,175.0);
                \draw[species border] (229.052,175.0) -- (229.052,226.86595) -- node[species label] {Y} (338.104,226.86595) -- (338.104,175.0);
                % gene branches
                \draw[branch] (24.8815,175.0) |- (135.302,81.25) -| (253.93349999999998,175.0);
                \draw[branch] (44.6445,175.0) |- (157.802,103.75) -| (273.6965,175.0);
                \draw[branch] (64.4075,175.0) |- (180.302,126.25) -| (293.4595,175.0);
                \draw[branch] (84.17049999999999,175.0) |- (202.802,148.75) -| (313.22249999999997,175.0);
                \draw[branch] (135.302,81.25) |- (146.552,48.75) -| (157.802,103.75);
                \draw[branch] (180.302,126.25) |- (191.552,48.75) -| (202.802,148.75);
                \draw[branch] (169.052,26.25) -- (169.052,0);
                \draw[branch] (146.552,48.75) |- (169.052,26.25) -| (191.552,48.75);
                \draw[branch] (24.8815,205.932975) -- (24.8815,175.0);
                \draw[branch] (44.6445,205.932975) -- (44.6445,175.0);
                \draw[branch] (64.4075,205.932975) -- (64.4075,175.0);
                \draw[branch] (84.17049999999999,205.932975) -- (84.17049999999999,175.0);
                \draw[branch] (253.93349999999998,205.932975) -- (253.93349999999998,175.0);
                \draw[branch] (273.6965,205.932975) -- (273.6965,175.0);
                \draw[branch] (293.4595,205.932975) -- (293.4595,175.0);
                \draw[branch] (313.22249999999997,205.932975) -- (313.22249999999997,175.0);
                % gene transfers
                % events
                \node[speciation] at (135.302,81.25) {};
                \node[speciation] at (157.802,103.75) {};
                \node[speciation] at (180.302,126.25) {};
                \node[speciation] at (202.802,148.75) {};
                \node[duplication] at (146.552,48.75) {};
                \node[duplication] at (191.552,48.75) {};
                \node[duplication] at (169.052,26.25) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,205.932975) {};
                \node[extant gene={x\textsubscript{2}}] at (44.6445,205.932975) {};
                \node[extant gene={x\textsubscript{3}}] at (64.4075,205.932975) {};
                \node[extant gene={x\textsubscript{4}}] at (84.17049999999999,205.932975) {};
                \node[extant gene={y\textsubscript{1}}] at (253.93349999999998,205.932975) {};
                \node[extant gene={y\textsubscript{2}}] at (273.6965,205.932975) {};
                \node[extant gene={y\textsubscript{3}}] at (293.4595,205.932975) {};
                \node[extant gene={y\textsubscript{4}}] at (313.22249999999997,205.932975) {};
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
                \draw[species border] (49.763,82.75) |- (149.526,32.5) -- (149.526,0);
                \draw[species border] (239.776,174.61595) |- (199.776,32.5) -- (199.776,0);
                \draw[species border] (99.763,82.75) |- (149.526,82.5) -| (199.776,174.61595);
                \draw[species border] (0,152.75) |- (49.763,102.75) -- (49.763,82.75);
                \draw[species border] (149.526,152.75) |- (99.763,102.75) -- (99.763,82.75);
                \draw[species border] (49.763,152.75) |- (49.763,152.75) -| (99.763,152.75);
                \draw[species border] (0,152.75) -- (0,204.61595) -- node[species label] {X} (49.763,204.61595) -- (49.763,152.75);
                \draw[species border] (99.763,152.75) -- (99.763,204.61595) -- node[species label] {Y} (149.526,204.61595) -- (149.526,152.75);
                \draw[species border] (199.776,174.61595) -- (199.776,204.61595) -- node[species label] {Z} (239.776,204.61595) -- (239.776,174.61595);
                % gene branches
                \draw[loss] (169.526,52.5) -- ++(20, 0);
                \draw[branch] (169.526,52.5) -| (69.763,82.75);
                \draw[loss] (179.526,62.5) -- ++(20, 0);
                \draw[branch] (179.526,62.5) -| (79.763,82.75);
                \draw[branch] (174.526,26.25) -- (174.526,0);
                \draw[branch] (169.526,52.5) |- (174.526,26.25) -| (179.526,62.5);
                \draw[branch] (69.763,122.75) -- (69.763,82.75);
                \draw[loss] (69.763,122.75) -- ++(-20, 0);
                \draw[branch] (69.763,122.75) -| (124.64450000000001,152.75);
                \draw[branch] (79.763,132.75) -- (79.763,82.75);
                \draw[loss] (79.763,132.75) -- ++(20, 0);
                \draw[branch] (79.763,132.75) -| (24.8815,152.75);
                \draw[branch] (24.8815,183.682975) -- (24.8815,152.75);
                \draw[branch] (124.64450000000001,183.682975) -- (124.64450000000001,152.75);
                % gene transfers
                % events
                \node[duplication] at (174.526,26.25) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,183.682975) {};
                \node[extant gene={y\textsubscript{1}}] at (124.64450000000001,183.682975) {};
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
                \draw[species border] (69.526,95.0) |- (214.052275,32.5) -- (214.052275,0);
                \draw[species border] (344.398275,202.500275) |- (276.552275,32.5) -- (276.552275,0);
                \draw[species border] (144.526275,95.0) |- (214.052275,95.0) -| (276.552275,202.500275);
                \draw[species border] (0,202.500275) |- (69.526,127.5) -- (69.526,95.0);
                \draw[species border] (214.052275,202.500275) |- (144.526275,127.5) -- (144.526275,95.0);
                \draw[species border] (69.526,202.500275) |- (69.526,202.5) -| (144.526275,202.500275);
                \draw[species border] (0,202.500275) -- (0,254.366225) -- node[species label] {X} (69.526,254.366225) -- (69.526,202.500275);
                \draw[species border] (144.526275,202.500275) -- (144.526275,254.366225) -- node[species label] {Y} (214.052275,254.366225) -- (214.052275,202.500275);
                \draw[species border] (276.552275,202.500275) -- (276.552275,254.366225) -- node[species label] {Z} (344.398275,254.366225) -- (344.398275,202.500275);
                % gene branches
                \draw[branch] (118.276,95.0) |- (240.302275,58.75) -| (319.936775,202.500275);
                \draw[loss] (256.552275,75.0) -- ++(20, 0);
                \draw[branch] (256.552275,75.0) -| (95.776,95.0);
                \draw[branch] (248.427275,26.25) -- (248.427275,0);
                \draw[branch] (256.552275,75.0) |- (248.427275,26.25) -| (240.302275,58.75);
                \draw[branch] (95.776,153.75) -- (95.776,95.0);
                \draw[branch] (24.8815,202.500275) |- (95.776,153.75) -| (169.407775,202.500275);
                \draw[branch] (44.6445,202.500275) |- (118.276,176.25) -| (189.170775,202.500275);
                \draw[branch] (118.276,121.25) -- (118.276,95.0);
                \draw[branch] (118.276,176.25) |- (118.276,121.25);
                \draw[branch] (24.8815,233.43325) -- (24.8815,202.500275);
                \draw[branch] (44.6445,233.43325) -- (44.6445,202.500275);
                \draw[branch] (169.407775,233.43325) -- (169.407775,202.500275);
                \draw[branch] (189.170775,233.43325) -- (189.170775,202.500275);
                \draw[branch] (301.013775,233.43325) -- (301.013775,202.500275);
                \draw[branch] (319.936775,233.43325) -- (319.936775,202.500275);
                % gene transfers
                \draw[transfer branch] (118.276,121.25) to[bend left=35] (301.013775,202.500275);
                % events
                \node[speciation] at (240.302275,58.75) {};
                \node[duplication] at (248.427275,26.25) {};
                \node[speciation] at (95.776,153.75) {};
                \node[speciation] at (118.276,176.25) {};
                \node[horizontal gene transfer] at (118.276,121.25) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,233.43325) {};
                \node[extant gene={x\textsubscript{2}}] at (44.6445,233.43325) {};
                \node[extant gene={y\textsubscript{1}}] at (169.407775,233.43325) {};
                \node[extant gene={y\textsubscript{2}}] at (189.170775,233.43325) {};
                \node[extant gene={z\textsubscript{1}}] at (301.013775,233.43325) {};
                \node[extant gene={z\textsubscript{2}}] at (319.936775,233.43325) {};
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
                \draw[species border] (0,145.000275) |- (49.763,20) -- (49.763,0);
                \draw[species border] (192.026275,59.999999999999986) |- (89.763,20) -- (89.763,0);
                \draw[species border] (49.763,145.000275) |- (49.763,40) -| (139.526,59.999999999999986);
                \draw[species border] (0,145.000275) -- (0,196.866225) -- node[species label] {X} (49.763,196.866225) -- (49.763,145.000275);
                \draw[species border] (89.763,145.000275) |- (139.526,92.49999999999999) -- (139.526,59.999999999999986);
                \draw[species border] (240.949275,145.000275) |- (192.026275,92.49999999999999) -- (192.026275,59.999999999999986);
                \draw[species border] (139.526,145.000275) |- (139.526,145.0) -| (192.026275,145.000275);
                \draw[species border] (89.763,145.000275) -- (89.763,196.866225) -- node[species label] {Y} (139.526,196.866225) -- (139.526,145.000275);
                \draw[species border] (192.026275,145.000275) -- (192.026275,196.866225) -- node[species label] {Z} (240.949275,196.866225) -- (240.949275,145.000275);
                % gene branches
                \draw[branch] (24.8815,175.93325) -- (24.8815,145.000275);
                \draw[branch] (114.64450000000001,145.000275) |- (165.776,118.74999999999999) -| (216.487775,145.000275);
                \draw[branch] (165.776,86.24999999999999) -- (165.776,59.999999999999986);
                \draw[branch] (165.776,118.74999999999999) |- (165.776,86.24999999999999);
                \draw[branch] (114.64450000000001,175.93325) -- (114.64450000000001,145.000275);
                \draw[branch] (216.487775,175.93325) -- (216.487775,145.000275);
                % gene transfers
                \draw[transfer branch] (165.776,86.24999999999999) to[bend right=35] (24.8815,145.000275);
                % events
                \node[extant gene={x\textsubscript{1}}] at (24.8815,175.93325) {};
                \node[speciation] at (165.776,118.74999999999999) {};
                \node[horizontal gene transfer] at (165.776,86.24999999999999) {};
                \node[extant gene={y\textsubscript{1}}] at (114.64450000000001,175.93325) {};
                \node[extant gene={z\textsubscript{1}}] at (216.487775,175.93325) {};
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
                \draw[species border] (283.341,95.0) |- (455.11,32.5) -- (455.11,0);
                \draw[species border] (665.219,245.0) |- (517.61,32.5) -- (517.61,0);
                \draw[species border] (368.341,95.0) |- (455.11,95.0) -| (612.719,245.0);
                \draw[species border] (89.289,212.5) |- (283.341,127.5) -- (283.341,95.0);
                \draw[species border] (455.11,362.5) |- (368.341,127.5) -- (368.341,95.0);
                \draw[species border] (174.289,212.5) |- (283.341,212.5) -| (368.341,362.5);
                \draw[species border] (0,362.5) |- (89.289,232.5) -- (89.289,212.5);
                \draw[species border] (283.341,317.5) |- (174.289,232.5) -- (174.289,212.5);
                \draw[species border] (89.289,362.5) |- (89.289,317.5) -| (174.289,317.5);
                \draw[species border] (0,362.5) -- (0,414.36595) -- node[species label] {X} (89.289,414.36595) -- (89.289,362.5);
                \draw[species border] (174.289,317.5) -- (174.289,414.36595) -- node[species label] {Y} (283.341,414.36595) -- (283.341,317.5);
                \draw[species border] (368.341,362.5) -- (368.341,414.36595) -- node[species label] {Z} (455.11,414.36595) -- (455.11,362.5);
                \draw[species border] (517.61,340.0) |- (612.719,265.0) -- (612.719,245.0);
                \draw[species border] (731.9649800000001,317.5) |- (665.219,265.0) -- (665.219,245.0);
                \draw[species border] (612.719,340.0) |- (612.719,317.5) -| (665.219,317.5);
                \draw[species border] (517.61,340.0) -- (517.61,414.36595) -- node[species label] {W} (612.719,414.36595) -- (612.719,340.0);
                \draw[species border] (665.219,317.5) -- (665.219,414.36595) -- node[species label] {T} (731.9649800000001,414.36595) -- (731.9649800000001,317.5);
                % gene branches
                \draw[branch] (340.216,95.0) |- (481.36,58.75) -| (638.969,245.0);
                \draw[loss] (497.61,75.0) -- ++(20, 0);
                \draw[branch] (497.61,75.0) -| (309.591,95.0);
                \draw[branch] (489.485,26.25) -- (489.485,0);
                \draw[branch] (497.61,75.0) |- (489.485,26.25) -| (481.36,58.75);
                \draw[branch] (154.289,212.5) |- (309.591,153.75) -| (392.8025,362.5);
                \draw[branch] (138.039,212.5) |- (332.091,176.25) -| (411.7255,362.5);
                \draw[branch] (309.591,121.25) -- (309.591,95.0);
                \draw[branch] (309.591,153.75) |- (309.591,121.25);
                \draw[loss] (348.341,192.5) -- ++(20, 0);
                \draw[branch] (348.341,192.5) -| (115.539,212.5);
                \draw[branch] (340.216,121.25) -- (340.216,95.0);
                \draw[branch] (348.341,192.5) |- (340.216,121.25) -| (332.091,176.25);
                \draw[branch] (115.539,258.75) -- (115.539,212.5);
                \draw[branch] (44.6445,362.5) |- (115.539,258.75) -| (199.17049999999998,317.5);
                \draw[branch] (138.039,281.25) -- (138.039,212.5);
                \draw[branch] (64.4075,362.5) |- (138.039,281.25) -| (233.75574999999998,317.5);
                \draw[branch] (154.289,297.5) -- (154.289,212.5);
                \draw[loss] (154.289,297.5) -- ++(20, 0);
                \draw[branch] (154.289,297.5) -| (24.8815,362.5);
                \draw[branch] (24.8815,393.432975) -- (24.8815,362.5);
                \draw[branch] (44.6445,393.432975) -- (44.6445,362.5);
                \draw[branch] (64.4075,393.432975) -- (64.4075,362.5);
                \draw[branch] (199.17049999999998,393.432975) -- (199.17049999999998,317.5);
                \draw[branch] (238.6965,393.432975) |- (248.57799999999997,366.25) -| (258.4595,393.432975);
                \draw[branch] (233.75574999999998,343.75) -- (233.75574999999998,317.5);
                \draw[branch] (218.93349999999998,393.432975) |- (233.75574999999998,343.75) -| (248.57799999999997,366.25);
                \draw[branch] (392.8025,393.432975) -- (392.8025,362.5);
                \draw[branch] (411.7255,393.432975) -- (411.7255,362.5);
                \draw[branch] (430.6485,393.432975) -- (430.6485,362.5);
                \draw[branch] (638.969,291.25) -- (638.969,245.0);
                \draw[branch] (586.8675000000001,340.0) |- (638.969,291.25) -| (698.59199,317.5);
                \draw[branch] (586.8675000000001,393.432975) -- (586.8675000000001,340.0);
                \draw[branch] (554.313,366.25) -- (554.313,340.0);
                \draw[branch] (543.4615,393.432975) |- (554.313,366.25) -| (565.1645,393.432975);
                \draw[branch] (689.4054950000001,393.432975) |- (698.59199,366.25) -| (707.778485,393.432975);
                \draw[branch] (698.59199,343.75) -- (698.59199,317.5);
                \draw[branch] (698.59199,366.25) |- (698.59199,343.75);
                % gene transfers
                \draw[transfer branch] (309.591,121.25) to[bend left=35] (554.313,340.0);
                \draw[transfer branch] (698.59199,343.75) to[bend right=35] (430.6485,362.5);
                % events
                \node[speciation] at (481.36,58.75) {};
                \node[duplication] at (489.485,26.25) {};
                \node[speciation] at (309.591,153.75) {};
                \node[speciation] at (332.091,176.25) {};
                \node[horizontal gene transfer] at (309.591,121.25) {};
                \node[duplication] at (340.216,121.25) {};
                \node[speciation] at (115.539,258.75) {};
                \node[speciation] at (138.039,281.25) {};
                \node[extant gene={x\textsubscript{1}}] at (24.8815,393.432975) {};
                \node[extant gene={x\textsubscript{2}}] at (44.6445,393.432975) {};
                \node[extant gene={x\textsubscript{3}}] at (64.4075,393.432975) {};
                \node[extant gene={y\textsubscript{4}}] at (199.17049999999998,393.432975) {};
                \node[extant gene={y\textsubscript{1}}] at (218.93349999999998,393.432975) {};
                \node[extant gene={y\textsubscript{2}}] at (238.6965,393.432975) {};
                \node[extant gene={y\textsubscript{3}}] at (258.4595,393.432975) {};
                \node[duplication] at (248.57799999999997,366.25) {};
                \node[duplication] at (233.75574999999998,343.75) {};
                \node[extant gene={z\textsubscript{1}}] at (392.8025,393.432975) {};
                \node[extant gene={z\textsubscript{2}}] at (411.7255,393.432975) {};
                \node[extant gene={z\textsubscript{3}}] at (430.6485,393.432975) {};
                \node[speciation] at (638.969,291.25) {};
                \node[extant gene={w\textsubscript{1}}] at (543.4615,393.432975) {};
                \node[extant gene={w\textsubscript{2}}] at (565.1645,393.432975) {};
                \node[extant gene={w\textsubscript{3}}] at (586.8675000000001,393.432975) {};
                \node[duplication] at (554.313,366.25) {};
                \node[extant gene={t\textsubscript{1}}] at (689.4054950000001,393.432975) {};
                \node[extant gene={t\textsubscript{2}}] at (707.778485,393.432975) {};
                \node[duplication] at (698.59199,366.25) {};
                \node[horizontal gene transfer] at (698.59199,343.75) {};
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
                \draw[species border] (360.84000000000003,132.8) |- (603.44,32.5) -- (603.44,0);
                \draw[species border] (878.2,347.06) |- (703.74,32.5) -- (703.74,0);
                \draw[species border] (492.6,132.8) |- (603.44,95.0) -| (807.64,347.06);
                \draw[species border] (116.4,297.06) |- (360.84000000000003,165.3) -- (360.84000000000003,132.8);
                \draw[species border] (603.44,482.62) |- (492.6,165.3) -- (492.6,132.8);
                \draw[species border] (236.96,297.06) |- (360.84000000000003,250.3) -| (492.6,482.62);
                \draw[species border] (0,482.62) |- (116.4,317.06) -- (116.4,297.06);
                \draw[species border] (360.84000000000003,437.62) |- (236.96,317.06) -- (236.96,297.06);
                \draw[species border] (116.4,482.62) |- (116.4,402.06) -| (236.96,437.62);
                \draw[species border] (0,482.62) -- (0,534.48595) -- node[species label] {X} (116.4,534.48595) -- (116.4,482.62);
                \draw[species border] (236.96,437.62) -- (236.96,534.48595) -- node[species label] {Y} (360.84000000000003,534.48595) -- (360.84000000000003,437.62);
                \draw[species border] (492.6,482.62) -- (492.6,534.48595) -- node[species label] {Z} (603.44,534.48595) -- (603.44,482.62);
                \draw[species border] (703.74,460.12) |- (807.64,367.06) -- (807.64,347.06);
                \draw[species border] (959.32,437.62) |- (878.2,367.06) -- (878.2,347.06);
                \draw[species border] (807.64,460.12) |- (807.64,419.56) -| (878.2,437.62);
                \draw[species border] (703.74,460.12) -- (703.74,534.48595) -- node[species label] {W} (807.64,534.48595) -- (807.64,460.12);
                \draw[species border] (878.2,437.62) -- (878.2,534.48595) -- node[species label] {T} (959.32,534.48595) -- (959.32,437.62);
                % gene branches
                \draw[branch] (451.68000000000006,132.8) |- (646.36,58.75) -| (842.92,347.06);
                \draw[loss] (679.2800000000001,75.0) -- ++(20, 0);
                \draw[branch] (679.2800000000001,75.0) -| (397.51000000000005,132.8);
                \draw[branch] (662.82,26.25) -- (662.82,0);
                \draw[branch] (679.2800000000001,75.0) |- (662.82,26.25) -| (646.36,58.75);
                \draw[branch] (216.96,297.06) |- (397.51000000000005,191.55) -| (523.02,482.62);
                \draw[branch] (191.96,297.06) |- (439.18000000000006,214.05) -| (549.4100000000001,482.62);
                \draw[branch] (397.51000000000005,159.05) -- (397.51000000000005,132.8);
                \draw[branch] (397.51000000000005,191.55) |- (397.51000000000005,159.05);
                \draw[loss] (464.18000000000006,230.3) -- ++(20, 0);
                \draw[branch] (464.18000000000006,230.3) -| (151.68,297.06);
                \draw[branch] (451.68000000000006,159.05) -- (451.68000000000006,132.8);
                \draw[branch] (464.18000000000006,230.3) |- (451.68000000000006,159.05) -| (439.18000000000006,214.05);
                \draw[branch] (151.68,343.31) -- (151.68,297.06);
                \draw[branch] (59.870000000000005,482.62) |- (151.68,343.31) -| (263.49,437.62);
                \draw[branch] (191.96,365.81) -- (191.96,297.06);
                \draw[branch] (87.65,482.62) |- (191.96,365.81) -| (303.695,437.62);
                \draw[branch] (216.96,382.06) -- (216.96,297.06);
                \draw[loss] (216.96,382.06) -- ++(20, 0);
                \draw[branch] (216.96,382.06) -| (30.42,482.62);
                \draw[branch] (30.42,513.552975) -- (30.42,482.62);
                \draw[branch] (59.870000000000005,513.552975) -- (59.870000000000005,482.62);
                \draw[branch] (87.65,513.552975) -- (87.65,482.62);
                \draw[branch] (263.49,513.552975) -- (263.49,437.62);
                \draw[branch] (309.18,513.552975) |- (321.4,486.37) -| (333.62,513.552975);
                \draw[branch] (303.695,463.87) -- (303.695,437.62);
                \draw[branch] (285.99,513.552975) |- (303.695,463.87) -| (321.4,486.37);
                \draw[branch] (523.02,513.552975) -- (523.02,482.62);
                \draw[branch] (549.4100000000001,513.552975) -- (549.4100000000001,482.62);
                \draw[branch] (574.4100000000001,513.552975) -- (574.4100000000001,482.62);
                \draw[branch] (842.92,393.31) -- (842.92,347.06);
                \draw[branch] (778.61,460.12) |- (842.92,393.31) -| (917.51,437.62);
                \draw[branch] (778.61,513.552975) -- (778.61,460.12);
                \draw[branch] (740.48,486.37) -- (740.48,460.12);
                \draw[branch] (729.02,513.552975) |- (740.48,486.37) -| (751.94,513.552975);
                \draw[branch] (904.73,513.552975) |- (917.51,486.37) -| (930.2900000000001,513.552975);
                \draw[branch] (917.51,463.87) -- (917.51,437.62);
                \draw[branch] (917.51,486.37) |- (917.51,463.87);
                % gene transfers
                \draw[transfer branch] (397.51000000000005,159.05) to[bend left=35] (740.48,460.12);
                \draw[transfer branch] (917.51,463.87) to[bend right=35] (574.4100000000001,482.62);
                % events
                \node[speciation] at (646.36,58.75) {abcdefg};
                \node[duplication] at (662.82,26.25) {abcdefg};
                \node[speciation] at (397.51000000000005,191.55) {abcd};
                \node[speciation] at (439.18000000000006,214.05) {cdef};
                \node[horizontal gene transfer] at (397.51000000000005,159.05) {abcd};
                \node[duplication] at (451.68000000000006,159.05) {abcdefg};
                \node[speciation] at (151.68,343.31) {defg};
                \node[speciation] at (191.96,365.81) {cdef};
                \node[extant gene={abcd}] at (30.42,513.552975) {};
                \node[extant gene={defg}] at (59.870000000000005,513.552975) {};
                \node[extant gene={cdef}] at (87.65,513.552975) {};
                \node[extant gene={def}] at (263.49,513.552975) {};
                \node[extant gene={cef}] at (285.99,513.552975) {};
                \node[extant gene={cde}] at (309.18,513.552975) {};
                \node[extant gene={cde}] at (333.62,513.552975) {};
                \node[duplication] at (321.4,486.37) {cde};
                \node[duplication] at (303.695,463.87) {cdef};
                \node[extant gene={abcd}] at (523.02,513.552975) {};
                \node[extant gene={cef}] at (549.4100000000001,513.552975) {};
                \node[extant gene={defg}] at (574.4100000000001,513.552975) {};
                \node[speciation] at (842.92,393.31) {defg};
                \node[extant gene={ab}] at (729.02,513.552975) {};
                \node[extant gene={abc}] at (751.94,513.552975) {};
                \node[extant gene={defg}] at (778.61,513.552975) {};
                \node[duplication] at (740.48,486.37) {abc};
                \node[extant gene={def}] at (904.73,513.552975) {};
                \node[extant gene={defg}] at (930.2900000000001,513.552975) {};
                \node[duplication] at (917.51,486.37) {defg};
                \node[horizontal gene transfer] at (917.51,463.87) {defg};
                """
            ).lstrip(),
        )
