import unittest
import textwrap
from ete3 import PhyloTree
from .tools import get_species_name, reconcile_leaves
from .draw import layout, render_to_tikz


class TestReconciliationDraw(unittest.TestCase):
    def test_speciations(self):
        gene_tree = PhyloTree(
            "((x_1,y_1)2,z_1)1;",
            sp_naming_function=get_species_name, format=1
        )

        species_tree = PhyloTree("((X,Y)XY,Z)XYZ;", format=1)

        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            gene_tree & "1": species_tree & "XYZ",
            gene_tree & "2": species_tree & "XY",
        }

        layout_info = layout(gene_tree, species_tree, rec)
        out = render_to_tikz(species_tree, rec, layout_info)

        self.assertEqual(out, textwrap.dedent(r"""
            \begin{tikzpicture}[
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
                    yshift=-10pt,
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
                    label={[label distance={0pt}]below:#1},
                },
                branch node/.style={
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={0.5pt},
                },
                speciation/.style={
                    branch node, circle,
                    minimum width={10pt},
                    minimum height={10pt},
                },
                duplication/.style={
                    branch node, rectangle,
                    minimum width={9pt},
                    minimum height={9pt},
                },
                horizontal gene transfer/.style={
                    branch node, diamond,
                    minimum width={10pt},
                    minimum height={10pt},
                },
            ]
            % species
            \draw[species border] (40,60) |- (120,20) -- (120,0);
            \draw[species border] (200,120) |- (160,20) -- (160,0);
            \draw[species border] (80,60) |- (120,60) -| (160,120);
            \draw[species border] (0,120) |- (40,80) -- (40,60);
            \draw[species border] (120,120) |- (80,80) -- (80,60);
            \draw[species border] (40,120) |- (40,120) -| (80,120);
            \draw[species border] (0,120) -- ([yshift=-16pt]0,140) -- node[species label] {X} ([yshift=-16pt]40,140) -- (40,120);
            \draw[species border] (80,120) -- ([yshift=-16pt]80,140) -- node[species label] {Y} ([yshift=-16pt]120,140) -- (120,120);
            \draw[species border] (160,120) -- ([yshift=-16pt]160,140) -- node[species label] {Z} ([yshift=-16pt]200,140) -- (200,120);
            % gene branches
            \draw[branch] (140,40) -- (140,0);
            \draw[branch] (60,60) |- (140,40) -| (180,120);
            \draw[branch] (60,100) -- (60,60);
            \draw[branch] (20,120) |- (60,100) -| (100,120);
            \draw[branch] (20,140) -- (20,120);
            \draw[branch] (100,140) -- (100,120);
            \draw[branch] (180,140) -- (180,120);
            % gene transfers
            % events
            \node[speciation] at (140,40) {};
            \node[speciation] at (60,100) {};
            \node[extant gene={\(x_1\)}] at (20,140) {};
            \node[extant gene={\(y_1\)}] at (100,140) {};
            \node[extant gene={\(z_1\)}] at (180,140) {};
            \end{tikzpicture}
        """).lstrip())

    def test_duplications(self):
        gene_tree = PhyloTree(
            "(((x_3,y_3)8,z_3)7,(((x_1,y_2)4,z_2)3,((x_2,y_1)6,z_1)5)2)1;",
            sp_naming_function=get_species_name, format=1
        )

        species_tree = PhyloTree("((X,Y)XY,Z)XYZ;", format=1)

        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            gene_tree & "1": species_tree & "XYZ",
            gene_tree & "2": species_tree & "XYZ",
            gene_tree & "3": species_tree & "XYZ",
            gene_tree & "4": species_tree & "XY",
            gene_tree & "5": species_tree & "XYZ",
            gene_tree & "6": species_tree & "XY",
            gene_tree & "7": species_tree & "XYZ",
            gene_tree & "8": species_tree & "XY",
        }

        layout_info = layout(gene_tree, species_tree, rec)
        out = render_to_tikz(species_tree, rec, layout_info)

        self.assertEqual(out, textwrap.dedent(r"""
            \begin{tikzpicture}[
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
                    yshift=-10pt,
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
                    label={[label distance={0pt}]below:#1},
                },
                branch node/.style={
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={0.5pt},
                },
                speciation/.style={
                    branch node, circle,
                    minimum width={10pt},
                    minimum height={10pt},
                },
                duplication/.style={
                    branch node, rectangle,
                    minimum width={9pt},
                    minimum height={9pt},
                },
                horizontal gene transfer/.style={
                    branch node, diamond,
                    minimum width={10pt},
                    minimum height={10pt},
                },
            ]
            % species
            \draw[species border] (70,105) |- (210,35) -- (210,0);
            \draw[species border] (350,195) |- (280,35) -- (280,0);
            \draw[species border] (140,105) |- (210,105) -| (280,195);
            \draw[species border] (0,195) |- (70,125) -- (70,105);
            \draw[species border] (210,195) |- (140,125) -- (140,105);
            \draw[species border] (70,195) |- (70,195) -| (140,195);
            \draw[species border] (0,195) -- ([yshift=-16pt]0,215) -- node[species label] {X} ([yshift=-16pt]70,215) -- (70,195);
            \draw[species border] (140,195) -- ([yshift=-16pt]140,215) -- node[species label] {Y} ([yshift=-16pt]210,215) -- (210,195);
            \draw[species border] (280,195) -- ([yshift=-16pt]280,215) -- node[species label] {Z} ([yshift=-16pt]350,215) -- (350,195);
            % gene branches
            \draw[branch] (120,105) |- (230,55) -| (300,195);
            \draw[branch] (90,105) |- (245,70) -| (315,195);
            \draw[branch] (105,105) |- (260,85) -| (330,195);
            \draw[branch] (245,70) |- (252.5,35) -| (260,85);
            \draw[branch] (241.25,20) -- (241.25,0);
            \draw[branch] (230,55) |- (241.25,20) -| (252.5,35);
            \draw[branch] (90,145) -- (90,105);
            \draw[branch] (35,195) |- (90,145) -| (175,195);
            \draw[branch] (105,160) -- (105,105);
            \draw[branch] (50,195) |- (105,160) -| (190,195);
            \draw[branch] (120,175) -- (120,105);
            \draw[branch] (20,195) |- (120,175) -| (160,195);
            \draw[branch] (20,215) -- (20,195);
            \draw[branch] (35,215) -- (35,195);
            \draw[branch] (50,215) -- (50,195);
            \draw[branch] (160,215) -- (160,195);
            \draw[branch] (175,215) -- (175,195);
            \draw[branch] (190,215) -- (190,195);
            \draw[branch] (300,215) -- (300,195);
            \draw[branch] (315,215) -- (315,195);
            \draw[branch] (330,215) -- (330,195);
            % gene transfers
            % events
            \node[speciation] at (230,55) {};
            \node[speciation] at (245,70) {};
            \node[speciation] at (260,85) {};
            \node[duplication] at (252.5,35) {};
            \node[duplication] at (241.25,20) {};
            \node[speciation] at (90,145) {};
            \node[speciation] at (105,160) {};
            \node[speciation] at (120,175) {};
            \node[extant gene={\(x_3\)}] at (20,215) {};
            \node[extant gene={\(x_1\)}] at (35,215) {};
            \node[extant gene={\(x_2\)}] at (50,215) {};
            \node[extant gene={\(y_3\)}] at (160,215) {};
            \node[extant gene={\(y_2\)}] at (175,215) {};
            \node[extant gene={\(y_1\)}] at (190,215) {};
            \node[extant gene={\(z_3\)}] at (300,215) {};
            \node[extant gene={\(z_2\)}] at (315,215) {};
            \node[extant gene={\(z_1\)}] at (330,215) {};
            \end{tikzpicture}
        """).lstrip())

    def test_speciations_losses(self):
        gene_tree = PhyloTree(
            "(x_1,z_1)1;",
            sp_naming_function=get_species_name, format=1
        )

        species_tree = PhyloTree("(X,(Y,(Z,W)ZW)YZW)XYZW;", format=1)

        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            gene_tree & "1": species_tree & "XYZW",
        }

        layout_info = layout(gene_tree, species_tree, rec)
        out = render_to_tikz(species_tree, rec, layout_info)

        self.assertEqual(out, textwrap.dedent(r"""
            \begin{tikzpicture}[
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
                    yshift=-10pt,
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
                    label={[label distance={0pt}]below:#1},
                },
                branch node/.style={
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={0.5pt},
                },
                speciation/.style={
                    branch node, circle,
                    minimum width={10pt},
                    minimum height={10pt},
                },
                duplication/.style={
                    branch node, rectangle,
                    minimum width={9pt},
                    minimum height={9pt},
                },
                horizontal gene transfer/.style={
                    branch node, diamond,
                    minimum width={10pt},
                    minimum height={10pt},
                },
            ]
            % species
            \draw[species border] (0,180) |- (40,20) -- (40,0);
            \draw[species border] (160,60) |- (80,20) -- (80,0);
            \draw[species border] (40,180) |- (40,60) -| (120,60);
            \draw[species border] (0,180) -- ([yshift=-16pt]0,200) -- node[species label] {X} ([yshift=-16pt]40,200) -- (40,180);
            \draw[species border] (80,180) |- (120,80) -- (120,60);
            \draw[species border] (240,120) |- (160,80) -- (160,60);
            \draw[species border] (120,180) |- (120,120) -| (200,120);
            \draw[species border] (80,180) -- ([yshift=-16pt]80,200) -- node[species label] {Y} ([yshift=-16pt]120,200) -- (120,180);
            \draw[species border] (160,180) |- (200,140) -- (200,120);
            \draw[species border] (280,180) |- (240,140) -- (240,120);
            \draw[species border] (200,180) |- (200,180) -| (240,180);
            \draw[species border] (160,180) -- ([yshift=-16pt]160,200) -- node[species label] {Z} ([yshift=-16pt]200,200) -- (200,180);
            \draw[species border] (240,180) -- ([yshift=-16pt]240,200) -- node[species label] {W} ([yshift=-16pt]280,200) -- (280,180);
            % gene branches
            \draw[branch] (60,40) -- (60,0);
            \draw[branch] (20,180) |- (60,40) -| (140,60);
            \draw[branch] (20,200) -- (20,180);
            \draw[branch] (140,100) -- (140,60);
            \draw[loss] (140,100) -- ++(-20pt, 0);
            \draw[branch] (140,100) -| (220,120);
            \draw[branch] (220,160) -- (220,120);
            \draw[loss] (220,160) -- ++(20pt, 0);
            \draw[branch] (220,160) -| (180,180);
            \draw[branch] (180,200) -- (180,180);
            % gene transfers
            % events
            \node[speciation] at (60,40) {};
            \node[extant gene={\(x_1\)}] at (20,200) {};
            \node[extant gene={\(z_1\)}] at (180,200) {};
            \end{tikzpicture}
        """).lstrip())

    def test_duplications_losses(self):
        gene_tree = PhyloTree(
            "(z_3,(((x_1,y_2)4,z_2)3,((x_2,y_1)6,z_1)5)2)1;",
            sp_naming_function=get_species_name, format=1
        )

        species_tree = PhyloTree("((X,Y)XY,Z)XYZ;", format=1)

        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            gene_tree & "1": species_tree & "XYZ",
            gene_tree & "2": species_tree & "XYZ",
            gene_tree & "3": species_tree & "XYZ",
            gene_tree & "4": species_tree & "XY",
            gene_tree & "5": species_tree & "XYZ",
            gene_tree & "6": species_tree & "XY",
        }

        layout_info = layout(gene_tree, species_tree, rec)
        out = render_to_tikz(species_tree, rec, layout_info)

        self.assertEqual(out, textwrap.dedent(r"""
            \begin{tikzpicture}[
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
                    yshift=-10pt,
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
                    label={[label distance={0pt}]below:#1},
                },
                branch node/.style={
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={0.5pt},
                },
                speciation/.style={
                    branch node, circle,
                    minimum width={10pt},
                    minimum height={10pt},
                },
                duplication/.style={
                    branch node, rectangle,
                    minimum width={9pt},
                    minimum height={9pt},
                },
                horizontal gene transfer/.style={
                    branch node, diamond,
                    minimum width={10pt},
                    minimum height={10pt},
                },
            ]
            % species
            \draw[species border] (55,105) |- (165,35) -- (165,0);
            \draw[species border] (305,180) |- (235,35) -- (235,0);
            \draw[species border] (110,105) |- (165,105) -| (235,180);
            \draw[species border] (0,180) |- (55,125) -- (55,105);
            \draw[species border] (165,180) |- (110,125) -- (110,105);
            \draw[species border] (55,180) |- (55,180) -| (110,180);
            \draw[species border] (0,180) -- ([yshift=-16pt]0,200) -- node[species label] {X} ([yshift=-16pt]55,200) -- (55,180);
            \draw[species border] (110,180) -- ([yshift=-16pt]110,200) -- node[species label] {Y} ([yshift=-16pt]165,200) -- (165,180);
            \draw[species border] (235,180) -- ([yshift=-16pt]235,200) -- node[species label] {Z} ([yshift=-16pt]305,200) -- (305,180);
            % gene branches
            \draw[branch] (75,105) |- (185,55) -| (270,180);
            \draw[branch] (90,105) |- (200,70) -| (285,180);
            \draw[branch] (185,55) |- (192.5,35) -| (200,70);
            \draw[loss] (215,85) -- ++(-20pt, 0);
            \draw[branch] (215,85) -| (255,180);
            \draw[branch] (203.75,20) -- (203.75,0);
            \draw[branch] (215,85) |- (203.75,20) -| (192.5,35);
            \draw[branch] (75,145) -- (75,105);
            \draw[branch] (20,180) |- (75,145) -| (130,180);
            \draw[branch] (90,160) -- (90,105);
            \draw[branch] (35,180) |- (90,160) -| (145,180);
            \draw[branch] (20,200) -- (20,180);
            \draw[branch] (35,200) -- (35,180);
            \draw[branch] (130,200) -- (130,180);
            \draw[branch] (145,200) -- (145,180);
            \draw[branch] (255,200) -- (255,180);
            \draw[branch] (270,200) -- (270,180);
            \draw[branch] (285,200) -- (285,180);
            % gene transfers
            % events
            \node[speciation] at (185,55) {};
            \node[speciation] at (200,70) {};
            \node[duplication] at (192.5,35) {};
            \node[duplication] at (203.75,20) {};
            \node[speciation] at (75,145) {};
            \node[speciation] at (90,160) {};
            \node[extant gene={\(x_1\)}] at (20,200) {};
            \node[extant gene={\(x_2\)}] at (35,200) {};
            \node[extant gene={\(y_2\)}] at (130,200) {};
            \node[extant gene={\(y_1\)}] at (145,200) {};
            \node[extant gene={\(z_3\)}] at (255,200) {};
            \node[extant gene={\(z_2\)}] at (270,200) {};
            \node[extant gene={\(z_1\)}] at (285,200) {};
            \end{tikzpicture}
        """).lstrip())

    def test_nested_duplications(self):
        gene_tree = PhyloTree(
            "(((x_1,y_1)4,(x_2,y_2)5)2,((x_3,y_3)6,(x_4,y_4)7)3)1;",
            sp_naming_function=get_species_name, format=1
        )

        species_tree = PhyloTree("(X,Y)XY;", format=1)

        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            gene_tree & "1": species_tree & "XY",
            gene_tree & "2": species_tree & "XY",
            gene_tree & "3": species_tree & "XY",
            gene_tree & "4": species_tree & "XY",
            gene_tree & "5": species_tree & "XY",
            gene_tree & "6": species_tree & "XY",
            gene_tree & "7": species_tree & "XY",
        }

        layout_info = layout(gene_tree, species_tree, rec)
        out = render_to_tikz(species_tree, rec, layout_info)

        self.assertEqual(out, textwrap.dedent(r"""
            \begin{tikzpicture}[
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
                    yshift=-10pt,
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
                    label={[label distance={0pt}]below:#1},
                },
                branch node/.style={
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={0.5pt},
                },
                speciation/.style={
                    branch node, circle,
                    minimum width={10pt},
                    minimum height={10pt},
                },
                duplication/.style={
                    branch node, rectangle,
                    minimum width={9pt},
                    minimum height={9pt},
                },
                horizontal gene transfer/.style={
                    branch node, diamond,
                    minimum width={10pt},
                    minimum height={10pt},
                },
            ]
            % species
            \draw[species border] (0,120) |- (85,35) -- (85,0);
            \draw[species border] (255,120) |- (170,35) -- (170,0);
            \draw[species border] (85,120) |- (85,120) -| (170,120);
            \draw[species border] (0,120) -- ([yshift=-16pt]0,140) -- node[species label] {X} ([yshift=-16pt]85,140) -- (85,120);
            \draw[species border] (170,120) -- ([yshift=-16pt]170,140) -- node[species label] {Y} ([yshift=-16pt]255,140) -- (255,120);
            % gene branches
            \draw[branch] (20,120) |- (105,55) -| (190,120);
            \draw[branch] (35,120) |- (120,70) -| (205,120);
            \draw[branch] (50,120) |- (135,85) -| (220,120);
            \draw[branch] (65,120) |- (150,100) -| (235,120);
            \draw[branch] (105,55) |- (112.5,35) -| (120,70);
            \draw[branch] (135,85) |- (142.5,35) -| (150,100);
            \draw[branch] (127.5,20) -- (127.5,0);
            \draw[branch] (112.5,35) |- (127.5,20) -| (142.5,35);
            \draw[branch] (20,140) -- (20,120);
            \draw[branch] (35,140) -- (35,120);
            \draw[branch] (50,140) -- (50,120);
            \draw[branch] (65,140) -- (65,120);
            \draw[branch] (190,140) -- (190,120);
            \draw[branch] (205,140) -- (205,120);
            \draw[branch] (220,140) -- (220,120);
            \draw[branch] (235,140) -- (235,120);
            % gene transfers
            % events
            \node[speciation] at (105,55) {};
            \node[speciation] at (120,70) {};
            \node[speciation] at (135,85) {};
            \node[speciation] at (150,100) {};
            \node[duplication] at (112.5,35) {};
            \node[duplication] at (142.5,35) {};
            \node[duplication] at (127.5,20) {};
            \node[extant gene={\(x_1\)}] at (20,140) {};
            \node[extant gene={\(x_2\)}] at (35,140) {};
            \node[extant gene={\(x_3\)}] at (50,140) {};
            \node[extant gene={\(x_4\)}] at (65,140) {};
            \node[extant gene={\(y_1\)}] at (190,140) {};
            \node[extant gene={\(y_2\)}] at (205,140) {};
            \node[extant gene={\(y_3\)}] at (220,140) {};
            \node[extant gene={\(y_4\)}] at (235,140) {};
            \end{tikzpicture}
        """).lstrip())

    def test_transfers(self):
        gene_tree = PhyloTree(
            "((x_1,y_1)2,(((x_2,y_2)5,z_1)4,z_2)3)1;",
            sp_naming_function=get_species_name, format=1
        )

        species_tree = PhyloTree("((X,Y)XY,Z)XYZ;", format=1)

        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            gene_tree & "1": species_tree,
            gene_tree & "2": species_tree & "XY",
            gene_tree & "3": species_tree,
            gene_tree & "4": species_tree & "XY",
            gene_tree & "5": species_tree & "XY",
        }

        layout_info = layout(gene_tree, species_tree, rec)
        out = render_to_tikz(species_tree, rec, layout_info)

        self.assertEqual(out, textwrap.dedent(r"""
            \begin{tikzpicture}[
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
                    yshift=-10pt,
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
                    label={[label distance={0pt}]below:#1},
                },
                branch node/.style={
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={0.5pt},
                },
                speciation/.style={
                    branch node, circle,
                    minimum width={10pt},
                    minimum height={10pt},
                },
                duplication/.style={
                    branch node, rectangle,
                    minimum width={9pt},
                    minimum height={9pt},
                },
                horizontal gene transfer/.style={
                    branch node, diamond,
                    minimum width={10pt},
                    minimum height={10pt},
                },
            ]
            % species
            \draw[species border] (55,75) |- (165,20) -- (165,0);
            \draw[species border] (275,150) |- (220,20) -- (220,0);
            \draw[species border] (110,75) |- (165,75) -| (220,150);
            \draw[species border] (0,150) |- (55,95) -- (55,75);
            \draw[species border] (165,150) |- (110,95) -- (110,75);
            \draw[species border] (55,150) |- (55,150) -| (110,150);
            \draw[species border] (0,150) -- ([yshift=-16pt]0,170) -- node[species label] {X} ([yshift=-16pt]55,170) -- (55,150);
            \draw[species border] (110,150) -- ([yshift=-16pt]110,170) -- node[species label] {Y} ([yshift=-16pt]165,170) -- (165,150);
            \draw[species border] (220,150) -- ([yshift=-16pt]220,170) -- node[species label] {Z} ([yshift=-16pt]275,170) -- (275,150);
            % gene branches
            \draw[branch] (90,75) |- (185,40) -| (255,150);
            \draw[loss] (200,55) -- ++(20pt, 0);
            \draw[branch] (200,55) -| (75,75);
            \draw[branch] (192.5,20) -- (192.5,0);
            \draw[branch] (200,55) |- (192.5,20) -| (185,40);
            \draw[branch] (75,115) -- (75,75);
            \draw[branch] (20,150) |- (75,115) -| (130,150);
            \draw[branch] (35,150) |- (90,130) -| (145,150);
            \draw[branch] (90,95) -- (90,75);
            \draw[branch] (90,130) |- (90,95);
            \draw[branch] (20,170) -- (20,150);
            \draw[branch] (35,170) -- (35,150);
            \draw[branch] (130,170) -- (130,150);
            \draw[branch] (145,170) -- (145,150);
            \draw[branch] (240,170) -- (240,150);
            \draw[branch] (255,170) -- (255,150);
            % gene transfers
            \draw[transfer branch] (90,95) to[bend left=35] (240,150);
            % events
            \node[speciation] at (185,40) {};
            \node[duplication] at (192.5,20) {};
            \node[speciation] at (75,115) {};
            \node[speciation] at (90,130) {};
            \node[horizontal gene transfer] at (90,95) {};
            \node[extant gene={\(x_1\)}] at (20,170) {};
            \node[extant gene={\(x_2\)}] at (35,170) {};
            \node[extant gene={\(y_1\)}] at (130,170) {};
            \node[extant gene={\(y_2\)}] at (145,170) {};
            \node[extant gene={\(z_1\)}] at (240,170) {};
            \node[extant gene={\(z_2\)}] at (255,170) {};
            \end{tikzpicture}
        """).lstrip())

    def test_all(self):
        gene_tree = PhyloTree("""
            (
                ((x_1,z_1)3,(w_1,w_2)4)2,
                (
                    (
                        (x_2,y_4)7,
                        ((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8
                    )6,
                    (w_3,(z_3,(t_1,t_2)14)13)12
                )5
            )1;
        """, sp_naming_function=get_species_name, format=1)

        species_tree = PhyloTree("(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;", format=1)

        rec = {
            **reconcile_leaves(gene_tree, species_tree),
            gene_tree & "1": species_tree,
            gene_tree & "2": species_tree & "XYZ",
            gene_tree & "3": species_tree & "XYZ",
            gene_tree & "4": species_tree & "W",
            gene_tree & "5": species_tree,
            gene_tree & "6": species_tree & "XYZ",
            gene_tree & "7": species_tree & "XY",
            gene_tree & "8": species_tree & "XYZ",
            gene_tree & "9": species_tree & "XY",
            gene_tree & "10": species_tree & "Y",
            gene_tree & "11": species_tree & "Y",
            gene_tree & "12": species_tree & "WT",
            gene_tree & "13": species_tree & "T",
            gene_tree & "14": species_tree & "T",
        }

        layout_info = layout(gene_tree, species_tree, rec)
        out = render_to_tikz(species_tree, rec, layout_info)

        self.assertEqual(out, textwrap.dedent(r"""
            \begin{tikzpicture}[
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
                    yshift=-10pt,
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
                    label={[label distance={0pt}]below:#1},
                },
                branch node/.style={
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={0.5pt},
                },
                speciation/.style={
                    branch node, circle,
                    minimum width={10pt},
                    minimum height={10pt},
                },
                duplication/.style={
                    branch node, rectangle,
                    minimum width={9pt},
                    minimum height={9pt},
                },
                horizontal gene transfer/.style={
                    branch node, diamond,
                    minimum width={10pt},
                    minimum height={10pt},
                },
            ]
            % species
            \draw[species border] (225,75) |- (365,20) -- (365,0);
            \draw[species border] (530,195) |- (420,20) -- (420,0);
            \draw[species border] (295,75) |- (365,75) -| (490,195);
            \draw[species border] (70,165) |- (225,95) -- (225,75);
            \draw[species border] (365,285) |- (295,95) -- (295,75);
            \draw[species border] (140,165) |- (225,165) -| (295,285);
            \draw[species border] (0,285) |- (70,185) -- (70,165);
            \draw[species border] (225,255) |- (140,185) -- (140,165);
            \draw[species border] (70,285) |- (70,255) -| (140,255);
            \draw[species border] (0,285) -- ([yshift=-16pt]0,305) -- node[species label] {X} ([yshift=-16pt]70,305) -- (70,285);
            \draw[species border] (140,255) -- ([yshift=-16pt]140,305) -- node[species label] {Y} ([yshift=-16pt]225,305) -- (225,255);
            \draw[species border] (295,285) -- ([yshift=-16pt]295,305) -- node[species label] {Z} ([yshift=-16pt]365,305) -- (365,285);
            \draw[species border] (420,270) |- (490,215) -- (490,195);
            \draw[species border] (585,255) |- (530,215) -- (530,195);
            \draw[species border] (490,270) |- (490,255) -| (530,255);
            \draw[species border] (420,270) -- ([yshift=-16pt]420,305) -- node[species label] {W} ([yshift=-16pt]490,305) -- (490,270);
            \draw[species border] (530,255) -- ([yshift=-16pt]530,305) -- node[species label] {T} ([yshift=-16pt]585,305) -- (585,255);
            % gene branches
            \draw[branch] (267.5,75) |- (385,40) -| (510,195);
            \draw[loss] (400,55) -- ++(20pt, 0);
            \draw[branch] (400,55) -| (245,75);
            \draw[branch] (392.5,20) -- (392.5,0);
            \draw[branch] (400,55) |- (392.5,20) -| (385,40);
            \draw[branch] (120,165) |- (245,115) -| (315,285);
            \draw[branch] (105,165) |- (260,130) -| (330,285);
            \draw[branch] (245,95) -- (245,75);
            \draw[branch] (245,115) |- (245,95);
            \draw[loss] (275,145) -- ++(20pt, 0);
            \draw[branch] (275,145) -| (90,165);
            \draw[branch] (267.5,95) -- (267.5,75);
            \draw[branch] (275,145) |- (267.5,95) -| (260,130);
            \draw[branch] (90,205) -- (90,165);
            \draw[branch] (35,285) |- (90,205) -| (160,255);
            \draw[branch] (105,220) -- (105,165);
            \draw[branch] (50,285) |- (105,220) -| (186.25,255);
            \draw[branch] (120,235) -- (120,165);
            \draw[loss] (120,235) -- ++(20pt, 0);
            \draw[branch] (120,235) -| (20,285);
            \draw[branch] (20,305) -- (20,285);
            \draw[branch] (35,305) -- (35,285);
            \draw[branch] (50,305) -- (50,285);
            \draw[branch] (160,305) -- (160,255);
            \draw[branch] (190,305) |- (197.5,290) -| (205,305);
            \draw[branch] (186.25,275) -- (186.25,255);
            \draw[branch] (175,305) |- (186.25,275) -| (197.5,290);
            \draw[branch] (315,305) -- (315,285);
            \draw[branch] (330,305) -- (330,285);
            \draw[branch] (345,305) -- (345,285);
            \draw[branch] (510,235) -- (510,195);
            \draw[branch] (470,270) |- (510,235) -| (557.5,255);
            \draw[branch] (470,305) -- (470,270);
            \draw[branch] (447.5,290) -- (447.5,270);
            \draw[branch] (440,305) |- (447.5,290) -| (455,305);
            \draw[branch] (550,305) |- (557.5,290) -| (565,305);
            \draw[branch] (557.5,275) -- (557.5,255);
            \draw[branch] (557.5,290) |- (557.5,275);
            % gene transfers
            \draw[transfer branch] (245,95) to[bend left=35] (447.5,270);
            \draw[transfer branch] (557.5,275) to[bend right=35] (345,285);
            % events
            \node[speciation] at (385,40) {};
            \node[duplication] at (392.5,20) {};
            \node[speciation] at (245,115) {};
            \node[speciation] at (260,130) {};
            \node[horizontal gene transfer] at (245,95) {};
            \node[duplication] at (267.5,95) {};
            \node[speciation] at (90,205) {};
            \node[speciation] at (105,220) {};
            \node[extant gene={\(x_1\)}] at (20,305) {};
            \node[extant gene={\(x_2\)}] at (35,305) {};
            \node[extant gene={\(x_3\)}] at (50,305) {};
            \node[extant gene={\(y_4\)}] at (160,305) {};
            \node[extant gene={\(y_1\)}] at (175,305) {};
            \node[extant gene={\(y_2\)}] at (190,305) {};
            \node[extant gene={\(y_3\)}] at (205,305) {};
            \node[duplication] at (197.5,290) {};
            \node[duplication] at (186.25,275) {};
            \node[extant gene={\(z_1\)}] at (315,305) {};
            \node[extant gene={\(z_2\)}] at (330,305) {};
            \node[extant gene={\(z_3\)}] at (345,305) {};
            \node[speciation] at (510,235) {};
            \node[extant gene={\(w_1\)}] at (440,305) {};
            \node[extant gene={\(w_2\)}] at (455,305) {};
            \node[extant gene={\(w_3\)}] at (470,305) {};
            \node[duplication] at (447.5,290) {};
            \node[extant gene={\(t_1\)}] at (550,305) {};
            \node[extant gene={\(t_2\)}] at (565,305) {};
            \node[duplication] at (557.5,290) {};
            \node[horizontal gene transfer] at (557.5,275) {};
            \end{tikzpicture}
        """).lstrip())
