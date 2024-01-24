from typing import TypeVar, Sequence, NamedTuple
from collections import defaultdict
from immutables import Map
from superrec2.utils.algebras import (
    Semiring,
    tuple_ordered_magma,
    min_plus,
    max_plus,
    pareto,
    viterbi,
    count,
    boolean,
    make_unit_magma,
    make_single_selector,
    make_multiple_selector,
    make_generator,
)
from sowing.node import Node


T = TypeVar("T")


def edit_distance(
    word1: Sequence[str],
    word2: Sequence[str],
    structure: type[Semiring[T]],
) -> T:
    table = defaultdict(structure.null)
    table[(0, 0)] = structure.unit()

    n = len(word1)
    m = len(word2)

    for i in range(n + 1):
        for j in range(m + 1):
            letter1 = word1[i - 1] if i >= 1 else None
            letter2 = word2[j - 1] if j >= 1 else None

            if i >= 1 and j >= 1:
                table[(i, j)] += table[(i - 1, j - 1)] * structure.make(
                    letter1, letter2
                )

            if i >= 1:
                table[(i, j)] += table[(i - 1, j)] * structure.make(letter1, None)

            if j >= 1:
                table[(i, j)] += table[(i, j - 1)] * structure.make(None, letter2)

    return table[(n, m)].value


@tuple_ordered_magma
class EditVector(NamedTuple):
    insert: int = 0
    delete: int = 0
    change: int = 0

    @classmethod
    def make(cls, source, target):
        if source is None:
            return EditVector(insert=1)
        elif target is None:
            return EditVector(delete=1)
        elif source != target:
            return EditVector(change=1)
        else:
            return EditVector()


def test_edit_distance():
    words = ("elephant", "relevnat")

    # Count the total number of possible alignments
    alignment_count = count("alignment_count", lambda source, target: 1)
    assert edit_distance(*words, alignment_count) == 265729

    # Compute the minimum cost of any alignment
    edit_cost = min_plus(
        "edit_cost", lambda source, target: 1 if source != target else 0
    )
    assert edit_distance(*words, edit_cost) == 4

    # Compute the maximum score of any alignment
    edit_score = max_plus(
        "edit_score", lambda source, target: -1 if source != target else 1
    )
    assert edit_distance(*words, edit_score) == 1

    # Compute the set of Pareto-optimal costs of any alignment
    edit_pareto = pareto("edit_pareto", EditVector)
    assert edit_distance(*words, edit_pareto) == frozenset(
        {
            EditVector(insert=0, delete=0, change=7),
            EditVector(insert=1, delete=1, change=2),
            EditVector(insert=2, delete=2, change=1),
            EditVector(insert=3, delete=3, change=0),
        }
    )

    # Count the number of minimum-cost alignments
    count_min_cost = make_single_selector("count_min_cost", edit_cost, alignment_count)
    min_count = edit_distance(*words, count_min_cost)
    assert min_count.key.value == 4
    assert min_count.value.value == 1

    # Count the number of maximum-score alignments
    count_max_score = make_single_selector(
        "count_max_score", edit_score, alignment_count
    )
    max_count = edit_distance(*words, count_max_score)
    assert max_count.key.value == 1
    assert max_count.value.value == 1

    # Count the number of each type of Pareto-optimal alignment
    count_min_pareto = make_multiple_selector(
        "count_min_pareto", edit_pareto, alignment_count
    )
    assert edit_distance(*words, count_min_pareto) == Map(
        {
            EditVector(insert=0, delete=0, change=7): alignment_count(1),
            EditVector(insert=1, delete=1, change=2): alignment_count(1),
            EditVector(insert=2, delete=2, change=1): alignment_count(9),
            EditVector(insert=3, delete=3, change=0): alignment_count(10),
        }
    )

    # Generate all minimum-cost alignments
    alignment_builder = make_unit_magma(
        "alignment_builder",
        unit=((), ()),
        mul=lambda align1, align2: (align1[0] + align2[0], align1[1] + align2[1]),
        make=lambda letter1, letter2: ((letter1,), (letter2,)),
    )
    alignment_generator = make_generator("alignment_generator", alignment_builder)
    min_alignment_selector = make_single_selector(
        "min_alignment_selector",
        edit_cost,
        alignment_generator,
    )

    solutions = edit_distance(*words, min_alignment_selector)
    assert solutions.key.value == 4
    assert solutions.value.value == frozenset(
        {
            alignment_builder(
                (
                    (None, "e", "l", "e", "p", "h", "a", "n", "t"),
                    ("r", "e", "l", "e", "v", "n", "a", None, "t"),
                )
            ),
        }
    )

    # Generate all maximum-score alignments
    max_alignment_selector = make_single_selector(
        "max_alignment_selector",
        edit_score,
        alignment_generator,
    )
    solutions = edit_distance(*words, max_alignment_selector)
    assert solutions.key.value == 1
    assert solutions.value.value == frozenset(
        {
            alignment_builder(
                (
                    (None, "e", "l", "e", "p", "h", "a", "n", "t"),
                    ("r", "e", "l", "e", "v", "n", "a", None, "t"),
                )
            ),
        }
    )


def parse_grammar(
    grammar: tuple[tuple[str, tuple[str, str] | str, float]],
    word: Sequence[str],
    structure: type[Semiring[T]],
) -> T:
    table = defaultdict(structure.null)
    n = len(word)

    for start in range(len(word)):
        for rule in grammar:
            head, tail, *_ = rule
            if tail == word[start]:
                table[(start, 1, head)] = structure.make(rule)

    for size in range(2, len(word) + 1):
        for start in range(len(word) - size + 1):
            for cut in range(1, size):
                for rule in grammar:
                    head, tail, *_ = rule
                    if isinstance(tail, tuple):
                        left, right = tail
                        table[(start, size, head)] += (
                            structure.make(rule)
                            * table[(start, cut, left)]
                            * table[(start + cut, size - cut, right)]
                        )

    return table[(0, n, "S")].value


def test_parse_grammar():
    grammar = (
        ("S", ("NP", "VP"), 1),
        ("NP", ("NP", "PP"), 0.4),
        ("PP", ("P", "NP"), 1),
        ("VP", ("VBD", "NP"), 0.4),
        ("VP", ("VP", "PP"), 0.6),
        ("NP", "chopsticks", 1),
        ("NP", "i", 1),
        ("NP", "sushi", 1),
        ("P", "with", 1),
        ("VBD", "ate", 1),
    )
    word = ("i", "ate", "sushi", "with", "chopsticks")
    inv_word = ("chopsticks", "i", "sushi", "ate", "with")

    # Check whether a sentence can be generated by the grammar
    parsable = boolean("parsable", lambda rule: True)
    assert parse_grammar(grammar, word, parsable)
    assert not parse_grammar(grammar, inv_word, parsable)

    # Compute the best parsing probability of a sentence
    best_prob = viterbi("best_prob", lambda rule: rule[2])
    assert parse_grammar(grammar, word, best_prob) == 0.24

    # Count the number of possible parse trees for a sentence
    count_parses = count("count_parses", lambda rule: 1)
    assert parse_grammar(grammar, word, count_parses) == 2

    # Generate all possible parse trees for a sentence
    def join_trees(node1, node2):
        if node1.data is None:
            return node2

        if node2.data is None:
            return node1

        return node1.add(node2)

    def make_parse_tree(rule):
        head, tail, *_ = rule

        if isinstance(tail, str):
            return Node(head).add(Node(tail))
        else:
            return Node(head)

    parse_tree_builder = make_unit_magma(
        "parse_tree_builder",
        unit=Node(),
        mul=join_trees,
        make=make_parse_tree,
    )
    parse_tree_generator = make_generator("parse_tree_generator", parse_tree_builder)

    assert parse_grammar(grammar, word, parse_tree_generator) == frozenset(
        {
            parse_tree_builder(
                Node("S")
                .add(Node("NP").add(Node("i")))
                .add(
                    Node("VP")
                    .add(
                        Node("VP")
                        .add(Node("VBD").add(Node("ate")))
                        .add(Node("NP").add(Node("sushi")))
                    )
                    .add(
                        Node("PP")
                        .add(Node("P").add(Node("with")))
                        .add(Node("NP").add(Node("chopsticks")))
                    )
                )
            ),
            parse_tree_builder(
                Node("S")
                .add(Node("NP").add(Node("i")))
                .add(
                    Node("VP")
                    .add(Node("VBD").add(Node("ate")))
                    .add(
                        Node("NP")
                        .add(Node("NP").add(Node("sushi")))
                        .add(
                            Node("PP")
                            .add(Node("P").add(Node("with")))
                            .add(Node("NP").add(Node("chopsticks")))
                        )
                    )
                )
            ),
        }
    )

    # Generate the most probable parse trees for a sentence
    best_parse_tree_selector = make_single_selector(
        "best_parse_tree_selector",
        best_prob,
        parse_tree_generator,
    )

    solutions = parse_grammar(grammar, word, best_parse_tree_selector)
    assert solutions.key.value == 0.24
    assert solutions.value.value == frozenset(
        {
            parse_tree_builder(
                Node("S")
                .add(Node("NP").add(Node("i")))
                .add(
                    Node("VP")
                    .add(
                        Node("VP")
                        .add(Node("VBD").add(Node("ate")))
                        .add(Node("NP").add(Node("sushi")))
                    )
                    .add(
                        Node("PP")
                        .add(Node("P").add(Node("with")))
                        .add(Node("NP").add(Node("chopsticks")))
                    )
                )
            ),
        }
    )
