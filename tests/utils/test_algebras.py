import pytest
from math import inf
from typing import TypeVar, Sequence, NamedTuple
from collections import defaultdict
from immutables import Map
from superrec2.utils.algebras import (
    Structure,
    UnitMagma,
    Counter,
    MinPlus,
    MaxPlus,
    Viterbi,
    Boolean,
    pareto_of,
    generator_of,
    vector,
)
from sowing.node import Node


T = TypeVar("T")


def test_metaclasses():
    assert MinPlus._pool is not MaxPlus._pool
    assert MinPlus._pool is not Boolean._pool

    assert MinPlus.zero is MinPlus(MinPlus._zero)
    assert MinPlus.one is MinPlus(MinPlus._one)

    assert Boolean.zero is Boolean(Boolean._zero)
    assert Boolean.one is Boolean(Boolean._one)


def test_selection():
    minplus_counter = MinPlus * Counter

    assert minplus_counter.__name__ == "(MinPlus * Counter)"
    assert minplus_counter.zero == minplus_counter((inf, 0))
    assert minplus_counter.one == minplus_counter((0, 1))
    assert minplus_counter((2, 5)) + minplus_counter((3, 2)) == minplus_counter((2, 5))
    assert minplus_counter((3, 5)) + minplus_counter((2, 2)) == minplus_counter((2, 2))
    assert minplus_counter((2, 5)) + minplus_counter((2, 2)) == minplus_counter((2, 7))
    assert minplus_counter((2, 5)) * minplus_counter((3, 2)) == minplus_counter((5, 10))

    with pytest.raises(AssertionError, match="left argument must be ordered"):
        Counter * MinPlus


def edit_distance(
    word1: Sequence[str],
    word2: Sequence[str],
    structure: Structure,
) -> T:
    table = defaultdict(lambda: structure.zero)
    table[(0, 0)] = structure.one

    n = len(word1)
    m = len(word2)

    for i in range(n + 1):
        for j in range(m + 1):
            letter1 = word1[i - 1] if i >= 1 else None
            letter2 = word2[j - 1] if j >= 1 else None

            if i >= 1 and j >= 1:
                table[(i, j)] += table[(i - 1, j - 1)] * structure(letter1, letter2)

            if i >= 1:
                table[(i, j)] += table[(i - 1, j)] * structure(letter1, None)

            if j >= 1:
                table[(i, j)] += table[(i, j - 1)] * structure(None, letter2)

    return table[(n, m)].value


@vector
class EditVector(NamedTuple):
    insert: int = 0
    delete: int = 0
    change: int = 0


def test_edit_distance():
    words = ("elephant", "relevnat")

    # Count the total number of possible alignments
    alignment_count = Structure(Counter, lambda source, target: 1)
    assert edit_distance(*words, alignment_count) == 265729

    # Compute the minimum cost of any alignment
    edit_cost = Structure(MinPlus, lambda source, target: 1 if source != target else 0)
    assert edit_distance(*words, edit_cost) == 4

    # Compute the maximum score of any alignment
    edit_score = Structure(
        MaxPlus, lambda source, target: -1 if source != target else 1
    )
    assert edit_distance(*words, edit_score) == 1

    # Compute the set of Pareto-optimal costs of any alignment
    def vector_morphism(source, target):
        if source is None:
            return frozenset({EditVector(insert=1)})
        elif target is None:
            return frozenset({EditVector(delete=1)})
        elif source != target:
            return frozenset({EditVector(change=1)})
        else:
            return frozenset({EditVector()})

    EditPareto = pareto_of(EditVector())
    edit_pareto = Structure(EditPareto, vector_morphism)

    assert edit_distance(*words, edit_pareto) == frozenset(
        {
            EditVector(insert=0, delete=0, change=7),
            EditVector(insert=1, delete=1, change=2),
            EditVector(insert=2, delete=2, change=1),
            EditVector(insert=3, delete=3, change=0),
        }
    )

    # Count the number of minimum-cost alignments
    assert edit_distance(*words, edit_cost * alignment_count) == (4, 1)

    # Count the number of maximum-score alignments
    assert edit_distance(*words, edit_score * alignment_count) == (1, 1)

    # Count the number of each type of Pareto-optimal alignment
    assert edit_distance(*words, edit_pareto @ alignment_count) == Map(
        {
            EditVector(insert=0, delete=0, change=7): 1,
            EditVector(insert=1, delete=1, change=2): 1,
            EditVector(insert=2, delete=2, change=1): 9,
            EditVector(insert=3, delete=3, change=0): 10,
        }
    )

    # Generate all minimum-cost alignments
    class AlignBuilder:
        def __init__(self, word1, word2):
            self.value = (word1, word2)

        def __eq__(self, other):
            return self.value == other.value

        def __hash__(self):
            return hash(self.value)

        def __mul__(align1, align2):
            return AlignBuilder(
                align1.value[0] + align2.value[0], align1.value[1] + align2.value[1]
            )

    AlignGenerator = generator_of(AlignBuilder((), ()))
    align_generator = Structure(
        AlignGenerator,
        lambda source, target: frozenset({AlignBuilder((source,), (target,))}),
    )

    assert edit_distance(*words, edit_cost * align_generator) == (
        4,
        frozenset(
            {
                AlignBuilder(
                    (None, "e", "l", "e", "p", "h", "a", "n", "t"),
                    ("r", "e", "l", "e", "v", "n", "a", None, "t"),
                ),
            }
        ),
    )

    # Generate all maximum-score alignments
    assert edit_distance(*words, edit_score * align_generator) == (
        1,
        frozenset(
            {
                AlignBuilder(
                    (None, "e", "l", "e", "p", "h", "a", "n", "t"),
                    ("r", "e", "l", "e", "v", "n", "a", None, "t"),
                ),
            }
        ),
    )


def parse_grammar(
    grammar: tuple[tuple[str, tuple[str, str] | str, float]],
    word: Sequence[str],
    structure: Structure,
) -> T:
    table = defaultdict(lambda: structure.zero)
    n = len(word)

    for start in range(len(word)):
        for rule in grammar:
            head, tail, *_ = rule
            if tail == word[start]:
                table[(start, 1, head)] = structure(rule)

    for size in range(2, len(word) + 1):
        for start in range(len(word) - size + 1):
            for cut in range(1, size):
                for rule in grammar:
                    head, tail, *_ = rule
                    if isinstance(tail, tuple):
                        left, right = tail
                        table[(start, size, head)] += (
                            structure(rule)
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
    parsable = Structure(Boolean, lambda rule: True)
    assert parse_grammar(grammar, word, parsable)
    assert not parse_grammar(grammar, inv_word, parsable)

    # Compute the best parsing probability of a sentence
    best_prob = Structure(Viterbi, lambda rule: rule[2])
    assert parse_grammar(grammar, word, best_prob) == 0.24

    # Count the number of possible parse trees for a sentence
    count_parses = Structure(Counter, lambda rule: 1)
    assert parse_grammar(grammar, word, count_parses) == 2

    # Generate all possible parse trees for a sentence
    class ParseBuilder(UnitMagma[Node[str, None]]):
        _one = Node()

        def _mul(node1, node2):
            if node1.data is None:
                return node2

            if node2.data is None:
                return node1

            return node1.add(node2)

    def parse_morphism(rule):
        head, tail, *_ = rule

        if isinstance(tail, str):
            return frozenset({ParseBuilder(Node(head).add(Node(tail)))})
        else:
            return frozenset({ParseBuilder(Node(head))})

    ParseGenerator = generator_of(ParseBuilder(Node()))
    parse_generator = Structure(ParseGenerator, parse_morphism)

    assert parse_grammar(grammar, word, parse_generator) == frozenset(
        {
            ParseBuilder(
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
            ParseBuilder(
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
    assert parse_grammar(grammar, word, best_prob * parse_generator) == (
        0.24,
        frozenset(
            {
                ParseBuilder(
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
        ),
    )
