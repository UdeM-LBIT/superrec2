"""Text wrapping utilities."""
import textwrap
from typing import List


def _wrap_badness(lines: List[str]) -> int:
    """
    Measures the balance badness of a wrapping solution, as the
    squared difference of each line length to the maximum line length.

    :param lines: list of lines resulting from the wrapping
    :returns: badness measure
    """
    max_len = max(len(line) for line in lines)
    return sum((len(line) - max_len) ** 2 for line in lines)


def balanced_wrap(text: str, width: int) -> str:
    """
    Wrap a text to not exceed the given width, while making the lines as
    balanced as possible (i.e. each having a number of characters as close
    as possible to being equal to the others).

    :param text: input text
    :param width: target wrap width
    :returns: best solution
    """
    if not text:
        return ""

    best_result = textwrap.wrap(text, width, break_long_words=False)
    best_badness = _wrap_badness(best_result)
    line_count = len(best_result)

    # Reduce wrap width looking for the most balanced solution
    # while keeping the same number of lines (naive algorithm)
    while width > 1:
        width -= 1
        next_result = textwrap.wrap(text, width, break_long_words=False)

        if len(next_result) != line_count:
            break

        next_badness = _wrap_badness(next_result)

        if next_badness < best_badness:
            best_result = next_result
            best_badness = next_badness

    return "\n".join(best_result)
