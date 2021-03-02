import os
import subprocess
import shutil
import tempfile
from typing import List, NamedTuple


class TeXError(Exception):
    """Raised when a TeX compiler returns an error."""
    def __init__(self, code, message):
        self.code = code
        self.message = message


def xelatex_compile(source: str, dest: str = None) -> str:
    """
    Compile a LaTeX source file.

    :param source: string containing the full LaTeX source
    :param dest: if not None, will store the generated PDF
        at the given location
    :raises TeXError: if a TeX error occurs
    :returns: output logs from the TeX compiler
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        source_path = os.path.join(tmpdir, "input.tex")
        pdf_path = os.path.join(tmpdir, "input.pdf")

        with open(source_path, "w") as source_file:
            source_file.write(source)

        result = subprocess.run(
            ["xelatex", "-interaction", "batchmode", source_path],
            cwd=tmpdir,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
        )

        out = result.stdout.decode()

        if result.returncode != 0:
            raise TeXError(result.returncode, out)

        if dest is not None:
            shutil.move(pdf_path, dest)

        return out


class MeasureBox(NamedTuple):
    """
    Dimensions of a TeX box.

    :attr width: total width of the box in points
    :attr height: height of the box above its baseline in points
    :attr depth: height of the box below its baseline in points
    """
    width: int
    height: int
    depth: int


def measure(texts: List[str]) -> List[MeasureBox]:
    """
    Measure dimensions of text pieces.

    :param texts: list of text strings to measure
    :raises TeXError: if a TeX error occurs
    :returns: list of measurements for each input string
    """
    src = (
        r"\documentclass{standalone}" "\n"
        r"\newsavebox{\measurebox}" "\n"
        r"\scrollmode" "\n"
    )

    for text in texts:
        src += (
            rf"\savebox{{\measurebox}}{{{text}}}" "\n"
            r"\typeout{$$$"
            r"\the\wd\measurebox,\the\ht\measurebox,\the\dp\measurebox}" "\n"
        )

    src += (
        r"\scrollmode" "\n"
        r"\begin{document}\end{document}" "\n"
    )

    out = xelatex_compile(src)
    boxes = []

    for line in out.splitlines():
        if line.startswith("$$$"):
            boxes.append(MeasureBox(*line[3:].split(",")))

    return boxes
