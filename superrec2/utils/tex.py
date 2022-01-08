"""Bridge with the XeLaTeX compiler."""
import os
import subprocess
import shutil
import tempfile
from typing import Iterable, IO, List, NamedTuple
from .geometry import Size


class TeXError(Exception):
    """Raised when a TeX compiler returns an error."""

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code
        self.message = message


def xelatex_compile(source: str, dest: IO = None) -> str:
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
            check=False,
        )

        out = result.stdout.decode()

        if result.returncode != 0:
            raise TeXError(result.returncode, out)

        if dest is not None:
            with open(pdf_path, "rb") as output_file:
                shutil.copyfileobj(output_file, dest)

        return out


class MeasureBox(NamedTuple):
    """
    Dimensions of a TeX box.

    :attr width: total width of the box in points
    :attr height: height of the box above its baseline in points
    :attr depth: height of the box below its baseline in points
    """

    width: float
    height: float
    depth: float

    def overall_size(self) -> Size:
        """Get the overall dimensions of the box."""
        return Size(w=self.width, h=self.height + self.depth)


def measure(texts: Iterable[str], preamble="") -> List[MeasureBox]:
    """
    Measure dimensions of TeX boxes.

    :param texts: list of text strings to measure
    :param preamble: modules to load and macro definitions
    :raises TeXError: if a TeX error occurs
    :returns: list of measurements for each input string
    """
    src = (
        r"\documentclass{standalone}"
        "\n" + preamble + "\n"
        r"\newsavebox{\measurebox}"
        "\n"
        r"\scrollmode"
        "\n"
    )

    for text in texts:
        src += (
            rf"\savebox{{\measurebox}}{{{text}}}"
            "\n"
            r"\typeout{$$$"
            r"\the\wd\measurebox,\the\ht\measurebox,\the\dp\measurebox}"
            "\n"
        )

    src += r"\scrollmode" "\n" r"\begin{document}\end{document}" "\n"

    out = xelatex_compile(src)
    boxes = []

    for line in out.splitlines():
        if line.startswith("$$$"):
            boxes.append(
                MeasureBox(
                    *map(
                        lambda value: float(value.removesuffix("pt")),
                        line[3:].split(","),
                    )
                )
            )

    return boxes
