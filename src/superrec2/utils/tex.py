"""Bridge with the XeLaTeX compiler."""
import os
import subprocess
import shutil
import tempfile
from typing import Iterable, IO, List
from .geometry import Rect, Position, Size


class TeXError(Exception):
    """Raised when a TeX compiler returns an error."""

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code
        self.message = message


def _xelatex_compile(source: str, dest: IO = None) -> str:
    with tempfile.TemporaryDirectory() as tmpdir:
        source_path = os.path.join(tmpdir, "input.tex")
        pdf_path = os.path.join(tmpdir, "input.pdf")

        with open(source_path, "w", encoding="utf8") as source_file:
            source_file.write(source)

        result = subprocess.run(
            ["xelatex", "-interaction", "batchmode", source_path],
            cwd=tmpdir,
            encoding="utf8",
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )

        if result.returncode != 0:
            raise TeXError(result.returncode, result.stdout)

        if dest is not None:
            with open(pdf_path, "rb") as output_file:
                shutil.copyfileobj(output_file, dest)

        return result.stdout


def _tectonic_compile(source: str, dest: IO = None) -> str:
    with tempfile.TemporaryDirectory() as tmpdir:
        commandline = [
            "tectonic",
            "-",
            "--outdir",
            tmpdir,
            "--print",
            "--untrusted",
            "--chatter",
            "minimal",
            "--reruns",
            "0",
        ]

        if dest is None:
            commandline += [
                "--outfmt",
                "aux",
            ]

        result = subprocess.run(
            commandline,
            input=source,
            encoding="utf8",
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )

        if result.returncode != 0:
            raise TeXError(result.returncode, result.stdout)

        if dest is not None:
            with open(tmpdir + "/texput.pdf", "rb") as output_file:
                shutil.copyfileobj(output_file, dest)

        return result.stdout


def tex_compile(source: str, dest: IO = None) -> str:
    """
    Compile a TeX source file.

    :param source: string containing the full LaTeX source
    :param dest: if not None, will store the generated PDF
        at the given location
    :raises TeXError: if a TeX error occurs
    :returns: output logs from the TeX compiler
    """
    if shutil.which("tectonic") is not None:
        return _tectonic_compile(source, dest)

    if shutil.which("xelatex") is not None:
        return _xelatex_compile(source, dest)

    raise TeXError("No xelatex or tectonic binary found")


def escape(text: str) -> str:
    """Escape a string for inclusion in a TeX document."""
    return text.replace("\\", "\\\\").replace(r"_", r"\_")


def measure_tikz(nodes: Iterable[str], preamble="") -> List[Rect]:
    """
    Measure dimensions of TikZ nodes.

    :param nodes: list of TikZ nodes to measure
    :param preamble: modules to load and macro definitions
    :raises TeXError: if a TeX error occurs
    :returns: list of measurements for each input string
    """
    log_prefix = ">>>"
    src = (
        r"\documentclass{standalone}",
        preamble,
        r"\newcommand{\tikzmeasure}[1]{\tikz{"
        "#1"
        r"\path (current bounding box.north west);"
        r"\pgfgetlastxy{\boxleft}{\boxtop}"
        r"\path (current bounding box.south east);"
        r"\pgfgetlastxy{\boxright}{\boxbottom}"
        rf"\typeout{{{log_prefix}\boxleft,\boxtop,\boxright,\boxbottom}}"
        "}}",
        r"\begin{document}",
    )

    for node in nodes:
        src += (rf"\tikzmeasure{{{node}}}",)

    src += (r"\end{document}",)

    out = tex_compile("\n".join(src))
    rects = []

    for line in out.splitlines():
        if line.startswith(log_prefix):
            left, top, right, bottom = map(
                lambda value: float(value.removesuffix("pt")),
                line.removeprefix(log_prefix).split(","),
            )
            rects.append(
                Rect(
                    Position(left, -top),
                    Size(right - left, top - bottom),
                )
            )

    return rects
