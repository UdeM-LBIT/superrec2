import importlib.resources
import subprocess
import os
import tempfile
from superrec2.render.model import Orientation


def test_fixtures():
    fixtures_path = importlib.resources.files(__package__) / "fixtures"

    with tempfile.TemporaryDirectory() as cache_dir:
        for fixture in fixtures_path.iterdir():
            run_fixture_test(fixture, cache_dir)


def run_fixture_test(fixture, cache_dir):
    input_path = fixture / "input.json"

    for orientation in Orientation:
        orientation_name = orientation.name.lower()
        output_path = fixture / f"output-{orientation_name}.tex"

        with open(output_path, "rb") as output_expect_file:
            output_expect = output_expect_file.read()
            result = subprocess.run(
                [
                    "python",
                    "-m",
                    "superrec2.cli",
                    "draw",
                    "--orientation",
                    orientation_name,
                    "--input",
                    input_path.as_posix(),
                ],
                env={
                    **os.environ,
                    "TECTONIC_CACHE_DIR": cache_dir,
                },
                capture_output=True,
                check=True,
            )
            assert result.stdout.decode() == output_expect.decode()
            assert result.stderr.decode() == ""
