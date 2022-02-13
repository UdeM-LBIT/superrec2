import unittest
import importlib.resources
import subprocess
from .model import Orientation


class TestDraw(unittest.TestCase):
    def test_fixtures(self):
        with importlib.resources.path(__package__, "fixtures") as fixtures_path:
            for fixture in fixtures_path.iterdir():
                self.do_test_fixture(fixture)

    def do_test_fixture(self, fixture):
        input_path = fixture / "input.json"

        for orientation in Orientation:
            orientation_name = orientation.name.lower()
            output_path = fixture / f"output-{orientation_name}.tex"

            with open(output_path, "rb") as output_expect_file:
                output_expect = output_expect_file.read()
                result = subprocess.run(
                    [
                        "./draw.py",
                        "--orientation",
                        orientation_name,
                        "--input",
                        input_path
                    ],
                    capture_output=True,
                    check=True,
                )
                self.assertEqual(
                    result.stdout.decode(),
                    output_expect.decode(),
                    f"{fixture.name} ({orientation}): valid output"
                )
                self.assertEqual(
                    result.stderr.decode(),
                    "",
                    f"{fixture.name} ({orientation}): no error"
                )
