import os
import sys
import tempfile
import pytest
from superrec2.cli.util import open_std


def test_open_std():
    with open_std("-", "r") as handle:
        assert handle == sys.stdin

    assert not sys.stdin.closed

    with open_std("-", "rb") as handle:
        assert handle == sys.stdin.buffer

    assert not sys.stdin.closed

    with open_std("-", "w") as handle:
        assert handle == sys.stdout

    assert not sys.stdout.closed

    # The following test does not work because unittest runner replaces sys.stdout
    # with a StringIO instance that does not support accessing a raw buffer
    #
    # with open_std("-", "wb") as handle:
    #    assert handle == sys.stdout.buffer
    #
    # assert not sys.stdout.closed

    with pytest.raises(RuntimeError) as error:
        with open_std("-", "r+") as handle:
            pass

    assert "Invalid mode 'r+' for standard stream" in str(error)

    path = tempfile.mkstemp()[1]

    with open_std(path, "w", encoding="utf8") as handle:
        assert handle.name == path
        assert not handle.closed

    assert handle.closed
    os.remove(path)
