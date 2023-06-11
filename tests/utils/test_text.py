from superrec2.utils.text import balanced_wrap


def test_balanced_wrap():
    assert (
        balanced_wrap("this is a short test with short words", 32)
        == "this is a short test\nwith short words"
    )
    assert (
        balanced_wrap("elongated beautiful unbreakable wording", 8)
        == "elongated\nbeautiful\nunbreakable\nwording"
    )
    assert balanced_wrap("", 8) == ""
