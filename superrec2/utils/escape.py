def escape(name: str, escape_chars: str) -> str:
    """
    Escape special characters from a name.

    :param name: original name
    :param escape_chars: set of characters to escape
    :returns: name with escaped characters
    """
    for escape_char in escape_chars:
        name = name.replace(escape_char, "\\" + escape_char)

    return name


def unescape(name: str) -> str:
    """
    Unescape special characters from a name.

    Replaces each \X sequence with X (where X is any character).

    :param name: escaped name
    :returns: original name
    """
    escaped = False
    result = ""

    for char in name:
        if escaped:
            result += char
            escaped = False
        elif char == "\\":
            escaped = True
        else:
            result += char

    return result
