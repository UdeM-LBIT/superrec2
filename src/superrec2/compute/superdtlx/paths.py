from typing import TypeVar
from sowing.indexed import IndexedTree
from ...utils.algebras import Semiring
from ...model.history import Host, Codiverge, Diverge, Gain, Loss, Extant
from .contents import Contents, EXTRA_CONTENTS


T = TypeVar("T")


def make_codiv_path(
    start_host: str,
    end_host: str,
    contents: Contents,
    host_index: IndexedTree[Host, None],
    structure: type[Semiring[T]],
    path: Semiring[T],
) -> Semiring[T]:
    """
    Try to link two hosts using codivergence and loss events.

    :param start_host: origin host
    :param end_host: end host
    :param contents: event contents
    :param host_index: indexed host tree
    :param structure: event value semiring
    :param path: path below
    :returns: possible paths, if any
    """
    if not host_index.is_ancestor_of(start_host, end_host):
        return structure.null()

    host = host_index[end_host]
    target = host_index[start_host]

    while host != target:
        last = host

        host = host.up()
        host_data = host.node.data

        other = last.sibling()
        other_data = other.node.data

        subpath = structure.make(Codiverge(host=host_data.name, contents=contents))

        if other_data.sampled:
            # Add loss of associate in sampled host
            event = Loss(host=other_data.name, contents=contents, segment=contents)
        else:
            # Add extant associate in unsampled host
            event = Extant(host=other_data.name, contents=contents)

        path = subpath * structure.make(event) * path

    return path


def make_transfer_path(
    start_host: str,
    end_host: str,
    start_contents: Contents,
    end_contents: Contents,
    host_index: IndexedTree[Host, None],
    structure: type[Semiring[T]],
    path: Semiring[T],
) -> Semiring[T]:
    """
    Try to link two hosts using codivergences, losses,
    and exactly one transfer event.

    :param start_host: origin host
    :param end_host: end host
    :param contents: event contents
    :param host_index: indexed host tree
    :param structure: event value semiring
    :param path: path below
    :returns: possible paths, if any
    """
    if host_index.is_ancestor_of(end_host, start_host):
        return structure.null()

    if not (end_contents <= start_contents):
        return structure.null()

    host = host_index[start_host]

    # If needed, first look for the closest separate host to start a transfer from
    if host_index.is_comparable(start_host, end_host):
        left = host.down(0)
        right = host.down(1)

        sep = right if host_index.is_ancestor_of(left, end_host) else left
        sep_host = sep.node.data.name

        subpath = make_transfer_path(
            sep_host,
            end_host,
            start_contents,
            end_contents,
            host_index,
            structure,
            path,
        )
        return make_codiv_path(
            start_host,
            sep_host,
            start_contents,
            host_index,
            structure,
            subpath,
        )

    copy = structure.make(
        Diverge(
            host=start_host,
            contents=start_contents,
            segment=end_contents,
            cut=False,
            transfer=True,
            result=1,
        )
    )

    if start_contents == end_contents:
        # Complete cut transfer
        cut = structure.make(
            Diverge(
                host=start_host,
                contents=start_contents,
                segment=end_contents,
                cut=True,
                transfer=True,
                result=0,
            )
        )
    else:
        cut = structure.make(
            Diverge(
                host=start_host,
                contents=start_contents,
                segment=end_contents,
                cut=True,
                transfer=True,
                result=1,
            )
        )

    sampled = host.node.data.sampled

    if sampled:
        # Lose remaining contents
        copy *= structure.make(
            Loss(
                host=start_host,
                contents=start_contents,
                segment=start_contents,
            )
        )

        if start_contents != end_contents:
            cut *= structure.make(
                Loss(
                    host=start_host,
                    contents=start_contents - end_contents,
                    segment=start_contents - end_contents,
                )
            )
    else:
        # Keep as extant in unsampled host
        copy *= structure.make(
            Extant(
                host=start_host,
                contents=start_contents,
            )
        )

        if start_contents != end_contents:
            cut *= structure.make(
                Extant(
                    host=start_host,
                    contents=start_contents - end_contents,
                )
            )

    return (copy + cut) * path


def make_gain_path(
    host: str,
    start_contents: Contents,
    end_contents: Contents,
    structure: type[Semiring[T]],
    path: Semiring[T],
) -> Semiring[T]:
    """
    Add a gain event, if required, to transition between two gene contents sets.

    :param host: event host
    :param start_contents: root contents
    :param end_contents: target contents with genes to gain
    :param structure: event value semiring
    :param path: path below
    :returns: extended path
    """
    to_gain = end_contents - start_contents

    if to_gain:
        path = (
            structure.make(
                Gain(
                    host=host,
                    contents=start_contents,
                    gained=to_gain,
                )
            )
            * path
        )

    return path


def make_path(
    start_host: str,
    end_host: str,
    start_contents: Contents,
    end_contents: Contents,
    host_index: IndexedTree[Host, None],
    structure: type[Semiring[T]],
    path: Semiring[T],
) -> Semiring[T]:
    """
    Link two hosts and contents using a compressible path.

    :param start_host: origin host
    :param end_host: end host
    :param start_contents: origin contents
    :param end_contents: end contents
    :param host_index: indexed host tree
    :param structure: event value semiring
    :param path: path below
    :returns: possible paths, if any
    """
    contents = end_contents

    if (
        start_contents <= contents
        and EXTRA_CONTENTS not in start_contents
        and EXTRA_CONTENTS in contents
    ):
        # No extra contents can be sent towards the child
        return structure.null()

    # Add gain events at the end, if needed
    without_gains = contents & (start_contents | {EXTRA_CONTENTS})
    path = make_gain_path(
        end_host,
        without_gains,
        contents,
        structure,
        path,
    )

    contents = without_gains

    # Try reaching using no transfers
    to_lose = start_contents - contents

    if to_lose and EXTRA_CONTENTS not in contents:
        ## TODO: Loss using duplication or cut in unsampled species ##
        ##############################################################
        codiv_path = (
            structure.make(
                Loss(
                    host=end_host,
                    contents=start_contents,
                    segment=to_lose,
                ),
            )
            * path
        )
    else:
        codiv_path = path

    codiv_path = make_codiv_path(
        start_host,
        end_host,
        start_contents,
        host_index,
        structure,
        codiv_path,
    )

    # Try reaching using one transfer
    transfer_path = make_transfer_path(
        start_host,
        end_host,
        start_contents,
        contents,
        host_index,
        structure,
        path,
    )

    ## TODO: Reach host using two transfers
    #######################################

    return codiv_path + transfer_path
