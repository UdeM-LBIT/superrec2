"""Compute a minimum-cost (super-)reconciliation of two trees."""

import argparse
import json
import inspect
from ast import literal_eval
from .util import add_arg_input, add_arg_output
from ..model.history import Reconciliation, History, graft_unsampled_hosts
from ..compute import synesth

from sowing.node import Node
from collections import Counter


methods = {}


def register_method(method):
    pretty_name = method.__name__.replace("_", "-")
    methods[pretty_name] = method


@register_method
def min_cost_single(setting, costs, output):
    """
    Report a single arbitrary minimum-cost history along with the total number
    of co-optimal solutions.
    """
    cost, count, history = synesth.min_cost_single(setting, costs)
    print(f"cost={cost}", file=output)
    print(f"count={count}", file=output)
    json.dump(history.to_mapping(), output)


@register_method
def min_cost_all(setting, costs, output):
    """Report all minimum-cost histories (may be exponentially slow!)."""
    cost, histories = synesth.min_cost_all(setting, costs)
    print(f"cost={cost}", file=output)
    print(f"count={len(histories)}", file=output)

    for history in histories:
        json.dump(history.to_mapping(), output)
        print(file=output)

@register_method
def general_solution(setting, _, output):
    """Report a representative solution in each class."""
    def get_sub_history(hist) :
        data, edge = get_data_edges(hist)
        if data.apparent:
            return synesth.HistoryBuilder(Node("continue")), hist
        else:
            h = synesth.HistoryBuilder(Node(data))
            next = None
            for i in edge:
                t, n = get_sub_history(i)
                h = h * t
                if n != None:
                    next = n
            return h, next

    def calc_freq(hist, frequence):
        h, n = get_sub_history(hist)

        if n == None :
            return frequence

        temp_n = n.node.data
        if temp_n not in frequence:
            frequence[temp_n] = {}
        if h not in frequence[temp_n]:
            frequence[temp_n][h] = 0
        frequence[temp_n][h] += 1

        for i in n.node.edges :
            frequence = calc_freq(i, frequence)

        return frequence

    def add_2_history(genral, node, hist):
        data, edge = get_data_edges(genral)
        if data == 'continue':
            newNode = synesth.HistoryBuilder(Node(node))
            for i in hist:
                newNode = newNode * i
            return newNode
        else:
            newNode = synesth.HistoryBuilder(Node(data))
            for i in edge:
                newNode = newNode * add_2_history(i, node, hist)
            return newNode


    def construct_representative(genral, node, edge, frequence):
        hist = []
        for i in edge:
            n, h = get_most_freq(i, frequence)
            _, edge1 = get_data_edges(i)
            hist.append(construct_representative(h, n, edge1, frequence))
        genral = add_2_history(genral, node, hist)

        return genral

    def get_most_freq(classe, frequence):
        data, edge = get_data_edges(classe)

        max = [-1, None]
        for i in frequence[data]:
            if frequence[data][i] > max[0]:
                max = [frequence[data][i], i]

        return data, max[1]

    def before_construct_representative(hist, frequence):
        n, h = get_most_freq(hist, frequence)
        edge = hist.node.edges if 'node' in dir(hist) else hist.value.edges
        return construct_representative(h, n, edge, frequence)

    def get_data_edges(hist):
        d = dir(hist) 
        if 'node' in d :
            return hist.node.data, hist.node.edges
        elif 'value' in d :
            return hist.value.data, hist.value.edges
        else :
            return hist.data, hist.edges


    import time

    t0 = time.time()
    result_number = synesth.solve(setting, synesth.partial_history_generator @ (synesth.event_counts_pareto @ synesth.history_generator)).value
    t1 = time.time()
    print(f"time run : {t1-t0}")
    nbClass = []
    nbClassOrigin = []

    for h in result_number:
        t = time.time()
        print(f"Classe :\n{h.value}")
        frequence = {}
        list_pareto = set()

        s = 0
        for h2 in result_number[h]:
            s += len(result_number[h][h2])
            for h3 in result_number[h][h2]:
                list_pareto.add(synesth.propagate_contents(History(setting.host_tree, h3.value).prune_unsampled()))

        nbClassOrigin.append(s)
        nbClass.append(len(list_pareto))

        for h3 in list_pareto:
            # print(f"Contient :\n{h3.event_tree}")
            frequence = calc_freq(h3.event_tree, frequence)

        general = before_construct_representative(h, frequence)

        print(f"La fréquence est {frequence}")
        print(f"La classe a {nbClass[-1]} histoires")
        print(f"L'histoire générale :\n{general.value}")
        print(f"time calcul : {time.time()-t}\n\n")

    print(f"""On a un total de {len(nbClass)} classes, pour {sum(nbClass)}, {nbClass} histoires et pour {sum(nbClassOrigin)} histoires à l'origine :
En moyenne {round(sum(nbClass)/len(nbClass),1)} histoires.
Les extrêmes : {min(nbClass)} et {max(nbClass)}.
Plus précisément : {"".join(f"\n- on a {y} classes qui ont {i} histoires" for i,y in dict(Counter(nb for nb in nbClass)).items())}""")

    t2 = time.time()
    print(f"time run all general : {t2-t1}")
    print(f"time global : {t2-t0}")


@register_method
def general_event_solution(setting, _, output):
    """Report a representative solution in each class."""
    import time

    t0 = time.time()
    result_number = synesth.solve(setting, event_vector_pareto @ history_generator)
    t1 = time.time()
    print(f"time run : {t1-t0}")
    nbClass = []
    nbClassOrigin = []
    for h in result_number:
        t = time.time()
        print(f"Classe :\n{h}")
        list_pareto = set()

        for h2 in result_number[h]:
            list_pareto.add(superdtlx.propagate_contents(History(setting.host_tree, h2.value).prune_unsampled()))
        nbClassOrigin.append(len(result_number[h]))
        nbClass.append(len(list_pareto))

        print(f"La classe a {nbClass[-1]} histoires")
        print(f"time calcul : {time.time()-t}\n\n")
    print(f"""On a un total de {len(nbClass)} classes, pour {sum(nbClass)} histoires et pour {sum(nbClassOrigin)} histoires à l'origine :
En moyenne {round(sum(nbClass)/len(nbClass),1)} histoires.
Les extrêmes : {min(nbClass)} et {max(nbClass)}.
Plus précisément : {"".join(f"\n- on a {y} classes qui ont {i} histoires" for i,y in dict(Counter(nb for nb in nbClass)).items())}""")
    t2 = time.time()
    print(f"time run all general : {t2-t1}")
    print(f"time global : {t2-t0}")


@register_method
def pareto_single(setting, _, output):
    """
    Report Pareto-optimal event count vectors along with the number of
    histories having that event count and an arbitrarily-selected history
    for each event count vector.
    """
    result = synesth.pareto_single(setting)

    for key in sorted(result.keys(), key=tuple):
        count, history = result[key]
        print(f"events={key}", file=output)
        print(f"count={count}", file=output)
        json.dump(history.to_mapping(), output)
        print("\n", file=output)


@register_method
def pareto_all(setting, _, output):
    """
    Report Pareto-optimal event count vectors along with all histories
    having that event count vector (may be exponentially slow!).
    """
    result = synesth.pareto_all(setting)

    for key in sorted(result.keys(), key=tuple):
        histories = result[key]
        print(f"events={key}", file=output)
        print(f"count={len(histories)}", file=output)

        for history in histories:
            json.dump(history.to_mapping(), output)
            print(file=output)

        print(file=output)


def reconcile(args):
    """Run the reconcile subcommand with the given arguments."""
    costs_dict: dict[str, float] = {}

    if args.cost:
        for cost_list in args.cost:
            for cost_entry in cost_list:
                kind, value = cost_entry.split("=", maxsplit=1)
                costs_dict[kind.replace("-", "_")] = literal_eval(value)

    costs = synesth.EventCosts(**costs_dict)
    setting = Reconciliation.from_mapping(json.load(args.input))

    if args.allow_unsampled:
        setting = Reconciliation(
            host_tree=graft_unsampled_hosts(setting.host_tree),
            associate_tree=setting.associate_tree,
        )

    setting.validate()
    methods[args.method](setting, costs, args.output)


def add_args(parser):
    """Add the reconcile subcommand to a command-line argument parser."""
    desc = __doc__

    subparser = parser.add_parser(
        "reconcile",
        description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    add_arg_input(subparser, "file defining the reconciliation problem")
    add_arg_output(subparser, "file where the results will be stored")

    subparser.add_argument(
        "--allow-unsampled",
        "-u",
        action="store_true",
        help="augment the host tree with candidate unsampled species",
    )

    event_types = [field.replace("_", "-") for field in synesth.EventCosts._fields]
    subparser.add_argument(
        "--cost",
        "-c",
        nargs="*",
        action="append",
        metavar="TYPE=VALUE",
        help=(
            "set the cost of an event type. repeat to set multiple costs. available "
            f"event types: {', '.join(event_types)} "
            "(default: use unit costs)"
        ),
    )

    methods_help = []

    for name, method in methods.items():
        raw_doc = inspect.getdoc(method)

        if raw_doc is None:
            doc = "no documentation"
        else:
            doc = inspect.getdoc(method).replace("\n", " ")[:-1].lower()

        methods_help.append(name + " (" + doc + ")")

    subparser.add_argument(
        "--method",
        "-m",
        metavar="METHOD",
        default="min-cost-single",
        choices=list(methods.keys()),
        help=(
            "select what to compute (default: %(default)s). "
            "available methods: " + ", ".join(methods_help)
        ),
    )

    subparser.add_argument(
        "--processes",
        "-p",
        metavar="NPROC",
        default=1,
        type=int,
        help="number of processes to spawn for computing (default: %(default)s)",
    )

    subparser.set_defaults(func=reconcile)
