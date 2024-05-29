"""Compute a minimum-cost (super-)reconciliation of two trees."""

import argparse
import json
import inspect
from ast import literal_eval
from tqdm import tqdm
from functools import partial
from pathos.multiprocessing import Pool
from .util import add_arg_input, add_arg_output
from ..model.history import Reconciliation, History, graft_unsampled_hosts
from ..utils.algebras import Structure, MinPlus
from ..compute.util import (
    DummyProgress,
    DummyPool,
    EventCosts,
    event_vector_pareto,
    history_counter,
    history_projector,
    history_generator,
    partial_history_projector,
    partial_history_generator,
    HistoryBuilder,
)
from ..compute import superdtlx
from sowing.node import Node
from collections import Counter

methods = {}


def register_method(method):
    pretty_name = method.__name__.replace("_", "-")
    methods[pretty_name] = method


@register_method
def single_solution(run, costs, output):
    """
    Report a single arbitrary minimum-cost solution
    along with the total number of optimal solutions.
    """
    min_cost = Structure(MinPlus, costs.event_cost_morphism)
    result = run(structure=min_cost * (history_counter + history_projector))
    cost, (count, history) = result

    print(f"cost={cost}", file=output)
    print(f"count={count}", file=output)
    yield history.value


@register_method
def all_solutions(run, costs, output):
    """Report all minimum-cost solutions."""
    min_cost = Structure(MinPlus, costs.event_cost_morphism)
    result = run(structure=min_cost * history_generator)
    cost, histories = result

    print(f"cost={cost}", file=output)
    print(f"count={len(histories)}", file=output)

    for solution in histories:
        yield solution.value

@register_method
def general_solution(run, a, output):
    """Report a representative solution in each class."""

    list_result = []
    while True :
        a,b,c = general_solution_copy(run, a, output)
        print(len(list_result), c)
        for i in list_result:
            if i[0] != a:
                print(i)
                print("\n\n\n")
                print(a,b)
                print(len(list_result))
                return
                yield
            if i[1] !=b:
                print("numA",i[1],"-",b)
            else:
                print("numB",i[1],"-",b)
        list_result.append([a,b])

    return
    yield

def general_solution_copy(run, _, output):
    """Report a representative solution in each class."""
    def get_sub_history(hist) :
        data, edge = get_data_edges(hist)
        if data.apparent:
            return HistoryBuilder(Node("continue")), hist
        else:
            h = HistoryBuilder(Node(data))
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
            newNode = HistoryBuilder(Node(node))
            for i in hist:
                newNode = newNode * i
            return newNode
        else:
            newNode = HistoryBuilder(Node(data))
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
    result_number = run(structure=partial_history_generator @ (event_vector_pareto @ history_generator))
    t1 = time.time()
    #print(f"time run : {t1-t0}")
    nbClass = []
    nbClassOrigin = []
    result_number_new = {}

    for clas in result_number:
        t = superdtlx.propagate_contents(History(run.keywords['setting'].host_tree, clas.value).prune_unsampled())
        if t not in result_number_new :
            result_number_new[t] = {}
        for vector in result_number[clas]:
            if vector not in result_number_new[t]:
                result_number_new[t][vector] = set()
            result_number_new[t][vector] = result_number_new[t][vector] | result_number[clas][vector]


    print("aa",len(result_number), len(result_number_new))

    for clas in result_number_new:
        t = time.time()
        #print(f"Classe :\n{h.value}")
        frequence = {}
        list_pareto = set()

        s = 0
        for vector in result_number_new[clas]:
            s += len(result_number_new[clas][vector])
            for hist in result_number_new[clas][vector]:
                #breakpoint()
                #list_pareto.add(History(run.keywords['setting'].host_tree, hist.value))
                list_pareto.add(superdtlx.propagate_contents(History(run.keywords['setting'].host_tree, hist.value).prune_unsampled()))

        nbClassOrigin.append(s)
        nbClass.append(len(list_pareto))

        #for h3 in list_pareto:
            # print(f"Contient :\n{h3.event_tree}")
        #    frequence = calc_freq(h3.event_tree, frequence)

        #general = before_construct_representative(h, frequence)

        #print(f"La fréquence est {frequence}")
        #print(f"La classe a {nbClass[-1]} histoires")
        #print(f"L'histoire générale :\n{general.value}")
        #print(f"time calcul : {time.time()-t}\n\n")

    #print(f"""On a un total de {len(nbClass)} classes, pour {sum(nbClass)}, {sorted(nbClass)} histoires et pour {sum(nbClassOrigin)} histoires à l'origine :
#En moyenne {round(sum(nbClass)/len(nbClass),1)} histoires.
#Les extrêmes : {min(nbClass)} et {max(nbClass)}.
#Plus précisément : {"".join(f"\n- on a {y} classes qui ont {i} histoires" for i,y in dict(Counter(nb for nb in nbClass)).items())}""")

    #t2 = time.time()
    #print(f"time run all general : {t2-t1}")
    #print(f"time global : {t2-t0}")

    #import hashlib
    #print(hashlib.sha256(str(result_number).encode()).hexdigest())

    return result_number, sum(nbClass), sum(nbClassOrigin)
    #yield

@register_method
def general_event_solution(run, _, output):
    """Report a representative solution in each class."""
    import time

    t0 = time.time()
    result_number = run(structure=event_vector_pareto @ history_generator)
    t1 = time.time()
    print(f"time run : {t1-t0}")
    nbClass = []
    nbClassOrigin = []
    for h in result_number:
        t = time.time()
        print(f"Classe :\n{h}")
        list_pareto = set()

        for h2 in result_number[h]:
            list_pareto.add(superdtlx.propagate_contents(History(run.keywords['setting'].host_tree, h2.value).prune_unsampled()))
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
    return
    yield


@register_method
def pareto(run, _, output):
    """
    Compute all Pareto-optimal event count vectors and
    the number of corresponding solutions for each vector.
    """
    result = run(structure=event_vector_pareto @ history_counter)

    for key in sorted(result.keys(), key=tuple):
        print(f"{key}: {result[key]}")

    return
    yield


def reconcile(args):
    """Run the reconcile subcommand with the given arguments."""
    costs_dict: dict[str, float] = {}

    if args.cost:
        for cost_list in args.cost:
            for cost_entry in cost_list:
                kind, value = cost_entry.split("=", maxsplit=1)
                costs_dict[kind.replace("-", "_")] = literal_eval(value)

    costs = EventCosts(**costs_dict)
    setting = Reconciliation.from_mapping(json.load(args.input))

    if args.allow_unsampled:
        setting = Reconciliation(
            host_tree=graft_unsampled_hosts(setting.host_tree),
            associate_tree=setting.associate_tree,
        )

    setting.validate()

    for event_tree in methods[args.method](
        partial(
            superdtlx.reconcile,
            setting=setting,
            progress=DummyProgress,
            pool=DummyPool(),
            # FIXME: Restore multiprocess support
            # progress=tqdm,
            # pool=Pool(args.processes),
        ),
        costs,
        args.output,
    ):
        history = History(setting.host_tree, event_tree)
        history = superdtlx.finalize_history(history)
        history.validate()
        json.dump(history.to_mapping(), args.output)
        print(file=args.output)


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

    event_types = [field.replace("_", "-") for field in EventCosts._fields]
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
        doc = inspect.getdoc(method).replace("\n", " ")[:-1].lower()
        methods_help.append(name + " (" + doc + ")")

    subparser.add_argument(
        "--method",
        "-m",
        metavar="METHOD",
        default="single-solution",
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
