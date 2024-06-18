from .recurrence import reconcile
from .contents import propagate_contents
from ...model.history import History


def finalize_history(history: History):
    return propagate_contents(history.prune_unsampled())
