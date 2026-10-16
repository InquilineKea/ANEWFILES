"""Starting point: the default parameters of the Twisted Kagome Moire Bands page.

solution() must return a dict of model parameters (see task.md for names and bounds).
kagome_core.py in this folder is an editable copy of the solver; the evaluator uses its own
sealed copy and re-solves whatever parameters you return.
"""


def solution():
    return {"m": 3, "t": 1.0, "breath": 0.0, "flux": 0.0, "tp": 0.35, "xi": 0.25, "rc": 0.8,
            "V": 0.0, "sx": 0.0, "sy": 0.0}
