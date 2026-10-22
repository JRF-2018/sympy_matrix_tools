# -*- coding: utf-8 -*-
__version__ = '0.7.0' # Time-stamp: <2026-10-22T17:31:49Z>

"""Execute the code cells of the notebooks in examples/."""

import json
import os
import pytest

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def run_notebook (name):
    """Execute the code cells of examples/name.ipynb and return the
    namespace.  Shell commands (lines starting with '!') are skipped."""
    with open(os.path.join(EXAMPLES, name + ".ipynb"), encoding="utf-8") as f:
        nb = json.load(f)
    ns = {"__name__": "__main__"}
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"])
        src = "\n".join(l for l in src.split("\n")
                        if not l.lstrip().startswith("!"))
        exec(compile(src, name + ".ipynb", "exec"), ns)
    return ns


@pytest.mark.parametrize("name", ["lq_control", "chain_rule",
                                  "first_order_logic"])
def test_notebook (name):
    run_notebook(name)


def test_lq_control_ledger ():
    d = run_notebook("lq_control")["d"]
    assert d.check(trials=3)
    s = d.summary()
    assert s["symbolic_failures"] == 0
    assert s["numeric_failures"] == 0
    assert s["open_obligations"] == 0
    assert s["asserted"] == 0
    assert "x = -K**(-1)*B.T*(W*(A*y0 + z) + p)" in d.report()
    assert d._repr_latex_().startswith(r"$\begin{array}")
