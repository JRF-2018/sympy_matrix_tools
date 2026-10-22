# Development

[Back to README](../README.md)

## Tests

You need SymPy >= 1.14 and pytest.

```sh
python -m pytest
```

`pytest.ini` runs the tests in `tests`, the doctests in the modules of
`sympy_matrix_tools`, and the doctests in `README.md` and `docs/*.md`.
The notebooks in `examples` are executed by `tests/test_examples.py`.

The tests of the chain rule (`tests/test_matrix_function.py`)
substitute concrete functions for `MatrixFunction`'s and compare the
results with the direct differentiation by `derive_by_array`.  The
tests of the first-order logic (`tests/test_fol.py`) check also that
non-theorems are rejected.

## Scheduled tests on GitHub

`.github/workflows/sympy-latest.yml` runs the tests with the latest
SymPy every Monday (and by hand from the "Actions" tab).  If they
fail, it opens an issue labeled `sympy-latest`, unless such an issue
is already open; when they pass again, it closes the issue with a
comment.  GitHub disables scheduled workflows of a repository without
activity for 60 days; enable it again in the "Actions" tab.

To try the path that opens an issue, run the workflow by hand with
"force_failure" checked: the tests fail on purpose and an issue is
opened (it says that the failure was forced).  A following normal run
closes the issue.


## Conventions

  * Each source file has `__version__` and a `Time-stamp` near the top
    (after the module docstring).  Update both when the file is
    changed.
  * Each module has a docstring, and the main functions have
    docstrings with doctests, so that `help()` is useful to humans and
    AI agents.
  * Write a space between a function name and its parameters in
    definitions: `def name (args):`.
  * Source files use CRLF line endings.
  * Everything public (code, comments, documents and examples) is
    written in English.
  * Record changes in [CHANGELOG.md](../CHANGELOG.md).
