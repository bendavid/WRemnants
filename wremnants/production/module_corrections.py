"""Python booking helpers for the C++ action helpers in module_corrections.hpp.

Each booking function constructs the corresponding C++ helper, books it on
the supplied RDataFrame, and patches the resulting RResultPtr so that
``GetValue()`` / ``GetPtr()`` lazily materializes a python-friendly object
(numpy arrays for the dense helpers, a scipy CSR array for the sparse one).
"""

import numpy as np
import ROOT

import narf
import narf.clingutils

narf.clingutils.Declare('#include "module_corrections.hpp"')


def _wrap_lazy_result(res, builder):
    """Patch ``res`` so its dereference methods return ``builder(cpp_result)``.

    The conversion is performed exactly once, on first access. The wrapped
    object still triggers the underlying RDataFrame action via ``GetPtr``.
    """
    res._GetPtr = res.GetPtr
    res._py_result = None

    def get():
        if res._py_result is None:
            res._py_result = builder(res._GetPtr())
        return res._py_result

    ret_null = lambda: None
    res.__deref__ = get
    res.__follow__ = get
    res.GetPtr = get
    res.GetValue = get
    res.begin = ret_null
    res.end = ret_null
    return res


def book_grad_helper(df, nparms, cols):
    """Book a wrem::GradHelper on ``df`` and return a result whose value is
    a numpy array of length ``nparms``.

    Parameters
    ----------
    df : ROOT.RDataFrame
    nparms : int
        Number of parameters (length of the gradient).
    cols : sequence of str
        Either ``[grad_vec, idx_vec]`` or ``[grad_vec, idx_vec, weight]``.
    """
    helper = ROOT.wrem.GradHelper(int(nparms))
    res = df.Book(helper, list(cols))

    def to_numpy(cpp_vec):
        return np.array(cpp_vec, dtype=np.float64, copy=True)

    return _wrap_lazy_result(res, to_numpy)


def book_hess_helper(df, nparms, cols):
    """Book a wrem::HessHelper on ``df`` and return a result whose value is a
    full symmetric ``(nparms, nparms)`` numpy array.

    The underlying ``SymMatrixAtomic`` only stores the upper triangle; this
    helper materializes the dense symmetric form.
    """
    helper = ROOT.wrem.HessHelper(int(nparms))
    res = df.Book(helper, list(cols))
    n = int(nparms)

    def to_numpy(cpp_mat):
        out = np.zeros((n, n), dtype=np.float64)
        row = np.zeros(n, dtype=np.float64)
        for i in range(n):
            cpp_mat.fill_row(i, row)
            out[i] = row
        # fill_row only populates row[i:]; mirror to the lower triangle.
        i_idx, j_idx = np.tril_indices(n, k=-1)
        out[i_idx, j_idx] = out[j_idx, i_idx]
        return out

    return _wrap_lazy_result(res, to_numpy)


def book_hess_helper_sparse(df, nparms, cols):
    """Book a wrem::HessHelperSparse on ``df`` and return a result whose value
    is a scipy CSR array of shape ``(nparms, nparms)``.

    The underlying ``SparseMatrixAtomic`` only stores the upper triangle of a
    symmetric matrix; off-diagonal entries are mirrored to the lower triangle
    in the returned CSR.
    """
    import scipy.sparse

    helper = ROOT.wrem.HessHelperSparse(int(nparms))
    res = df.Book(helper, list(cols))
    n = int(nparms)

    def to_scipy(cpp_mat):
        cpp_idxvals = cpp_mat.index_values()
        i = np.array(cpp_idxvals.idxs0(), dtype=np.int64, copy=True)
        j = np.array(cpp_idxvals.idxs1(), dtype=np.int64, copy=True)
        v = np.array(cpp_idxvals.vals(), dtype=np.float64, copy=True)
        # Mirror off-diagonal entries (only upper triangle is stored).
        diag = i == j
        i_full = np.concatenate([i, j[~diag]])
        j_full = np.concatenate([j, i[~diag]])
        v_full = np.concatenate([v, v[~diag]])
        return scipy.sparse.coo_array(
            (v_full, (i_full, j_full)), shape=(n, n)
        ).tocsr()

    return _wrap_lazy_result(res, to_scipy)
