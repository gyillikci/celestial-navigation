''' Minimal numerical-derivative helper for validating analytic Jacobians.

    The GTSAM 4.2.1 pip wheel does not expose gtsam.numericalDerivative*, so this
    is a trimmed port of the develop-branch pure-Python implementation
    (gtsam/python/gtsam/utils/numerical_derivative.py, Joel Truher & Frank
    Dellaert) — enough to check a scalar-or-vector function of a single manifold
    argument against central finite differences.

    © GTSAM 2010-2019 (BSD); trimmed for this project.  MIT-compatible use.
'''

import numpy as np


def _local(a, b):
    if isinstance(a, np.ndarray):
        return b - a
    if isinstance(a, (float, int, np.floating, np.integer)):
        return np.array([float(b) - float(a)])
    return a.localCoordinates(b)


def _retract(a, xi):
    if isinstance(a, (np.ndarray, float, int)):
        return a + xi
    return a.retract(xi)


def numericalDerivative11(h, x, delta=1e-5):
    ''' Central-difference Jacobian of y=h(x) w.r.t. the tangent of x. '''
    hx = h(x)
    m = _local(hx, hx).shape[0]
    n = _local(x, x).shape[0]
    dx = np.zeros(n)
    jac = np.zeros((m, n))
    f = 1.0 / (2.0 * delta)
    for j in range(n):
        dx[j] = delta
        dy1 = _local(hx, h(_retract(x, dx)))
        dx[j] = -delta
        dy2 = _local(hx, h(_retract(x, dx)))
        dx[j] = 0.0
        jac[:, j] = (dy1 - dy2) * f
    return jac
