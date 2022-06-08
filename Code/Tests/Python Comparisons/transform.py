def transform_spf(spf, C):
    return np.tensordot(C, spf, axes=((0), (0)))

def transform_one_body_elements(h, C, np, C_tilde=None):
    if C_tilde is None:
        C_tilde = C.conj().T

    return np.dot(C_tilde, np.dot(h, C))

def transform_two_body_elements(u, C, C_tilde=None):
    if C_tilde is None:
        C_tilde = C.conj().T

    # abcd, ds -> abcs
    _u = np.tensordot(u, C, axes=(3, 0))
    # abcs, cr -> absr -> abrs
    _u = np.tensordot(_u, C, axes=(2, 0)).transpose(0, 1, 3, 2)
    # abrs, qb -> arsq -> aqrs
    _u = np.tensordot(_u, C_tilde, axes=(1, 1)).transpose(0, 3, 1, 2)
    # pa, aqrs -> pqrs
    _u = np.tensordot(C_tilde, _u, axes=(1, 0))

    return _u