def getP(C, n):
    l = C.shape[0]
    P = np.zeros((l, l), dtype="c16")
    
    for a in range(l):
        for b in range(l):
            for i in range(n):
                P[b, a] += np.conj(C[a, i]) * C[b, i]
    return P

def getF(P, h, u, n):
    l = P.shape[0]
    F = np.zeros((l, l), dtype="c16")
    
    F += h
    for a in range(l):
        for b in range(l):
            for c in range(l):
                for d in range(l):
                    F[a, b] += P[d, c] * u[a, c, b, d]
    return F