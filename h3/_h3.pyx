cdef extern from "h3api.h":
    double edgeLengthKm(int res)

DEF MAX_RES=15

def edge_length_km(int res):
    if not 0 <= res <= MAX_RES:
        raise ValueError(
            f"resolution {res} not within bounds (0, 15)."
        )
    return edgeLengthKm(res)
