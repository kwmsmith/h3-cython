from libc cimport stdlib
from libc cimport stdint

DEF MAX_RES=15
DEF DEGS_TO_RADS= 0.017453292519943295
DEF RADS_TO_DEGS= 57.29577951308232
DEF MAX_CELL_BNDRY_VERTS = 10

cdef extern from "h3api.h":

    ctypedef stdint.uint64_t H3Index

    ctypedef struct GeoCoord:
        double lat  # latitude in radians
        double lon  # longitude in radians

    ctypedef struct GeoBoundary:
        int numVerts
        GeoCoord verts[MAX_CELL_BNDRY_VERTS]

    H3Index geoToH3(const GeoCoord *, int)
    void h3ToGeo(H3Index, GeoCoord *)
    void h3ToGeoBoundary(H3Index, GeoBoundary *)
    double edgeLengthKm(int)
    double degsToRads(double)
    double radsToDegs(double)
    int h3GetBaseCell(H3Index)
    int h3IsValid(H3Index)
    int h3IsResClassIII(H3Index)
    int h3IsPentagon(H3Index)
    void kRing(H3Index, int, H3Index*)
    int maxKringSize(int)
    void kRingDistances(H3Index, int, H3Index*, int*)
    int hexRange(H3Index, int, H3Index*)
    int hexRangeDistances(H3Index, int, H3Index*, int*)


cdef int _check_res(int res) except -1:
    if not 0 <= res <= MAX_RES:
        raise ValueError(
            f"resolution {res} not within bounds (0, 15)."
        )
    return res


cpdef inline double degs_to_rads(double degrees):
    return DEGS_TO_RADS * degrees


cpdef inline double rads_to_degs(double radians):
    return RADS_TO_DEGS * radians


cdef inline double _mercator_lat(double lat):
    """Helper coerce lat range"""
    return lat - 180 if lat > 90 else lat


cdef inline double _mercator_lon(double lon):
    """Helper coerce lng range"""
    return lon - 360 if lon > 180 else lon


def edge_length_km(int res):
    return edgeLengthKm(_check_res(res))


def h3_to_string(H3Index h3_int):
    return hex(h3_int)[2:]


def string_to_h3(str h3_address):
    return int(h3_address, 16)


def geo_to_h3(double lat, double lon, int res):
    """Index a geo-coordinate at a resolution into an h3 address"""
    cdef GeoCoord geo_coord = GeoCoord(
        lat=degs_to_rads(_mercator_lat(lat)),
        lon=degs_to_rads(_mercator_lon(lon))
    )
    return h3_to_string(geoToH3(&geo_coord, _check_res(res)))


def h3_to_geo(str h3_address):
    """Reverse lookup an h3 address into a geo-coordinate"""
    cdef GeoCoord geo_coord
    h3ToGeo(string_to_h3(h3_address), &geo_coord)
    return (
        _mercator_lat(rads_to_degs(geo_coord.lat)),
        _mercator_lon(rads_to_degs(geo_coord.lon))
    )


def h3_to_geo_boundary(str h3_address, geo_json=False):
    """Compose an array of geo-coordinates that outlines a hexagonal cell"""
    cdef GeoBoundary geo_boundary
    cdef int i
    h3ToGeoBoundary(string_to_h3(h3_address), &geo_boundary)
    out = []
    for i in range(geo_boundary.numVerts):
        out.append([
            _mercator_lon(rads_to_degs(geo_boundary.verts[i].lon)),
            _mercator_lat(rads_to_degs(geo_boundary.verts[i].lat))
        ]) if geo_json else out.append([
            _mercator_lat(rads_to_degs(geo_boundary.verts[i].lat)),
            _mercator_lon(rads_to_degs(geo_boundary.verts[i].lon))
        ])
    if geo_json:
        out.append(out[0])
    return out


cpdef int h3_get_resolution(str h3_address):
    """Returns the resolution of an `h3_address`

    :return: nibble (0-15)
    """
    return int(h3_address[1], 16)


cpdef int h3_get_base_cell(h3_address):
    return h3GetBaseCell(string_to_h3(h3_address))


cpdef bint h3_is_valid(h3_address):
    """Validates an `h3_address`

    :returns: boolean
    """
    try:
        return h3IsValid(string_to_h3(h3_address))
    except Exception:
        return False


cpdef bint h3_is_res_class_III(h3_address):
    return h3IsResClassIII(string_to_h3(h3_address))

cpdef bint h3_is_pentagon(h3_address):
    return h3IsPentagon(string_to_h3(h3_address))


cdef _hexagon_c_array_to_set(H3Index *h3_addresses, int array_len):
    cdef set s = set()
    cdef int i
    cdef H3Index h3_address
    for i in range(array_len):
        h3_address = h3_addresses[i]
        if h3_address != 0:
            s.add(h3_to_string(h3_address))
    return s
        

def k_ring(h3_address, int ring_size):
    """Get K-Rings for a given hexagon"""
    if ring_size < 0:
        raise ValueError(f"ring_size must be >= 0, given {ring_size}")
    cdef int array_len = maxKringSize(ring_size)
    cdef H3Index *kring_array = <H3Index *> stdlib.calloc(array_len, sizeof(H3Index))
    if not kring_array:
        raise MemoryError()
    try:
        kRing(string_to_h3(h3_address), ring_size, kring_array)
        return _hexagon_c_array_to_set(kring_array, array_len)
    finally:
        if kring_array:
            stdlib.free(kring_array)


def k_ring_distances(h3_address, int ring_size):
    """Get K-Rings for a given hexagon properly split by ring"""
    cdef:
        int array_len, i
        int *distance_array
        H3Index *kring_array
        list out
    if ring_size < 0:
        raise ValueError(f"ring_size must be >= 0, given {ring_size}.")
    try:
        array_len = maxKringSize(ring_size)
        kring_array = <H3Index *>stdlib.calloc(array_len, sizeof(H3Index))
        if not kring_array:
            raise MemoryError()
        distance_array = <int *>stdlib.calloc(array_len, sizeof(int))
        if not distance_array:
            raise MemoryError()
        kRingDistances(
            string_to_h3(h3_address),
            ring_size,
            kring_array,
            distance_array,
        )
        out = []
        for i in range(0, ring_size + 1):
            out.append(set())
        for i in range(0, array_len):
            ring_index = distance_array[i]
            out[ring_index].add(h3_to_string(kring_array[i]))
        return out
    finally:
        if kring_array:
            stdlib.free(kring_array)
        if distance_array:
            stdlib.free(distance_array)


def hex_range(h3_address, int ring_size):
    cdef:
        int array_len
        H3Index *kring_array
    try:
        array_len = maxKringSize(ring_size)
        kring_array = <H3Index *>stdlib.calloc(array_len, sizeof(H3Index))
        if not kring_array:
            raise MemoryError()
        success = hexRange(string_to_h3(h3_address), ring_size, kring_array)
        if success != 0:
            raise ValueError('Specified hexagon range contains a pentagon')
        return _hexagon_c_array_to_set(kring_array, array_len)
    finally:
        if kring_array:
            stdlib.free(kring_array)


# def hex_range_distances(h3_address, int ring_size):
#     """
#     Get K-Rings for a given hexagon properly split by ring,
#     aborting if a pentagon is reached
#     """
#     cdef:
#         int i, array_len, success
#         int *distance_array
#         H3Index *kring_array
#         list out
#     if ring_size < 0:
#         raise ValueError(f"ring_size must be >= 0, given {ring_size}.")
#     try:
#         array_len = maxKringSize(ring_size)
#         kring_array = <H3Index *>stdlib.calloc(array_len, sizeof(H3Index))
#         if not kring_array:
#             raise MemoryError()
#         distance_array = <int *>stdlib.calloc(array_len, sizeof(int))
#         if not distance_array:
#             raise MemoryError()
#         success = hexRangeDistances(
#             string_to_h3(h3_address),
#             ring_size,
#             kring_array,
#             distance_array,
#         )
#         if success != 0:
#             raise ValueError('Specified hexagon range contains a pentagon')
#         out = []
#         for i in range(0, ring_size + 1):
#             out.append(set())
#         for i in range(0, array_len):
#             ring_index = distance_array[i]
#             out[ring_index].add(h3_to_string(kring_array[i]))
#         return out
#     finally:
#         if kring_array:
#             stdlib.free(kring_array)
#         if distance_array:
#             stdlib.free(distance_array)



# def polyfill(geo_json, res, geo_json_conformant=False):
#     """
#     Get hexagons for a given GeoJSON region

#     :param geo_json dict: A GeoJSON dictionary
#     :param res int: The hexagon resolution to use (0-15)
#     :param geo_json_conformant bool: Determines (lat, lng) vs (lng, lat)
#         ordering Default is false, which is (lat, lng) ordering, violating
#         the spec http://geojson.org/geojson-spec.html#id2 which is (lng, lat)

#     :returns: Set of hex addresses
#     """
#     geo_json_lite = _geo_json_to_geo_json_lite(geo_json, geo_json_conformant)
#     array_len = libh3.maxPolyfillSize(byref(geo_json_lite), res)
#     HexagonArray = c_long * array_len
#     hexagons = HexagonArray()
#     libh3.polyfill(byref(geo_json_lite), res, hexagons)
#     return hexagon_c_array_to_set(hexagons)

