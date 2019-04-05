from libc cimport stdlib
from libc cimport stdint
cimport h3api

include "_h3.pxi"

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


def h3_to_string(h3api.H3Index h3_int):
    return hex(h3_int)[2:]


def string_to_h3(str h3_address):
    return int(h3_address, 16)


def geo_to_h3(double lat, double lon, int res):
    """Index a geo-coordinate at a resolution into an h3 address"""
    cdef h3api.GeoCoord geo_coord = h3api.GeoCoord(
        lat=degs_to_rads(_mercator_lat(lat)),
        lon=degs_to_rads(_mercator_lon(lon))
    )
    return h3_to_string(h3api.geoToH3(&geo_coord, _check_res(res)))


def h3_to_geo(str h3_address):
    """Reverse lookup an h3 address into a geo-coordinate"""
    cdef h3api.GeoCoord geo_coord
    h3api.h3ToGeo(string_to_h3(h3_address), &geo_coord)
    return (
        _mercator_lat(rads_to_degs(geo_coord.lat)),
        _mercator_lon(rads_to_degs(geo_coord.lon))
    )


def h3_to_geo_boundary(str h3_address, geo_json=False):
    """Compose an array of geo-coordinates that outlines a hexagonal cell"""
    cdef h3api.GeoBoundary geo_boundary
    cdef int i
    h3api.h3ToGeoBoundary(string_to_h3(h3_address), &geo_boundary)
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
    return h3api.h3GetBaseCell(string_to_h3(h3_address))


cpdef bint h3_is_valid(h3_address):
    """Validates an `h3_address`

    :returns: boolean
    """
    try:
        return h3api.h3IsValid(string_to_h3(h3_address))
    except Exception:
        return False


cpdef bint h3_is_res_class_III(h3_address):
    return h3api.h3IsResClassIII(string_to_h3(h3_address))

cpdef bint h3_is_pentagon(h3_address):
    return h3api.h3IsPentagon(string_to_h3(h3_address))


cdef _hexagon_c_array_to_list(h3api.H3Index *h3_addresses, int array_len):
    cdef:
        out = []
        int i
        h3api.H3Index h3_address
    for i in range(array_len):
        h3_address = h3_addresses[i]
        if h3_address != 0:
            out.append(h3_to_string(h3_address))
    return out


cdef _hexagon_c_array_to_set(h3api.H3Index *h3_addresses, int array_len):
    cdef set s = set()
    cdef int i
    cdef h3api.H3Index h3_address
    for i in range(array_len):
        h3_address = h3_addresses[i]
        if h3_address != 0:
            s.add(h3_to_string(h3_address))
    return s


def k_ring(h3_address, int ring_size):
    """Get K-Rings for a given hexagon"""
    cdef:
        int array_len
        h3api.H3Index *kring_array = NULL

    if ring_size < 0:
        raise ValueError(f"ring_size must be >= 0, given {ring_size}")

    try:
        array_len = h3api.maxKringSize(ring_size)
        kring_array = <h3api.H3Index *> stdlib.calloc(
            array_len,
            sizeof(h3api.H3Index)
        )
        if not kring_array:
            raise MemoryError()
        h3api.kRing(string_to_h3(h3_address), ring_size, kring_array)
        return _hexagon_c_array_to_set(kring_array, array_len)
    finally:
        if kring_array:
            stdlib.free(kring_array)


def k_ring_distances(h3_address, int ring_size):
    """Get K-Rings for a given hexagon properly split by ring"""
    cdef:
        int array_len, i
        int *distance_array = NULL
        h3api.H3Index *kring_array = NULL
        list out
    if ring_size < 0:
        raise ValueError(f"ring_size must be >= 0, given {ring_size}.")
    try:
        array_len = h3api.maxKringSize(ring_size)
        kring_array = <h3api.H3Index *>stdlib.calloc(
            array_len,
            sizeof(h3api.H3Index)
        )
        if not kring_array:
            raise MemoryError()
        distance_array = <int *>stdlib.calloc(array_len, sizeof(int))
        if not distance_array:
            raise MemoryError()
        h3api.kRingDistances(
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
        h3api.H3Index *kring_array = NULL
    try:
        array_len = h3api.maxKringSize(ring_size)
        kring_array = <h3api.H3Index *>stdlib.calloc(
            array_len,
            sizeof(h3api.H3Index)
        )
        if not kring_array:
            raise MemoryError()
        success = h3api.hexRange(
            string_to_h3(h3_address),
            ring_size,
            kring_array
        )
        if success != 0:
            raise ValueError('Specified hexagon range contains a pentagon')
        return _hexagon_c_array_to_set(kring_array, array_len)
    finally:
        if kring_array:
            stdlib.free(kring_array)


def hex_range_distances(h3_address, int ring_size):
    """
    Get K-Rings for a given hexagon properly split by ring,
    aborting if a pentagon is reached
    """
    cdef:
        int i, array_len, success
        int *distance_array = NULL
        h3api.H3Index *kring_array = NULL
        list out
    if ring_size < 0:
        raise ValueError(f"ring_size must be >= 0, given {ring_size}.")
    try:
        array_len = h3api.maxKringSize(ring_size)
        kring_array = <h3api.H3Index *>stdlib.calloc(
            array_len,
            sizeof(h3api.H3Index)
        )
        if not kring_array:
            raise MemoryError()
        distance_array = <int *>stdlib.calloc(array_len, sizeof(int))
        if not distance_array:
            raise MemoryError()
        success = h3api.hexRangeDistances(
            string_to_h3(h3_address),
            ring_size,
            kring_array,
            distance_array,
        )
        if success != 0:
            raise ValueError('Specified hexagon range contains a pentagon')
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


def hex_ranges(h3_address_list, int ring_size):
    """
    Get K-Rings for all hexagons properly split by ring,
    aborting if a pentagon is reached
    """
    cdef:
        int num_hexagons = len(h3_address_list)
        int array_len, i, ring_index, ring_end, range_size, success
        h3api.H3Index *hex_array = NULL
        h3api.H3Index *kring_array = NULL
        dict out
        list hex_range_list

    if ring_size < 0:
        raise ValueError(f"ring_size must be >= 0, given {ring_size}.")

    try:

        array_len = num_hexagons * h3api.maxKringSize(ring_size)

        kring_array = <h3api.H3Index *>stdlib.calloc(
            array_len,
            sizeof(h3api.H3Index)
        )
        if not kring_array:
            raise MemoryError()

        hex_array = <h3api.H3Index *>stdlib.calloc(
            num_hexagons,
            sizeof(h3api.H3Index)
        )
        if not hex_array:
            raise MemoryError()

        for i in range(num_hexagons):
            hex_array[i] = string_to_h3(h3_address_list[i])

        success = h3api.hexRanges(
            hex_array,
            num_hexagons,
            ring_size,
            kring_array,
        )

        if success != 0:
            raise ValueError(
                'One of the specified hexagon ranges contains a pentagon')

        out = {}
        for i in range(0, num_hexagons):

            h3_address = h3_address_list[i]
            hex_range_list = []
            out[h3_address] = hex_range_list

            for j in range(0, ring_size + 1):
                hex_range_list.append(set())

            ring_index = 0
            ring_end = 0
            range_size = int(array_len / num_hexagons)

            for j in range(0, range_size):
                if j > ring_end:
                    ring_index = ring_index + 1
                    ring_end = ring_end + 6 * ring_index
                # hexRanges doesn't return distance array
                hex_range_list[ring_index].add(
                    h3_to_string(kring_array[i * range_size + j]))

        return out

    finally:
        if kring_array:
            stdlib.free(kring_array)
        if hex_array:
            stdlib.free(hex_array)


def hex_ring(h3_address, int ring_size):
    """
    Get a hexagon ring for a given hexagon.
    Returns individual rings, unlike `k_ring`.

    If a pentagon is reachable, falls back to a
    MUCH slower form based on `k_ring`.
    """

    # This technically should be defined in the C code,
    # but this is much faster
    cdef:
        int array_len = 6 * ring_size
        int success
        h3api.H3Index *hex_ring_array = NULL

    try:
        hex_ring_array = <h3api.H3Index *>stdlib.calloc(
            array_len,
            sizeof(h3api.H3Index)
        )
        if not hex_ring_array:
            raise MemoryError()
        success = h3api.hexRing(
            string_to_h3(h3_address), ring_size, hex_ring_array
        )
        if success != 0:
            raise Exception(
                f'Failed to get hexagon ring for pentagon {h3_address}'
            )

        return _hexagon_c_array_to_set(hex_ring_array, array_len)

    finally:
        if hex_ring_array:
            stdlib.free(hex_ring_array)


def h3_line(start, end):
    cdef:
        int line_size
        h3api.H3Index *line = NULL

    try:
        line_size = h3api.h3LineSize(string_to_h3(start), string_to_h3(end))
        line = <h3api.H3Index *>stdlib.calloc(line_size, sizeof(h3api.H3Index))
        if not line:
            raise MemoryError()
        success = h3api.h3Line(string_to_h3(start), string_to_h3(end), line)
        if success != 0:
            raise Exception(
                f"Unable to find line between {start} and {end}"
            )
        return _hexagon_c_array_to_list(line, line_size)

    finally:
        if line:
            stdlib.free(line)


def h3_distance(h3_address_origin, h3_address_h3):
    return h3api.h3Distance(
        string_to_h3(h3_address_origin),
        string_to_h3(h3_address_h3)
    )


def h3_to_parent(h3_address, int res):
    return h3_to_string(
        h3api.h3ToParent(
            string_to_h3(h3_address),
            res
        )
    )


def h3_to_children(h3_address, int res):
    cdef:
        h3api.H3Index h3_index = string_to_h3(h3_address)
        int max_children

    res = _check_res(res)
    try:
        max_children = h3api.maxH3ToChildrenSize(h3_index, res)
        children_array = <h3api.H3Index *>stdlib.calloc(
            max_children, sizeof(h3api.H3Index)
        )
        if not children_array:
            raise MemoryError()
        h3api.h3ToChildren(h3_index, res, children_array)
        return _hexagon_c_array_to_set(children_array, max_children)
    finally:
        if children_array:
            stdlib.free(children_array)


def compact(h3_addresses):
    cdef:
        int num_hexagons, i
        h3api.H3Index *hex_set_array = NULL
        h3api.H3Index *compacted_hex_set_array = NULL

    if not h3_addresses:
        return set()

    try:
        num_hexagons = len(h3_addresses)
        hex_set_array = <h3api.H3Index *>stdlib.calloc(
            num_hexagons,
            sizeof(h3api.H3Index)
        )
        if not hex_set_array:
            raise MemoryError()

        compacted_hex_set_array = <h3api.H3Index *>stdlib.calloc(
            num_hexagons,
            sizeof(h3api.H3Index)
        )
        if not compacted_hex_set_array:
            raise MemoryError()

        for i, address in enumerate(h3_addresses):
            hex_set_array[i] = int(address, 16)

        ret_val = h3api.compact(
            hex_set_array,
            compacted_hex_set_array,
            num_hexagons
        )

        if ret_val != 0:
            raise Exception(
                'Failed to compact, malformed input data (duplicate hexagons?)'
            )

        return _hexagon_c_array_to_set(compacted_hex_set_array, num_hexagons)

    finally:
        if hex_set_array:
            stdlib.free(hex_set_array)
        if compacted_hex_set_array:
            stdlib.free(compacted_hex_set_array)


def uncompact(h3_addresses, int res):
    cdef:
        int i, num_hexagons, max_uncompacted_num
        h3api.H3Index *hex_set_array = NULL
        h3api.H3Index *hex_uncompacted_set_array = NULL

    res = _check_res(res)

    if not h3_addresses:
        return set()

    try:
        num_hexagons = len(h3_addresses)
        hex_set_array = <h3api.H3Index *>stdlib.calloc(
            num_hexagons,
            sizeof(h3api.H3Index)
        )
        if not hex_set_array:
            raise MemoryError()

        for i, address in enumerate(h3_addresses):
            hex_set_array[i] = int(address, 16)

        max_uncompacted_num = h3api.maxUncompactSize(
            hex_set_array,
            num_hexagons,
            res
        )
        if max_uncompacted_num < 0:
            raise Exception(
                'Failed to determine max uncompact '
                'output size (bad resolution?)'
            )

        hex_uncompacted_set_array = <h3api.H3Index *>stdlib.calloc(
            max_uncompacted_num,
            sizeof(h3api.H3Index)
        )
        if not hex_uncompacted_set_array:
            raise MemoryError()

        ret_val = h3api.uncompact(
            hex_set_array,
            num_hexagons,
            hex_uncompacted_set_array,
            max_uncompacted_num,
            res
        )
        if ret_val != 0:
            raise Exception('Failed to uncompact (bad resolution?)')

        return _hexagon_c_array_to_set(
            hex_uncompacted_set_array,
            max_uncompacted_num
        )

    finally:
        if hex_set_array:
            stdlib.free(hex_set_array)
        if hex_uncompacted_set_array:
            stdlib.free(hex_uncompacted_set_array)


cdef h3api.GeoCoord _coord_array_to_geo_coord(
    list coord_array,
    bint geo_json_conformant=False
) except *:
    cdef:
        double lat, lon

    if len(coord_array) != 2:
        raise ValueError("invalid coord array, expecting length 2.")

    if geo_json_conformant:
        lon = coord_array[0]
        lat = coord_array[1]
    else:
        lat = coord_array[0]
        lon = coord_array[1]

    return h3api.GeoCoord(
        degs_to_rads(_mercator_lat(lat)),
        degs_to_rads(_mercator_lon(lon))
    )


cdef h3api.Geofence _polygon_array_to_geofence(
    list polygon_array,
    bint geo_json_conformant=False
) except *:

    cdef:
        int i, num_verts
        h3api.GeoCoord *geo_coord_array = NULL

    try:
        num_verts = len(polygon_array)
        geo_coord_array = <h3api.GeoCoord *>stdlib.calloc(
            num_verts,
            sizeof(h3api.GeoCoord)
        )
        if not geo_coord_array:
            raise MemoryError()

        for i in range(num_verts):
            geo_coord_array[i] = _coord_array_to_geo_coord(
                polygon_array[i],
                geo_json_conformant
            )

    except: # noqa
        if geo_coord_array:
            stdlib.free(geo_coord_array)
        raise

    else:
        return h3api.Geofence(num_verts, geo_coord_array)


cdef void _free_geofence(h3api.Geofence geofence):
    if geofence.verts:
        stdlib.free(geofence.verts)


cdef h3api.GeoPolygon _geo_json_to_geo_polygon(
    geo_json,
    geo_json_conformant=False
) except *:

    cdef:
        h3api.Geofence geofence
        int num_holes, i
        h3api.Geofence *holes = NULL

    if geo_json['type'] != 'Polygon':
        raise TypeError('Only Polygon GeoJSON supported')

    try:
        num_holes = len(geo_json['coordinates']) - 1
        geofence = _polygon_array_to_geofence(
            geo_json['coordinates'][0],
            geo_json_conformant
        )

        if num_holes > 0:
            holes = <h3api.Geofence *>stdlib.calloc(
                num_holes, sizeof(h3api.Geofence)
            )
            if not holes:
                raise MemoryError()
            for i in range(num_holes):
                holes[i] = _polygon_array_to_geofence(
                    geo_json['coordinates'][i + 1],
                    geo_json_conformant
                )

    except: # noqa
        _free_geofence(geofence)
        if holes:
            for i in range(num_holes):
                _free_geofence(holes[i])
            stdlib.free(holes)
        raise

    else:
        return h3api.GeoPolygon(geofence, num_holes, holes)


def polyfill(geo_json, int res, bint geo_json_conformant=False):
    """
    Get hexagons for a given GeoJSON region

    :param geo_json dict: A GeoJSON dictionary
    :param res int: The hexagon resolution to use (0-15)
    :param geo_json_conformant bool: Determines (lat, lng) vs (lng, lat)
        ordering Default is false, which is (lat, lng) ordering, violating
        the spec http://geojson.org/geojson-spec.html#id2 which is (lng, lat)

    :returns: Set of hex addresses
    """

    cdef:
        h3api.GeoPolygon geo_polygon = h3api.GeoPolygon(
            h3api.Geofence(0, NULL), 0, NULL
        )
        int array_len
        h3api.H3Index *hex_array = NULL

    res = _check_res(res)

    try:
        geo_polygon = _geo_json_to_geo_polygon(geo_json, geo_json_conformant)
        array_len = h3api.maxPolyfillSize(&geo_polygon, res)
        hex_array = <h3api.H3Index *>stdlib.calloc(
            array_len, sizeof(h3api.H3Index)
        )
        if not hex_array:
            raise MemoryError()
        h3api.polyfill(&geo_polygon, res, hex_array)
        return _hexagon_c_array_to_set(hex_array, array_len)

    finally:
        if hex_array:
            stdlib.free(hex_array)

        _free_geofence(geo_polygon.geofence)

        if geo_polygon.holes:
            for i in range(geo_polygon.numHoles):
                _free_geofence(geo_polygon.holes[i])
            stdlib.free(geo_polygon.holes)


cpdef double hex_area_km2(int res) except -1:
    return h3api.hexAreaKm2(_check_res(res))


cpdef double hex_area_m2(int res) except -1:
    return h3api.hexAreaM2(_check_res(res))


cpdef double edge_length_km(int res) except -1:
    return h3api.edgeLengthKm(_check_res(res))


cpdef double edge_length_m(int res) except -1:
    return h3api.edgeLengthM(_check_res(res))


cpdef stdint.int64_t num_hexagons(int res) except -1:
    return h3api.numHexagons(_check_res(res))


def get_res_zero_indexes():
    cdef:
        int idx_count
        h3api.H3Index *h3_addresses = NULL

    try:
        idx_count = h3api.res0IndexCount()
        h3_addresses = <h3api.H3Index *>stdlib.calloc(
            idx_count, sizeof(h3api.H3Index)
        )
        if not h3_addresses:
            raise MemoryError()
        h3api.getRes0Indexes(h3_addresses)
        return _hexagon_c_array_to_set(
            h3_addresses, idx_count
        )

    finally:
        if h3_addresses:
            stdlib.free(h3_addresses)


def h3_set_to_multi_polygon(h3_addresses, geo_json=False):
    """
    Get the outlines of a set of H3 hexagons, returned in GeoJSON MultiPolygon
    format (an array of polygons, each with an array of loops, each an array of
    coordinates). Coordinates are returned as [lat, lng] pairs unless GeoJSON
    is requested.

    :param h3_addresses string[]: H3 addresses to get outlines for
    :param geo_json bool: Whether to follow GeoJSON: [lng, lat], closed loops
    :returns: MultiPolygon-style output.
    """

    cdef:
        int address_count
        h3api.H3Index *hexagon_array = NULL
        h3api.LinkedGeoPolygon polygon = h3api.LinkedGeoPolygon(
            NULL, NULL, NULL
        )
        h3api.LinkedGeoPolygon *polygon_ptr = NULL
        h3api.LinkedGeoPolygon *original_polygon = NULL
        h3api.LinkedGeoLoop *loop_ptr = NULL
        list output
        int lat_index, lng_index

    # Early exit on empty input
    if not h3_addresses:
        return []

    try:
        # Set up input set
        address_count = len(h3_addresses)
        hexagon_array = <h3api.H3Index *>stdlib.calloc(
            address_count, sizeof(h3api.H3Index)
        )
        if not hexagon_array:
            raise MemoryError()
        for i, address in enumerate(h3_addresses):
            hexagon_array[i] = string_to_h3(address)

        # Store a reference to the first polygon - that's the one we need for
        # memory deallocation
        original_polygon = &polygon
        polygon_ptr = &polygon
        h3api.h3SetToLinkedGeo(hexagon_array, address_count, &polygon)

        # Loop through the linked structure, building the output
        output = []
        lat_index = 1 if geo_json else 0
        lng_index = 0 if geo_json else 1

        while polygon_ptr:
            loops = []
            output.append(loops)
            # Follow ->first pointer
            loop_ptr = polygon_ptr.first
            while loop_ptr:

                coords = []
                loops.append(coords)
                # Follow ->first pointer
                coord_ptr = loop_ptr.first

                while coord_ptr:

                    pair = [None, None]
                    coords.append(pair)
                    pair[lat_index] = _mercator_lat(
                        rads_to_degs(coord_ptr.vertex.lat)
                    )
                    pair[lng_index] = _mercator_lon(
                        rads_to_degs(coord_ptr.vertex.lon)
                    )
                    # Follow ->next pointer
                    coord_ptr = coord_ptr.next

                if geo_json:
                    # Close loop if GeoJSON is requested
                    coords.append(coords[0])

                # Follow ->next pointer
                loop_ptr = loop_ptr.next

            # Follow ->next pointer
            polygon_ptr = polygon_ptr.next

        # Clean up
        h3api.destroyLinkedPolygon(original_polygon)
        return output

    finally:
        if hexagon_array:
            stdlib.free(hexagon_array)
