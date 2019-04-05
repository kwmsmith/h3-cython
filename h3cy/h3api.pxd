from libc cimport stdlib
from libc cimport stdint

include "_h3.pxi"

cdef extern from "h3api.h":

    ctypedef stdint.uint64_t H3Index

    ctypedef struct GeoCoord:
        double lat  # latitude in radians
        double lon  # longitude in radians

    ctypedef struct GeoBoundary:
        int numVerts
        GeoCoord verts[MAX_CELL_BNDRY_VERTS]

    ctypedef struct Geofence:
        int numVerts
        GeoCoord *verts

    ctypedef struct GeoPolygon:
        Geofence geofence
        int numHoles
        Geofence *holes

    ctypedef struct LinkedGeoCoord:
        GeoCoord vertex
        LinkedGeoCoord *next

    ctypedef struct LinkedGeoLoop:
        LinkedGeoCoord *first
        LinkedGeoCoord *last
        LinkedGeoLoop *next

    ctypedef struct LinkedGeoPolygon:
        LinkedGeoLoop *first
        LinkedGeoLoop *last
        LinkedGeoPolygon *next


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
    int hexRanges(H3Index*, int, int, H3Index*)
    int hexRing(H3Index, int, H3Index*)
    int h3Line(H3Index, H3Index, H3Index*)
    int h3LineSize(H3Index, H3Index)
    int h3Distance(H3Index, H3Index)
    H3Index h3ToParent(H3Index, int)
    void h3ToChildren(H3Index, int, H3Index *)
    int maxH3ToChildrenSize(H3Index, int)
    int compact(const H3Index *, H3Index *, const int)
    int uncompact(const H3Index *, const int, H3Index *, const int, const int)
    int maxUncompactSize(const H3Index *, const int, const int)
    void polyfill(const GeoPolygon*, int, H3Index*)
    int maxPolyfillSize(const GeoPolygon*, int)
    double hexAreaKm2(int)
    double hexAreaM2(int)
    double edgeLengthKm(int)
    double edgeLengthM(int)
    stdint.int64_t numHexagons(int)
    int res0IndexCount()
    void getRes0Indexes(H3Index *)
    void h3SetToLinkedGeo(const H3Index*, const int, LinkedGeoPolygon*)
    void destroyLinkedPolygon(LinkedGeoPolygon*)
    int h3IndexesAreNeighbors(H3Index, H3Index)
    H3Index getH3UnidirectionalEdge(H3Index, H3Index)
    int h3UnidirectionalEdgeIsValid(H3Index)
    H3Index getOriginH3IndexFromUnidirectionalEdge(H3Index)
    H3Index getDestinationH3IndexFromUnidirectionalEdge(H3Index)
    void getH3IndexesFromUnidirectionalEdge(H3Index, H3Index*)
    void getH3UnidirectionalEdgesFromHexagon(H3Index, H3Index*)
    void getH3UnidirectionalEdgeBoundary(H3Index, GeoBoundary*)
    int h3IsResClassIII(H3Index)
