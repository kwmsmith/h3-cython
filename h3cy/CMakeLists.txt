# TODO: Version stuff -- ensure version matches H3's version.

set(LIB_SOURCE_FILES
    h3lib/include/bbox.h
    h3lib/include/polygon.h
    h3lib/include/polygonAlgos.h
    h3lib/include/h3Index.h
    h3lib/include/h3UniEdge.h
    h3lib/include/geoCoord.h
    h3lib/include/vec2d.h
    h3lib/include/vec3d.h
    h3lib/include/linkedGeo.h
    h3lib/include/localij.h
    h3lib/include/baseCells.h
    h3lib/include/faceijk.h
    h3lib/include/vertexGraph.h
    h3lib/include/mathExtensions.h
    h3lib/include/constants.h
    h3lib/include/coordijk.h
    h3lib/include/algos.h
    h3lib/include/stackAlloc.h
    h3lib/lib/algos.c
    h3lib/lib/coordijk.c
    h3lib/lib/bbox.c
    h3lib/lib/polygon.c
    h3lib/lib/h3Index.c
    h3lib/lib/vec2d.c
    h3lib/lib/vec3d.c
    h3lib/lib/linkedGeo.c
    h3lib/lib/localij.c
    h3lib/lib/geoCoord.c
    h3lib/lib/h3UniEdge.c
    h3lib/lib/mathExtensions.c
    h3lib/lib/vertexGraph.c
    h3lib/lib/faceijk.c
    h3lib/lib/baseCells.c)

set(UNCONFIGURED_API_HEADER h3lib/include/h3api.h.in)
set(CONFIGURED_API_HEADER h3lib/include/h3api.h)
configure_file(${UNCONFIGURED_API_HEADER} ${CONFIGURED_API_HEADER})

add_cython_target(_h3)
# add_library(_h3 MODULE ${_h3})
add_library(_h3 MODULE ${_h3} ${LIB_SOURCE_FILES} ${CONFIGURED_API_HEADER})
python_extension_module(_h3)

set_property(TARGET _h3 PROPERTY C_STANDARD 99)

find_library(M_LIB m)
if(M_LIB)
    # target_link_libraries(_h3 PUBLIC ${M_LIB})
    target_link_libraries(_h3 ${M_LIB})
endif()
if(have_alloca)
    target_compile_definitions(_h3 PUBLIC H3_HAVE_ALLOCA)
endif()
if(have_vla)
    target_compile_definitions(_h3 PUBLIC H3_HAVE_VLA)
endif()
target_include_directories(_h3 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/h3lib/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/h3lib/include>)

install(TARGETS _h3 LIBRARY DESTINATION h3cy)
