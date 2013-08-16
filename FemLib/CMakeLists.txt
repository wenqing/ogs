# Source files
GET_SOURCE_FILES(SOURCES_FEMLIB)
SET ( SOURCES ${SOURCES_FEMLIB})
GET_SOURCE_FILES(SOURCES_Element FiniteElement)
SET ( SOURCES ${SOURCES} ${SOURCES_Element})
GET_SOURCE_FILES(SOURCES_Integration Integration)
SET ( SOURCES ${SOURCES} ${SOURCES_Integration})
GET_SOURCE_FILES(SOURCES_ShapeFunction ShapeFunction)
SET ( SOURCES ${SOURCES} ${SOURCES_ShapeFunction})
GET_SOURCE_FILES(SOURCES_Extrapolation Extrapolation)
SET ( SOURCES ${SOURCES} ${SOURCES_Extrapolation})
GET_SOURCE_FILES(SOURCES_CoordinatesMapping CoordinatesMapping)
SET ( SOURCES ${SOURCES} ${SOURCES_CoordinatesMapping})

# Create the library
ADD_LIBRARY(FemLib STATIC ${SOURCES})

include_directories(
	.
	${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/BaseLib
    ${CMAKE_SOURCE_DIR}/GeoLib
    ${CMAKE_SOURCE_DIR}/MathLib
    ${CMAKE_SOURCE_DIR}/MeshLib
)

INCLUDE_DIRECTORIES (${EIGEN_INCLUDE_DIR})

target_link_libraries (FemLib
	BaseLib
	GeoLib
	MathLib
	MeshLib
)

