# Source files
get_source_files(SOURCES)
get_source_files(SOURCES_ADSORPTION Adsorption)
set(SOURCES ${SOURCES} ${SOURCES_ADSORPTION})

get_source_files(SOURCES_FLUID_DENSITY Fluid/Density)
set(SOURCES ${SOURCES} ${SOURCES_FLUID_DENSITY})

get_source_files(SOURCES_PERMEABILITY Fluid/Permeability)
set(SOURCES ${SOURCES} ${SOURCES_PERMEABILITY})

get_source_files(SOURCES_STORAGE Fluid/Storage)
set(SOURCES ${SOURCES} ${SOURCES_STORAGE})

get_source_files(SOURCES_VISICOSITY Fluid/Viscosity)
set(SOURCES ${SOURCES} ${SOURCES_VISICOSITY})

get_source_files(SOURCES_SOLID_DENSITY Solid/Density)
set(SOURCES ${SOURCES} ${SOURCES_SOLID_DENSITY})

add_library(MaterialLib ${SOURCES} )
target_link_libraries(MaterialLib
    BaseLib
)
