INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( sector.sif sector.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO sector DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_steady_periodic LABELS whitney mortar benchmark)
