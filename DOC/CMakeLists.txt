# build HTML documentation using Doxygen

find_package(Doxygen)
set_package_properties("Doxygen" PROPERTIES
  DESCRIPTION "Documentation generator"
  URL "www.doxygen.org"
  PURPOSE "Generate HTML documentation from C sources")

# set Doxygen options
set(DOXYGEN_WARN_LOGFILE "doxygen.log")

doxygen_add_docs(doc
                 "${CMAKE_CURRENT_SOURCE_DIR}/mainpage.txt"
                 "${PROJECT_SOURCE_DIR}/EXAMPLE"
                 "${PROJECT_SOURCE_DIR}/SRC"
                 "${PROJECT_SOURCE_DIR}/TESTING"
                 COMMENT "Generate HTML documentation with Doxygen")
