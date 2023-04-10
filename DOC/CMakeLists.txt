# build HTML documentation using Doxygen

find_package(Doxygen)
set_package_properties("Doxygen" PROPERTIES
  DESCRIPTION "Documentation generator"
  URL "www.doxygen.org"
  PURPOSE "Generate HTML documentation from C sources")

doxygen_add_docs(doc
                 "${PROJECT_SOURCE_DIR}/SRC"
                 COMMENT "Generate HTML documentation with Doxygen")
