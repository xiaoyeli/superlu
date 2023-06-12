# build HTML documentation using Doxygen

find_package(Doxygen)
set_package_properties("Doxygen" PROPERTIES
  DESCRIPTION "Documentation generator"
  URL "www.doxygen.org"
  PURPOSE "Generate HTML documentation from C sources")

# set Doxygen options
set(DOXYGEN_WARN_LOGFILE "doxygen.log")

# build documentation for target "all" if enable_doc is set
if(enable_doc)
  set(_DEPENDENCY_ALL "ALL")
endif()

if(DOXYGEN_FOUND)
  doxygen_add_docs(doc
                  "${CMAKE_CURRENT_SOURCE_DIR}/mainpage.txt"
                  "${CMAKE_CURRENT_SOURCE_DIR}/modules.txt"
                  "${PROJECT_SOURCE_DIR}/EXAMPLE"
                  "${PROJECT_SOURCE_DIR}/SRC"
                  "${PROJECT_SOURCE_DIR}/TESTING"
                  "${_DEPENDENCY_ALL}"
                  COMMENT "Generate HTML documentation with Doxygen")
elseif(enable_doc)
  message(WARNING "Disabling building documentation as Doxygen is missing.")
endif()
