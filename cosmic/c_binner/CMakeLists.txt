CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(c_py_bindings)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(PythonLibs_FIND_VERSION 2.7)
SET(PythonInterp_FIND_VERSION 2.7)

FIND_PACKAGE(SWIG)

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE(${SWIG_USE_FILE})

################
#HELPER FOR GETTING DETAILS OF PY FILE
FUNCTION(find_python_module module)
    STRING(TOUPPER ${module} module_upper)
    IF(NOT PY_${module_upper})
        IF(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
            SET(${module}_FIND_REQUIRED TRUE)
        ENDIF()
        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        EXECUTE_PROCESS(COMMAND 
            "${PYTHON_EXECUTABLE}" "-c" "import re, ${module}; print(re.compile('/${module}/__init__.py.*').sub('',${module}.__file__))"
            RESULT_VARIABLE _${module}_status 
            OUTPUT_VARIABLE _${module}_location
            ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        IF(NOT _${module}_status)
            SET(PY_${module_upper} ${_${module}_location} CACHE STRING 
                "Location of Python module ${module}")
        ENDIF(NOT _${module}_status)
    ENDIF(NOT PY_${module_upper})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PY_${module} DEFAULT_MSG PY_${module_upper})
ENDFUNCTION(find_python_module)
################

FIND_PYTHON_MODULE(numpy)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${PY_NUMPY}/numpy/core/include/)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

SET_SOURCE_FILES_PROPERTIES(c_binner.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(c_binner.i PROPERTIES SWIG_FLAGS "-threads" )
SWIG_ADD_MODULE(c_binner python c_binner.i)
SWIG_LINK_LIBRARIES(c_binner ) #add external libs here

