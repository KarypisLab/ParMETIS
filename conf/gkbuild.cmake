# Helper modules.
include(CheckFunctionExists)
include(CheckIncludeFile)

# Setup options.
option(OPENMP "enable OpenMP support" OFF)
option(PCRE "enable PCRE support" OFF)
option(GKREGEX "enable GKREGEX support" OFF)
option(GKRAND "enable GKRAND support" OFF)


# Add compiler flags.
if(MSVC)
  set(GKlib_COPTIONS WIN32 MSC _CRT_SECURE_NO_DEPRECATE USE_GKREGEX)
elseif(WIN32)
  set(GKlib_COPTIONS USE_GKREGEX)
else()
  set(GKlib_COPTIONS LINUX FILE_OFFSET_BITS=64)
endif(MSVC)

if(CYGWIN)
  list(APPEND GKlib_COPTIONS CYGWIN)
endif()

if(APPLE)
  list(APPEND GKlib_COPTIONS MACOS)
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  list(APPEND GKlib_COPTS -fno-strict-aliasing -Werror -Wall -pedantic -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unknown-pragmas -Wno-unused-label)
endif()

if(UNIX)
include(CheckPIESupported)
check_pie_supported()
set(CMAKE_POSITION_INDEPENDENT_CODE true)
endif()

# Find OpenMP if it is requested.
if(OPENMP)
  find_package(OpenMP REQUIRED)
  list(APPEND GKlib_COPTIONS __OPENMP__)
endif()

# Set the CPU type
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  list(APPEND GKlib_COPTIONS NO_X86=1)
endif()

# Add various options
if(PCRE)
  list(APPEND GKlib_COPTIONS __WITHPCRE__)
endif()

if(GKREGEX)
list(APPEND GKlib_COPTIONS USE_GKREGEX)
endif()

if(GKRAND)
list(APPEND GKlib_COPTIONS USE_GKRAND)
endif()


# Check for features.
check_include_file(execinfo.h HAVE_EXECINFO_H)
if(HAVE_EXECINFO_H)
  list(APPEND GKlib_COPTIONS HAVE_EXECINFO_H)
endif(HAVE_EXECINFO_H)

check_function_exists(getline HAVE_GETLINE)
if(HAVE_GETLINE)
 list(APPEND GKlib_COPTIONS HAVE_GETLINE)
endif(HAVE_GETLINE)


# Custom check for TLS.
if(MSVC)
 list(APPEND GKlib_COPTIONS __thread=__declspec(thread))

  # This if checks if that value is cached or not.
  if("${HAVE_THREADLOCALSTORAGE}" MATCHES "^${HAVE_THREADLOCALSTORAGE}$")
    message(CHECK_START "checking for thread-local storage")
    try_compile(HAVE_THREADLOCALSTORAGE
      ${CMAKE_BINARY_DIR}
      ${CMAKE_CURRENT_LIST_DIR}/check_thread_storage.c)
    if(HAVE_THREADLOCALSTORAGE)
      message(CHECK_PASS "found")
    else()
      message(CHECK_FAIL "not found")
    endif()
  endif()
  if(NOT HAVE_THREADLOCALSTORAGE)
    list(APPEND GKlib_COPTIONS __thread=)
  endif()
endif()

# Finally set the official C flags.
add_compile_options("$<$<COMPILE_LANGUAGE:C>:${GKlib_COPTS}>")
add_compile_definitions("$<$<COMPILE_LANGUAGE:C>:${GKlib_COPTIONS}>")
