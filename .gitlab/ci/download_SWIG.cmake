cmake_minimum_required(VERSION 3.12)

# Download link only on windows, currently we only support windows for SWIG
# if (NOT "$ENV{CMAKE_CONFIGURATION}" MATCHES "windows")
#   message(FATAL_ERROR
#     "Unknown platform for SWIG")
# endif ()

set(swig_version "4.0.2")
set(swig_url_root "http://prdownloads.sourceforge.net/swig/")
set(subdir_name "swigwin-${swig_version}")
set(filename "swigwin-${swig_version}.zip")
set(sha256sum "daadb32f19fe818cb9b0015243233fc81584844c11a48436385e87c050346559") # \todo: Fix

# Download the file.
file(DOWNLOAD
  "${swig_url_root}/${filename}"
  ".gitlab/${filename}"
  STATUS download_status
  EXPECTED_HASH "SHA256=${sha256sum}")

# Check the download status.
list(GET download_status 0 res)
if (res)
  list(GET download_status 1 err)
  message(FATAL_ERROR
    "Failed to download ${filename}: ${err}")
endif ()

# Extract the file.
execute_process(
  COMMAND
    "${CMAKE_COMMAND}"
    -E tar
    xf "${filename}"
  WORKING_DIRECTORY ".gitlab"
  RESULT_VARIABLE res
  ERROR_VARIABLE err
  ERROR_STRIP_TRAILING_WHITESPACE)
if (res)
  message(FATAL_ERROR
    "Failed to extract ${filename}: ${err}")
endif ()

# Move to a predictable prefix.
file(RENAME
  ".gitlab/${subdir_name}"
  ".gitlab/swig")
