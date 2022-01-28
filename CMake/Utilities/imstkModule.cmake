# Provides the list of all imstk modules
set(iMSTK_MODULE_LIST "")

# Add an option to toggle a module, adds to global list
macro(define_module moduleName moduleEnabled)
  option(${PROJECT_NAME}_MODULE_ENABLE_${moduleName} "Build with ${moduleName}" ${moduleEnabled})
  mark_as_superbuild(${PROJECT_NAME}_MODULE_ENABLE_${moduleName}:BOOL)
  mark_as_advanced(${PROJECT_NAME}_MODULE_ENABLE_${moduleName})

  if (${PROJECT_NAME}_MODULE_ENABLE_${moduleName})
    list(APPEND iMSTK_MODULE_LIST ${moduleName})
  endif()
endmacro()

# Update the module global list given the current state of all modules in the list
macro(update_module_list)
    set(iMSTK_MODULE_LIST_NEW "")
    foreach(module ${iMSTK_MODULE_LIST})
        if (${PROJECT_NAME}_MODULE_ENABLE_${module})
            list(APPEND iMSTK_MODULE_LIST_NEW ${module})
        endif()
    endforeach()
    set(iMSTK_MODULE_LIST ${iMSTK_MODULE_LIST_NEW})
    list(REMOVE_DUPLICATES iMSTK_MODULE_LIST_NEW)
endmacro()