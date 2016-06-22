cmake_minimum_required(VERSION 2.8.9)

project(gimp-filmulator CXX C)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
set(CMAKE_POSITION_INDEPENDENT_CODE True)

find_package(PkgConfig)

find_package (Exiv2)
pkg_check_modules(GEXIV2 REQUIRED gexiv2)
MESSAGE( STATUS "GEXIV2_INCLUDE_DIRS:         " ${GEXIV2_INCLUDE_DIRS} )
MESSAGE( STATUS "GEXIV2_LIBRARY_DIRS:         " ${GEXIV2_LIBRARY_DIRS} )

MESSAGE(STATUS "Check for LCMS2 availability...")
#FIND_PACKAGE(LCMS2)
pkg_check_modules(LCMS2 REQUIRED lcms2) 

find_package (TIFF REQUIRED)

pkg_check_modules(PIXMAN pixman-1 REQUIRED) 
pkg_check_modules(PANGO pango) 
pkg_check_modules(PANGOFT2 pangoft2) 
pkg_check_modules(GLIB glib-2.0) 
pkg_check_modules(GMODULE gmodule-2.0) 
pkg_check_modules(GTHREAD gthread-2.0) 
pkg_check_modules(GOBJECT gobject-2.0) 
pkg_check_modules(GLIBMM glibmm-2.4) 
pkg_check_modules(GTKMM2 REQUIRED gtkmm-2.4) 

pkg_check_modules(GIMP_2_9 REQUIRED gimp-2.0>=2.9.0)
pkg_check_modules(GIMPUI REQUIRED gimpui-2.0>=2.9.0) 


SET(CMAKE_CXX_FLAGS_DEBUG " -Wall -O0 -g -DNDEBUG") 
SET(CMAKE_CXX_FLAGS_RELEASE " -O3 -ffast-math -ftree-vectorize -DNDEBUG -Wno-unused-parameter") 



include_directories(${EXIV2_INCLUDE_DIRS})
link_directories(${EXIV2_LIBRARY_DIRS})
include_directories(${GEXIV2_INCLUDE_DIRS})
link_directories(${GEXIV2_LIBRARY_DIRS})

include_directories(${TIFF_INCLUDE_DIRS})
link_directories(${TIFF_LIBRARY_DIRS})

include_directories(${LCMS22_INCLUDE_DIRS})
link_directories(${LCMS22_LIBRARY_DIRS})
include_directories(${GLIB_INCLUDE_DIRS})
link_directories(${GLIB_LIBRARY_DIRS})
include_directories(${GMODULE_INCLUDE_DIRS})
link_directories(${GMODULE_LIBRARY_DIRS}  )
include_directories(${GOBJECT_INCLUDE_DIRS})
link_directories(${GOBJECT_LIBRARY_DIRS}  )
include_directories(${GTHREAD_INCLUDE_DIRS})
link_directories(${GTHREAD_LIBRARY_DIRS}  )

include_directories(${GLIBMM_INCLUDE_DIRS})
link_directories(${GLIBMM_LIBRARY_DIRS})

include_directories(${PIXMAN_INCLUDE_DIRS})
link_directories(${PIXMAN_LIBRARY_DIRS})

link_directories(${GTKMM2_LIBRARY_DIRS}  )
include_directories(${GTKMM2_INCLUDE_DIRS}  )
set(GTKMM_LIBRARIES ${GTKMM2_LIBRARIES})


include_directories(${GIMP_2_9_INCLUDE_DIRS})
link_directories(${GIMP_2_9_LIBRARY_DIRS})
include_directories(${GIMPUI_INCLUDE_DIRS})
link_directories(${GIMPUI_LIBRARY_DIRS})

add_executable(gimp-filmulator gimp-filmulator.cc)

set_target_properties(gimp-filmulator PROPERTIES COMPILE_FLAGS " -DHAVE_GIMP_2_9=1 ")

target_link_libraries(gimp-filmulator 
    ${TIFF_LIBRARIES} ${LCMS2_LIBRARIES} 
    ${EXIV2_LIBRARIES}
    ${GEXIV2_LIBRARIES}
    ${PIXMAN_LIBRARIES} 
    ${GIMP_LIBS}
    ${GIMP_2_9_LIBRARIES}
    ${GIMPUI_LIBRARIES}
    ${GTKMM_LIBRARIES} 
    ${GLIBMM_LIBRARIES} 
    ${GLIB_LIBRARIES} 
    ${GMODULE_LIBRARIES} 
    ${GOBJECT_LIBRARIES} 
    ${GTHREAD_LIBRARIES} 
)
  
  
  INSTALL(TARGETS gimp-filmulator RUNTIME DESTINATION bin)
  