CMAKE_MINIMUM_REQUIRED( VERSION 2.8.4 )

#Accomodates a project that's built in or out of source
#Use the src directory for out of source builds
#Create a test directory called UnitTests in the top level directory
#Put the unit test for each .cpp file in here.
#Make sure the CMakeModules directory exists and contains CodeCoverage.cmake
#Every *.cpp file must have a *Test.cpp file in the Test folder

#Targets

#make clean_cmake
#   deletes all the files created by cmake but leaves the executables

#make coverage
#   uses gcov, lcov, and genhtml to make html files for each .cpp file
#   The final output is index.html, which is found in the folder
#   Coverage/<fileName>. Example: Coverage/main/index.html

#make man
#   Creates the man.1 files

#make doc
#   Creates the doxygen documentation files

#make <filename>_coverage
#   makes the index.html for just <filename> and puts it in Coverage/<filename>

#make Experimental, Nightly, or Continuous
#   ...Start
#   ...Update
#   ...Configure
#   ...Build
#   ...Test
#   ...Coverage
#   ...MemCheck
#   ...Submit
# See bookmark CTest Targets for more details
#
#make clean_outputs
#   removes all the files in the UnitTests/outputs folder

PROJECT ( geom )

SET ( headerDir geom )

SET ( staticUsrLibPath /usr/local/lib )

SET ( projLibs 
${staticUsrLibPath}/libjson.a
${staticUsrLibPath}/libjson.a
)

SET ( ${CMAKE_PROJECT_NAME}_VERSION_MAJOR 1 )

SET ( ${CMAKE_PROJECT_NAME}_VERSION_MINOR 0 )

SET ( CMAKE_BUILD_TYPE Debug ) # Debug, Release, RelWithDebInfo

SET ( testDir ${CMAKE_SOURCE_DIR}/UnitTests ) 

SET ( outType static ) # exe, static, shared, header

SET ( installation "TRUE" )

#CMAKE_POLICY ( SET CMP0049 OLD ) #expansion of source file variables

################################################

INCLUDE ( CTest )

SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/"  )

IF ( EXISTS ${CMAKE_SOURCE_DIR}/src )
    SET ( srcDir ${CMAKE_SOURCE_DIR}/src )
ELSE ()
    SET ( srcDir ${CMAKE_SOURCE_DIR} )
ENDIF ()

FILE ( GLOB children ${CMAKE_INSTALL_PREFIX}/* )
FOREACH( child ${children} )
  IF ( IS_DIRECTORY ${child} )
    INCLUDE_DIRECTORIES ( ${child}/impl )
  ENDIF()
ENDFOREACH()
INCLUDE_DIRECTORIES ( ${srcDir} )

FILE ( GLOB_RECURSE projSrc ${srcDir}/*.cpp )

FILE ( MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/Coverage )

FILE ( MAKE_DIRECTORY ${testDir} )

FILE ( MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/CMakeModules )

ADD_CUSTOM_TARGET ( clean_cmake
    COMMAND ${CMAKE_COMMAND} -E remove *.cmake
    COMMAND ${CMAKE_COMMAND} -E remove DartConfiguration.tcl
    COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove_directory Testing
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove_directory Coverage
    COMMAND ${CMAKE_COMMAND} -E remove ${testDir}/outputs/*.*
    COMMAND ${CMAKE_COMMAND} -E remove Makefile
    COMMAND ${CMAKE_COMMAND} -E remove install_manifest.txt
)

ADD_CUSTOM_TARGET ( clean_outputs 
    COMMAND ${CMAKE_COMMAND} -E remove ${testDir}/outputs/*.*
)

ADD_CUSTOM_TARGET ( clean_test 
    COMMAND ${CMAKE_COMMAND} -E remove ${testDir}/outputs/*.*
    COMMAND ${CMAKE_COMMAND} -E remove ${testDir}/*.exe*
    COMMAND ${CMAKE_COMMAND} -E remove ${testDir}/*.stackdump*
)

ADD_CUSTOM_TARGET ( clean_all DEPENDS clean_cmake 
    IF ( EXISTS ${CMAKE_SOURCE_DIR}/doc/html )
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/doc/html
    ENDIF ()
)

ADD_CUSTOM_TARGET ( doc
    COMMAND doxygen ${CMAKE_SOURCE_DIR}/doc/Doxyfile
    #go to the doc/html directory and double click on the index.hhp file to build the .chm file with hhc.exe
)

IF ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
  
    MESSAGE ( "Debug" )

    IF ( CMAKE_COMPILER_IS_GNUCXX )
#SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -std=c++11 -std=gnu++11")
        SET ( CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")
        SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=gnu++11")
        SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
        IF ( CYGWIN )
            SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static")
        ENDIF ()
        INCLUDE ( CodeCoverage )
        FOREACH ( srcFile ${projSrc}  )
            GET_FILENAME_COMPONENT ( fname ${srcFile} NAME_WE )
            IF ( EXISTS ${testDir}/${fname}Test.cpp )
              IF ( EXISTS ${testDir}/${fname}Test.cmakeList  )
                  FILE ( READ ${testDir}/${fname}Test.cmakeList testSrc )
                  ADD_EXECUTABLE ( ${fname}_test ${testSrc} )
              ELSE ()
                  ADD_EXECUTABLE ( ${fname}_test ${testDir}/${fname}Test.cpp ${srcFile} )
              ENDIF ()
              SET_TARGET_PROPERTIES ( ${fname}_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${testDir} )
              TARGET_INCLUDE_DIRECTORIES ( ${fname}_test PUBLIC /usr/local/include )
              FOREACH( stLib ${projLibs} )
                  TARGET_LINK_LIBRARIES ( ${fname}_test ${stLib} )
              ENDFOREACH()
          #Setting target properties isn't recognized by add_test
              ADD_TEST ( ${fname} ${testDir}/${fname}_test)
              IF( NOT CYGWIN )
                SETUP_TARGET_FOR_COVERAGE ( ${fname}_coverage ${fname}_test Coverage/${fname} )
                IF ( covDep )
                    SET ( covDep ${covDep} ${fname}_coverage )
                ELSE ()
                    SET ( covDep ${fname}_coverage )
                ENDIF ()
              ENDIF()
            ENDIF ()
            

        ENDFOREACH ()

        IF( NOT CYGWIN )
          ADD_CUSTOM_TARGET ( coverage DEPENDS ${covDep} ) 

          SET ( CTEST_MEMORYCHECK_COMMAND "valgrind" )
          SET ( CTEST_MEMORYCHECK_COMMAND_OPTIONS "--log-file=valgr.log --leak-check=full" )
        ENDIF()

    ENDIF ()
  
    SET ( MYDEBUG ON )

    FOREACH ( ot ${outType}  )
      IF ( "${ot}" STREQUAL "exe" )
          ADD_EXECUTABLE ( ${CMAKE_PROJECT_NAME} ${projSrc} )
      ELSEIF ( "${ot}" STREQUAL "static" )
          ADD_LIBRARY ( ${CMAKE_PROJECT_NAME} STATIC ${projSrc} )
      ELSEIF ( "${ot}" STREQUAL "shared" )
          ADD_LIBRARY ( ${CMAKE_PROJECT_NAME} SHARED ${projSrc} )
      ELSEIF ( "${ot}" STREQUAL "header" )
          ADD_LIBRARY ( ${CMAKE_PROJECT_NAME} STATIC ${projSrc} )
      ENDIF ()
    ENDFOREACH()
    
    IF ( CMAKE_COMPILER_IS_GNUCXX AND NOT BUILD_SHARED_LIBS AND NOT CYGWIN )
        TARGET_LINK_LIBRARIES ( ${PROJECT_NAME} gcov )
    ENDIF ()
    
    TARGET_INCLUDE_DIRECTORIES (  ${PROJECT_NAME} PUBLIC /usr/local/include )
    FOREACH( stLib ${projLibs} )
        TARGET_LINK_LIBRARIES (  ${PROJECT_NAME} ${stLib} )
    ENDFOREACH()
    
    IF ( EXISTS ${CMAKE_SOURCE_DIR}/build-debug )
        SET_TARGET_PROPERTIES ( ${CMAKE_PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build-debug )
        SET_TARGET_PROPERTIES ( ${CMAKE_PROJECT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build-debug )
        SET_TARGET_PROPERTIES ( ${CMAKE_PROJECT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build-debug )
        SET ( installDir ${CMAKE_SOURCE_DIR}/build-debug )
    ELSE ()
        SET ( installDir ${CMAKE_SOURCE_DIR} )
    ENDIF ()
    IF ( EXISTS ${srcDir}/impl/${CMAKE_PROJECT_NAME}Config.h.in )
        CONFIGURE_FILE ( ${srcDir}/impl/${CMAKE_PROJECT_NAME}Config.h.in ${srcDir}/impl/${CMAKE_PROJECT_NAME}Config.h )  
    ENDIF ()
    

  
ELSE ( ) #Release code

    MESSAGE ( "Release" )
    
    SET ( MYDEBUG OFF )
    
    FOREACH ( ot ${outType}  )
      IF ( "${ot}" STREQUAL "exe" )
          ADD_EXECUTABLE ( ${CMAKE_PROJECT_NAME} ${projSrc} )
      ELSEIF ( "${ot}" STREQUAL "static" )
          ADD_LIBRARY ( ${CMAKE_PROJECT_NAME} STATIC ${projSrc} )
      ELSEIF ( "${ot}" STREQUAL "shared" )
          ADD_LIBRARY ( ${CMAKE_PROJECT_NAME} SHARED ${projSrc} )
      ELSEIF ( "${ot}" STREQUAL "header" )
          ADD_LIBRARY ( ${CMAKE_PROJECT_NAME} STATIC ${projSrc} )
      ENDIF ()
    ENDFOREACH()
    
    TARGET_INCLUDE_DIRECTORIES (  ${PROJECT_NAME} PUBLIC /usr/local/include )
    FOREACH( stLib ${projLibs} )
        TARGET_LINK_LIBRARIES (  ${PROJECT_NAME} ${stLib} )
    ENDFOREACH ()
    
    IF ( EXISTS ${CMAKE_SOURCE_DIR}/build-release )
        SET_TARGET_PROPERTIES ( ${CMAKE_PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build-release )
        SET_TARGET_PROPERTIES ( ${CMAKE_PROJECT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build-release )
        SET_TARGET_PROPERTIES ( ${CMAKE_PROJECT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build-release )
        SET ( installDir ${CMAKE_SOURCE_DIR}/build-release )
    ELSE ()
        SET ( installDir ${CMAKE_SOURCE_DIR} )
    ENDIF ()
    IF ( EXISTS ${srcDir}/impl/${CMAKE_PROJECT_NAME}Config.h.in )
        CONFIGURE_FILE ( ${srcDir}/impl/${CMAKE_PROJECT_NAME}Config.h.in ${srcDir}/impl/${CMAKE_PROJECT_NAME}Config.h )  
    ENDIF ()
    IF ( CMAKE_COMPILER_IS_GNUCXX )
        SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=gnu++11")
        IF ( CYGWIN )
            SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static")
        ENDIF ()
    ENDIF ()

ENDIF ()

###############################################################
#The lines below are dependent on whether you want to install to the system
#and create install targets, man pages, documentation

IF ( "${installation}" STREQUAL "TRUE" )
    
    ADD_CUSTOM_TARGET ( man 
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/man/*.gz
        COMMAND help2man ${installDir}/${CMAKE_PROJECT_NAME} -o ${CMAKE_SOURCE_DIR}/man/${CMAKE_PROJECT_NAME}.1
        COMMAND gzip ${CMAKE_SOURCE_DIR}/man/${CMAKE_PROJECT_NAME}.1
    )
    
    FILE ( GLOB projHeader ${srcDir}/*.h )
    FOREACH ( ot ${outType}  )
      IF ( "${ot}" STREQUAL "exe" )
              INSTALL ( TARGETS ${CMAKE_PROJECT_NAME} RUNTIME DESTINATION bin )
      ELSEIF ( "${ot}" STREQUAL "static" )
          INSTALL ( TARGETS ${CMAKE_PROJECT_NAME} ARCHIVE DESTINATION lib )
          INSTALL ( DIRECTORY ${srcDir}/ DESTINATION include/${headerDir} CONFIGURATIONS Release FILES_MATCHING PATTERN "*.h" )
      ELSEIF ( "${ot}" STREQUAL "shared" )
          INSTALL ( TARGETS ${CMAKE_PROJECT_NAME} LIBRARY DESTINATION lib )
          INSTALL ( DIRECTORY ${srcDir}/ DESTINATION include/${headerDir} CONFIGURATIONS Release FILES_MATCHING PATTERN "*.h" )
      ELSEIF ( "${ot}" STREQUAL "header" )
          INSTALL ( DIRECTORY ${srcDir}/ DESTINATION include/${headerDir} CONFIGURATIONS Release FILES_MATCHING PATTERN "*.h" )
      ENDIF ()
    ENDFOREACH()

    IF ( EXISTS ${CMAKE_SOURCE_DIR}/share/${CMAKE_PROJECT_NAME} )
        INSTALL ( DIRECTORY ${CMAKE_SOURCE_DIR}/share/${CMAKE_PROJECT_NAME} DESTINATION share CONFIGURATIONS Release COMPONENT support )
    ENDIF ()
    
    IF ( EXISTS ${CMAKE_SOURCE_DIR}/man/${CMAKE_PROJECT_NAME}.1.gz )
        INSTALL ( FILES ${CMAKE_SOURCE_DIR}/man/${CMAKE_PROJECT_NAME}.1.gz DESTINATION man/man1 CONFIGURATIONS Release COMPONENT support )
    ENDIF ()

ENDIF ()
