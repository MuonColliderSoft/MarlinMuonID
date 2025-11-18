########################################################
# cmake file for building a Marlin Package
# @author Karol Krizka
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( MarlinMuonID )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )

cmake_policy(SET CMP0008 NEW)

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

SET( COMPILER_FLAGS "-Wno-effc++" )
MESSAGE( STATUS "FLAGS ${COMPILER_FLAGS}" )
FOREACH( FLAG ${COMPILER_FLAGS} )
  CHECK_CXX_COMPILER_FLAG( "${FLAG}" CXX_FLAG_WORKS_${FLAG} )
  IF( ${CXX_FLAG_WORKS_${FLAG}} )
    MESSAGE ( STATUS "Adding ${FLAG} to CXX_FLAGS" )
    ### We append the flag, so they overwrite things from somewhere else
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG} ")
  ELSE()
    MESSAGE ( STATUS "NOT Adding ${FLAG} to CXX_FLAGS" )
  ENDIF()
ENDFOREACH()

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDRec)
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )
INCLUDE( DD4hep )

FIND_PACKAGE( ROOT REQUIRED) ## COMPONENTS Hist)
SET( ROOT_COMPONENT_LIBRARIES Geom Reflex)

if(DD4HEP_USE_XERCESC)
  find_package(XercesC)
endif()
include(DD4hep_XML_setup)
INCLUDE_DIRECTORIES( BEFORE SYSTEM ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES}  )

INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

FIND_PACKAGE( RAIDA REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM ${RAIDA_INCLUDE_DIRS} )
LINK_LIBRARIES( ${RAIDA_LIBRARIES} )
ADD_DEFINITIONS( ${RAIDA_DEFINITIONS} )

FIND_PACKAGE( GSL REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM ${GSL_INCLUDE_DIRS} )

FIND_PACKAGE(PandoraSDK 03.00.00 REQUIRED)
INCLUDE_DIRECTORIES( SYSTEM ${PandoraSDK_INCLUDE_DIRS} )

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### DATA ###########################################################

### LIBRARY #################################################################

# Main library with Processors
ADD_LIBRARY(MarlinMuonID SHARED
  ./src/MarlinMuonID.cc
  )
TARGET_INCLUDE_DIRECTORIES(MarlinMuonID
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  PRIVATE
  ${Marlin_INCLUDE_DIRS}
  ${RAIDA_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIRS}
  ${PandoraSDK_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(MarlinMuonID
  PRIVATE
  ${RAIDA_LIBRARIES}
  ROOT::Hist
  ${ROOT_LIBRARIES}
  ${GSL_LIBRARIES}
  ${PandoraSDK_LIBRARIES}
  )
TARGET_COMPILE_DEFINITIONS(MarlinMuonID
  PRIVATE
  ${RAIDA_DEFINITIONS}
  ${GSL_DEFINITIONS}
  ${PandoraSDK_DEFINITIONS}
  )

INSTALL( TARGETS MarlinMuonID DESTINATION lib )
