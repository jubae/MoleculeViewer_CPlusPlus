cmake_minimum_required(VERSION 3.5)  # 3.5 is a bit old; consider 3.16+

project(MoleculeViewer_Qt VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- FIND PACKAGES (Corrected) ---
find_package(Qt6 COMPONENTS Core Gui Widgets OpenGL OpenGLWidgets REQUIRED) # Use Qt6
find_package(OpenGL REQUIRED)  # <--- THIS WAS MISSING!

set(PROJECT_SOURCES
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    MoleculeFactory.h
    MoleculeFactory.cpp
    MoleculeGLWidget.h
    MoleculeGLWidget.cpp
    SDFParser.h
    SDFParser.cpp
    DataStructs.h
)

# --- Qt6 Add Executable ---
qt_add_executable(MoleculeViewer_Qt
    ${PROJECT_SOURCES}
)

# --- Link Libraries (Corrected) ---
target_link_libraries(MoleculeViewer_Qt PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
    ${OPENGL_LIBRARIES}  # Now correctly used
    GLU32              # Add for Windows (harmless on other platforms)
)

# --- Bundle ID (Optional) ---
if(QT_VERSION LESS 6.1.0)
    set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.MoleculeViewer_Qt)
endif()

set_target_properties(MoleculeViewer_Qt PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

# --- Installation (Optional) ---
include(GNUInstallDirs)
install(TARGETS MoleculeViewer_Qt
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

#qt_finalize_executable(MoleculeViewer_Qt) #This isnt needed