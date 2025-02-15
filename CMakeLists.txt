cmake_minimum_required(VERSION 3.16)
project(MoleculeViewer_Qt VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 COMPONENTS Core Gui Widgets OpenGL OpenGLWidgets REQUIRED)
find_package(OpenGL REQUIRED)

set(PROJECT_SOURCES
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    MoleculeGLWidget.h
    MoleculeGLWidget.cpp
    SDFParser.h
    SDFParser.cpp
    DataStructs.h
)

qt_add_executable(MoleculeViewer_Qt ${PROJECT_SOURCES})

target_link_libraries(MoleculeViewer_Qt PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
    ${OPENGL_LIBRARIES}
    GLU32
)

set_target_properties(MoleculeViewer_Qt PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER com.example.MoleculeViewer_Qt
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE 
)

include(GNUInstallDirs)
install(TARGETS MoleculeViewer_Qt
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
    add_custom_command(TARGET MoleculeViewer_Qt
        POST_BUILD
        COMMAND $<TARGET_FILE:Qt6::windeployqt> --dir $<TARGET_FILE_DIR:MoleculeViewer_Qt> $<TARGET_FILE:MoleculeViewer_Qt>
        COMMENT "Running windeployqt on $<TARGET_FILE:MoleculeViewer_Qt>"
    )
endif()

