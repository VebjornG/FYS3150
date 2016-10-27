TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
#DEFINES += CATCH_CONFIG_MAIN

SOURCES += main.cpp \
    potential.cpp \
    systemsolver.cpp \
    unit_test.cpp

HEADERS += \
    potential.h \
    systemsolver.h \
    catch.hpp

LIBS += -llapack -lblas -larmadillo

