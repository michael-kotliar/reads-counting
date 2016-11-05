######################################################################
# pro file for complexity plot Wed Oct 19 14:59:16 2011
######################################################################

TEMPLATE = app
TARGET   = reads-counting

CONFIG   += console warn_on release
CONFIG   -= app_bundle

QT       -= gui


DEFINES     += D_USE_BAM      \
               _NORMAL


!win32{

DEFINES        += _APPNAME=\\\"$$TARGET\\\"
LIBS           += -lm -lz ../bamtools/libbamtools.a -lz
#QMAKE_CXXFLAGS += -Werror 
#-std=c++11

lib_bamtools.commands = cd ../bamtools/; qmake; $(MAKE) -j 8
QMAKE_EXTRA_TARGETS   = lib_bamtools
PRE_TARGETDEPS        = lib_bamtools

OBJECTS_DIR = GeneratedFiles
UI_DIR      = GeneratedFiles
MOC_DIR     = GeneratedFiles
RCC_DIR     = GeneratedFiles

INCLUDEPATH += /usr/local/include/

}

win32{

DEFINES        += _APPNAME=\"$$TARGET\"
LIBS           += -lbamtools

}

macx{

QMAKE_CFLAGS_X86_64 += -mmacosx-version-min=10.7
QMAKE_CXXFLAGS_X86_64 = $$QMAKE_CFLAGS_X86_64

}



DEPENDPATH +=  . \
               ./src \
               ../bamtools

INCLUDEPATH += . \
               ./src \
               ../bamtools



HEADERS     +=  ./src/Arguments.hpp \
                ./src/Reads.hpp \
		./src/Threads.hpp \
		./src/ReadsCounting.hpp \
                ./src/Matrix.hpp \
                ./src/Math.hpp
		
               
               
SOURCES     += \
                ./src/Arguments.cpp \
                ./src/Reads.cpp \
                ./src/Threads.cpp \
                ./src/ReadsCounting.cpp \
                ./src/main.cpp


QMAKE_CLEAN += $$TARGET logfile.log *~ ./src/*~ *.txt ../global/*~

