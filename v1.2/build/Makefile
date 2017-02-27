SNG_PATH = ../cases/single_profile
FILES = ../src/*.cpp \
 				../lib/*.cpp

CC = g++
CFLAGS = -std=c++11

single_profile:
	cp ${SNG_PATH}/src/main.cpp ../src/
	cp ${SNG_PATH}/src/pulsar_parameters.* ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: single_profile' >> ../src/log.txt
	mkdir -p ${SNG_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${SNG_PATH}/out/main.o
