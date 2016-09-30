BETA_PATH = ../cases/beta_change
BT_PATH = ../cases/Bt_change
RM_PATH = ../cases/r_em-change
SNG_PATH = ../cases/single_profile
LMBD_PATH = ../cases/lambda_change
GAM_PATH = ../cases/gamma_change
INT_PATH = ../cases/interpulse
FILES = ../src/*.cpp \
 				../lib/*.cpp

CC = g++
CFLAGS = -std=c++11

beta_change:
	cp ${BETA_PATH}/src/mainWIN.cpp ../src/
	cp ${BETA_PATH}/src/{pulsar_parameters.cpp,pulsar_parameters.h} ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: beta_change' >> ../src/log.txt
	mkdir -p ${BETA_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${BETA_PATH}/out/main.o

bt_change:
	cp ${BT_PATH}/src/mainWIN.cpp ../src/
	cp ${BT_PATH}/src/{pulsar_parameters.cpp,pulsar_parameters.h} ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: bt_change' >> ../src/log.txt
	mkdir -p ${BT_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${BT_PATH}/out/main.o

rm_change:
	cp ${RM_PATH}/src/mainWIN.cpp ../src/
	cp ${RM_PATH}/src/{pulsar_parameters.cpp,pulsar_parameters.h} ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: rm_change' >> ../src/log.txt
	mkdir -p ${RM_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${RM_PATH}/out/main.o

single_profile:
	cp ${SNG_PATH}/src/mainWIN.cpp ../src/
	cp ${SNG_PATH}/src/{pulsar_parameters.cpp,pulsar_parameters.h} ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: single_profile' >> ../src/log.txt
	mkdir -p ${SNG_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${SNG_PATH}/out/main.o

lambda_change:
	cp ${LMBD_PATH}/src/mainWIN.cpp ../src/
	cp ${LMBD_PATH}/src/{pulsar_parameters.cpp,pulsar_parameters.h} ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: lambda_change' >> ../src/log.txt
	mkdir -p ${LMBD_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${LMBD_PATH}/out/main.o

gamma_change:
	cp ${GAM_PATH}/src/mainWIN.cpp ../src/
	cp ${GAM_PATH}/src/{pulsar_parameters.cpp,pulsar_parameters.h} ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: gamma_change' >> ../src/log.txt
	mkdir -p ${GAM_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${GAM_PATH}/out/main.o

interpulse:
	cp ${INT_PATH}/src/mainWIN.cpp ../src/
	cp ${INT_PATH}/src/{pulsar_parameters.cpp,pulsar_parameters.h} ../lib/
	rm -f ../src/log.txt
	touch ../src/log.txt
	echo 'CURRENT: gamma_change' >> ../src/log.txt
	mkdir -p ${INT_PATH}/out
	${CC} ${CFLAGS} ${FILES} -o ${INT_PATH}/out/main.o
