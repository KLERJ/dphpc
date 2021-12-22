include ../../config.mk

EXTRA_FLAGS=
BIN_DIR=./

heat-3d: heat-3d.c heat-3d.h
	${VERBOSE} ${CC} -o ${BIN_DIR}/heat-3d heat-3d.c ${CFLAGS} -I. -I../../utilities ../../utilities/polybench.c ${EXTRA_FLAGS}

clean:
	@ rm -f heat-3d

