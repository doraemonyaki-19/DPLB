#CC = pgcc
#CC = gcc
CC = icc

#pgcc acc profile#CFLAGS = -DMKL_ILP64 -I${MKLROOT}/include -O0 -Minform=inform -Minfo -acc -acclibs -ta=multicore -g -pg -Mprof
#pgcc omp profile#CFLAGS = -DMKL_ILP64 -I${MKLROOT}/include -O0 -Minform=inform -Minfo -g -pg -Mprof -mp

#CFLAGS = -DMKL_ILP64 -I${MKLROOT}/include -std=gnu11 -O0 -Wall -g -qno-openmp -xHost#-qno-openmp #icc omp profile #static
CFLAGS = -I${MKLROOT}/include -std=gnu11 -O3 -Wall -ipo -qopenmp   #-ipo #-g for profile #icc omp profile #dynamic

#gcc omp profile#
#CFLAGS = -m64 -I${MKLROOT}/include -std=gnu11 -O2 -fipa-pta -flto -Wall -Wextra -fopenmp -march=native #-g -pg for profile

#LDFLAGS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group#static
LDFLAGS = -L${MKLROOT}/lib/intel64#dynamic

#gcc
#LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed

#LDLIBS = -lpthread -lm -ldl #static
LDLIBS = -lmkl_rt -lpthread -lm -ldl#dynamic

C_OBJ  = main.o driver.o update.o lbe_update.o \
	bnodes.o bnodes_init.o bnodes_dp.o bnodes_tube.o \
	implicit_force.o lbe.o modes_write.o initlbe.o \
	lub.o velcs_update.o hs3d.o n_list.o objects_init.o objects_map.o \
 	clusters.o cluster_force.o cluster_update.o cj_grad.o \
	global_sums.o utils.o output.o \
	init_sphere.o init_polymer.o init_cylinder.o ran_num.o get_forces.o verlet_update.o 

clean:
	rm -f *.o aslb
	rm -f ./data/*

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS) $(LDLIBS)

serial: $(C_OBJ)
	$(CC) $(CFLAGS) -c msg_ser.c $(LDFLAGS) $(LDLIBS)
	$(CC) $(CFLAGS) -o aslb $(C_OBJ) msg_ser.o $(LDFLAGS) $(LDLIBS)

#test:
#	$(CC) $(CFLAGS) -c lbe.c $(LDFLAGS) $(LDLIBS)
