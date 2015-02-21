CFLAGS = -g -idp /cygwin/ -ddp /cygwin/ -m 64

%.o : %.cpp
	nvcc -c $^ $(CFLAGS) -o $@
	nvcc -M $^ $(CFLAGS) > $@.dep

%.o : %.c
	nvcc -c $^ $(CFLAGS) -o $@
	nvcc -M $^ $(CFLAGS) > $@.dep

%.o : %.cu
	nvcc -c $^ $(CFLAGS) -o $@
	nvcc -M $^ $(CFLAGS) > $@.dep

include $(wildcard *.dep) /dev/null

CPPLIST = $(wildcard *.cpp)
CLIST = $(wildcard *.c)
OBJECTS = $(CLIST:%.c=%.o) $(CPPLIST:%.cpp=%.o)
APP = MDSerial.exe

$(APP) : $(OBJECTS)
	nvcc $(OBJECTS) $(LDFLAGS) -o $@

clean:
	$(RM) $(OBJECTS) *.dep *.pdb *.exe
