CFLAGS = -g -idp /cygwin/ -ddp /cygwin/ -m 64 -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\amd64"
LDFLAGS = -L"C:\freeglut\lib\x64" -l"freeglut" -I"C:\freeglut\include" -I ""

%.o : %.cpp
	nvcc -x cu -c $^ $(CFLAGS) $(LDFLAGS) -o $@
	nvcc -M $^ $(CFLAGS) $(LDFLAGS) > $@.dep

%.o : %.c
	nvcc -c $^ $(CFLAGS) $(LDFLAGS) -o $@
	nvcc -M $^ $(CFLAGS) $(LDFLAGS) > $@.dep

%.o : %.cu
	nvcc -c $^ $(CFLAGS) $(LDFLAGS) -o $@
	nvcc -M $^ $(CFLAGS) $(LDFLAGS) > $@.dep

include $(wildcard *.dep) /dev/null

CPPLIST = $(wildcard *.cpp)
CLIST = $(wildcard *.c)
OBJECTS = $(CLIST:%.c=%.o) $(CPPLIST:%.cpp=%.o)
APP = MDSerial.exe

$(APP) : $(OBJECTS)
	nvcc $(OBJECTS) $(LDFLAGS) -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\amd64" -o $@

clean:
	$(RM) $(OBJECTS) *.dep *.pdb *.exe
