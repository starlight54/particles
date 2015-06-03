CFLAGS = -g
LDFLAGS = -l"freeglut" -I "" -I"/usr/gcc/include"

%.o : %.cpp
	gcc -x c++ $^ $(CFLAGS) $(LDFLAGS) $@
	gcc -M $^ $(CFLAGS) $(LDFLAGS) > $@.dep

%.o : %.c
	gcc $^ $(CFLAGS) $(LDFLAGS) $@
	gcc -M $^ $(CFLAGS) $(LDFLAGS) > $@.dep

include $(wildcard *.dep) /dev/null

CPPLIST = $(wildcard *.cpp)
CLIST = $(wildcard *.c)
OBJECTS = $(CLIST:%.c=%.o) $(CPPLIST:%.cpp=%.o)
APP = mdserial.o

$(APP) : $(OBJECTS)
	gcc $(OBJECTS) $(LDFLAGS) $@

clean:
	$(RM) $(OBJECTS) *.dep *.pdb *.exe
