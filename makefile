objects = main.o Initializer.o LavasimCal.o

all: $(objects)
	nvcc -arch=sm_30 $(objects) -o build -std=c++11

%.o: %.cpp
	nvcc -x cu -arch=sm_30 -I. -dc $< -o $@ -std=c++11

clean:
	rm -f *.o app
