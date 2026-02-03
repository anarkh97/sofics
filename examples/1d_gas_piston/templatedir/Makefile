CXX = gcc
CXXFLAGS = -O3 -fPIC
LDFLAGS = -shared

INCLUDE_DIR = /home/anarkhede/tinkercliffs/m2c_global/

TARGET = StateCalculator.so
OBJECTS = StateCalculator.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $(TARGET) 
	rm $(OBJECTS)

StateCalculator.o: StateCalculator.cpp
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) StateCalculator.cpp -o StateCalculator.o

clean:
	rm $(TARGET)
