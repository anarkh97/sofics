CXX = g++
CXXFLAGS = -O3 -fPIC
LDFLAGS = -shared

# Override to build against a different copy.
INCLUDE_DIR = ../../../src/solvers/m2c

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
