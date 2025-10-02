CXX = g++
CXXFLAGS = -std=c++17 -Wall -I../SPH/include

TARGET = memtest
SRC = main.cpp
OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET) vaciado_tanque.txt

clean:
	rm -f $(OBJ) $(TARGET)
