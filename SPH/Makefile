CXX = g++
CXXFLAGS = -std=c++17 -Wall -Iinclude

SRC = $(wildcard src/*.cpp)
OBJDIR = build
OBJ = $(patsubst src/%.cpp,$(OBJDIR)/%.o,$(SRC))

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET) src/*.o

