CXX = g++
CXXFLAGS = -std=c++17 -Wall -Iinclude

# Archivos fuente (todos los .cpp dentro de src + main.cpp en raíz)
SRC = $(wildcard src/core/*.cpp src/io/*.cpp src/physics/*.cpp src/utils/*.cpp) main.cpp

OBJDIR = build
OBJ = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC))

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Regla genérica para compilar .cpp a .o
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET) test_results/*.output output/*.png output/*.txt test_results/ build/
