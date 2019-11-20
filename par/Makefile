CC = g++

CFLAGS = -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -fopenmp
TARGET = calculo

all: $(TARGET)

$(TARGET): $(TARGET).cpp 
	$(CC) $(CFLAGS) $(TARGET).cpp nasteroids-par.cpp -o nasteroids-par 

clean:
	$(RM) $(TARGET)