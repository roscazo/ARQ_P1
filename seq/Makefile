CC = g++

CFLAGS = -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3
TARGET = nasteroids-seq

all: $(TARGET)

$(TARGET): $(TARGET).cpp 
	$(CC) $(CFLAGS) $(TARGET).cpp calculo.cpp -o $(TARGET) 

clean:
	$(RM) $(TARGET)