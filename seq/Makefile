CC = g++

CFLAGS = -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3
TARGET = nasteroids-seq
AMIGO = calculo.cpp

all: $(TARGET)

$(TARGET): $(TARGET).cpp 
$(AMIGO): $(AMIGO).cpp 
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp $(AMIGO).cpp

clean:
	$(RM) $(TARGET)