C = g++
FLAGS = -std=c++17 -Wall -Wextra -pedantic

SRC = lab6.cpp
TARGET = lab6

$(TARGET): $(SRC)
	$(C) $(FLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)
