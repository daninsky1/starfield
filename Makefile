ifeq ($(OS),Windows_NT)
    EXE_EXT = .exe
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S), MINGW32_NT)
        EXE_EXT = .exe
    else ifeq ($(UNAME_S), MINGW64_NT)
        EXE_EXT = .exe
    else
        EXE_EXT =
    endif
endif

TARGET = starfield$(EXE_EXT)

.PHONY: all run clean video gif

all: $(TARGET)

$(TARGET): FORCE
	gcc -O2 -std=c11 -pedantic -Wall -Wextra -mavx starfield.c -o $(TARGET) -fopenmp -lm
	
run:
	./$(TARGET)

clean:
	rm -f $(TARGET) ./field*.png
    
FORCE:
    