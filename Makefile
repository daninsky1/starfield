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
	g++ -O2 -std=c2x -pedantic -Wall -Wextra -mavx stars.c -o $(TARGET) -fopenmp
	
run:
	./$(TARGET)

clean:
	rm -f $(TARGET) ./field*.jpg
    
# video:
#     ffmpeg -framerate 30 -i output/field%04d.jpg -c:v libx264 -pix_f mt yuv420p demo.mp4
    
# gif:
#     ffmpeg -framerate 10 -start_number 15 -i output/field%04d.jpg -vframes 30 -loop -0 demo.gif
    
FORCE:
    