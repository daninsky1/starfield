# Starfield C Port

This is a C port of Joel Yliluoma’s dithering algorithm demo [Starfield demo](https://bisqwit.iki.fi/story/howto/dither/jy/)

<p align="center">
  <img src="https://raw.githubusercontent.com/daninsky1/starfield/main/output.png" alt="Starfield output" width="600">
</p>

## About
I didn’t read the original paper at the time because I'am not a graphics programmer, it’s just a personal interest. I don't have the technical background to understand the math and theory behind it.  
For transparency, I hadn’t seen the C++ examples from the paper before or during the porting process. I skimmed through the paper while updating this README now and noticed them. Honestly, if I had read it back then, I might not have attempted the port at all.

My goal was to:
- Learn and tinker with C/C++.
- Practice porting code.
- Modify into:
    - Generate a standalone image output with no external dependencies (except `stb_image_write.h`)
    - Multithreading (partially successful)
    - I attempted to learn SIMD C/C++ intrinsics, but I couldn't achive any performace benefity.

This is not a real time demo, is the oposite, is still a pretty expensive implementation.

I've preserved the author comments in the code. I've added more simple didatic comments to help me with the migration.

The first iteration of the code was behaving weird, all stars was spawning to close to the camera, and after they passed the camera they respawn far away from the camera in a box. Those were fixed to spawn the stars evenly distributed.  

The computation still to much expensive even using multhreading, I need to analise if there are some caching issues, the threads seems to work fine, and there is a significant performance boost, maybe the code is just to unnoptimized.

SIMD operations was a faillure with no good resuts, I think that the branches in between SIMD operations kills any performance benefit, there was no improvement or regretion in performance.  
## Building

### Requirements

- A C compiler (GCC, Clang, MSVC — tested with MSVC 2019)
- `make` (for Linux/Unix) or `nmake` (for Windows, using `Makefile.vc`)

### Windows (MSVC)

```sh
nmake /f Makefile.vc
```

### MSYS2 and Linux
```sh
make
```
