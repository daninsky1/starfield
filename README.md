# Starfield C Port

This is port of [Joel Yliluoma’s QuickBasic dithering demo](https://bisqwit.iki.fi/story/howto/dither/jy/) to C, adding image export, multithreading.  

<p align="center">
  <img src="https://raw.githubusercontent.com/daninsky1/starfield/main/output.png" alt="Starfield output" width="600">
</p>  

I added the ability to export the final result as PNG or any format supported by the `stb_image_write` library.  
Multithreading support using OpenMP was also added.  
I even attempted to use SIMD operations didn’t yield any noticeable improvement. I suspect that branching between SIMD operations negates any benefits, resulting in no performance gain.  

The computation is still quite expensive even with multithreading. I need to analyze potential cache line optmization oportunities.  

At the time, I hadn’t read the paper and didn’t know he had already implemented this algorithm in C/C++. In hindsight, the port may not have been strictly necessary, but it makes no difference, this project was just for fun.
## Building

### Windows (MSVC)

```batch
nmake /f Makefile.vc
```

### MSYS2 and Linux
```sh
make
```
