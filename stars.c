#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <omp.h>

#include <immintrin.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#if defined(_WIN32)
#include <direct.h>
#define mkdir(dir) _mkdir(dir)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define mkdir(dir) mkdir(dir, 0755)
#endif


/*  Fundamental QuickBasic types
 *  variable%: sined 16-bit integer
 *  variable&: sined 32-bit integer
 *  variable!: single-precision 32-bit floating-point
 *  variable: variable with no sufix is a single-precision 32-bit floating-point
 *  variable#: double-precision 64-bit floating-point
 *  variable$: string of characters
 */
static const int N = 25000;         // Number of stars
static const int PAL_SIZE = 253;
static const size_t csz = 3;        // Channel size, 3 colors
static float pal[PAL_SIZE][csz];    // palette
static float palG[PAL_SIZE][csz];   // gamma correct palette

static const float gamma = 2.2f, ungamma = 1.0f / gamma;

// Define an optimal palette (generated with NeuQuant). The colors are sorted
// by luminosity (for reasons that become apparent later).
// NOTE(daniel): 0x FF FF FF FF - 0x PP RR GG BB. PP: padding
static const uint32_t quanti_colors[PAL_SIZE] = {
    0x00000000, 0x00000001, 0x00010101, 0x00010102, 0x00020101, 0x00010201, 0x00010202, 0x00030102, 0x00020202, 0x00020203,
    0x00030201, 0x00040104, 0x00020205, 0x00030302, 0x00010404, 0x00020402, 0x00060202, 0x00030305, 0x00050206, 0x00040304,
    0x00030309, 0x00040402, 0x00020505, 0x00050404, 0x00040407, 0x00080303, 0x00050502, 0x0006030A, 0x00030603, 0x00070405,
    0x0004040E, 0x00090307, 0x00040605, 0x00020707, 0x00060508, 0x00060605, 0x00060703, 0x000B030B, 0x000C0404, 0x0004070A,
    0x00030903, 0x00070707, 0x00090609, 0x00050614, 0x00090705, 0x0008060D, 0x000B0411, 0x0004090D, 0x00060A06, 0x00070909,
    0x000F0509, 0x0009080A, 0x00040B09, 0x000B0807, 0x0010040F, 0x00090B03, 0x00090B07, 0x000D080C, 0x00040E05, 0x0007081C,
    0x000F0517, 0x00150507, 0x000F0904, 0x00090A10, 0x00080C0C, 0x00050D0F, 0x000E0B08, 0x000B0C0B, 0x000F0910, 0x0011090B,
    0x000C0A14, 0x00080F09, 0x0014051E, 0x00140810, 0x000A0928, 0x000C0F05, 0x00080E16, 0x000E0D0D, 0x0005120B, 0x00180617,
    0x0009100E, 0x000C0E11, 0x00130A16, 0x001C070E, 0x000E0C1A, 0x00140D06, 0x00061212, 0x00100F0A, 0x00120D11, 0x00081505,
    0x000F100F, 0x001C0B07, 0x00150E0C, 0x00110F14, 0x000C140A, 0x00110D21, 0x00101305, 0x00081220, 0x001F0818, 0x00180D11,
    0x000F0B36, 0x001B0826, 0x000D1315, 0x000B1512, 0x00081619, 0x00141210, 0x00190D1C, 0x00260A0C, 0x00161017, 0x0009132D,
    0x00181306, 0x0013150A, 0x00101610, 0x000A1A0D, 0x001E100E, 0x0011141C, 0x00051C12, 0x00200930, 0x00260921, 0x00131515,
    0x00170F2C, 0x000F1B07, 0x0019150D, 0x001A1316, 0x000E1724, 0x00220F1A, 0x000E1A18, 0x002D0B16, 0x00171423, 0x00082108,
    0x00171813, 0x00131B0F, 0x000A1D1D, 0x00141919, 0x002C0A2B, 0x0019171C, 0x00251313, 0x00231608, 0x001F141F, 0x00191C07,
    0x00231127, 0x00081F26, 0x000D2110, 0x00131B1F, 0x00380C10, 0x001D1A10, 0x002A0B3C, 0x002E130A, 0x001F1817, 0x000E1A3C,
    0x00171F10, 0x00121C2D, 0x001D172B, 0x00191D17, 0x001E143B, 0x00122118, 0x00072817, 0x00132508, 0x002E131F, 0x00102222,
    0x001E1C1D, 0x001B1C25, 0x003A0D23, 0x0025191E, 0x00271531, 0x00380D33, 0x000B2D09, 0x0019221F, 0x001F2211, 0x002F1915,
    0x00271A26, 0x001F2408, 0x00291E0D, 0x000C282A, 0x00261F16, 0x00172329, 0x001F1E31, 0x00321629, 0x00182814, 0x0014291F,
    0x00232122, 0x003C1813, 0x0020251A, 0x00112738, 0x0020222A, 0x00112E15, 0x00083121, 0x003D113E, 0x00231F3B, 0x0030183E,
    0x002C1E2B, 0x003D1820, 0x002E201F, 0x001A2930, 0x001B300A, 0x001F2A24, 0x00242B13, 0x003E192E, 0x0028281C, 0x00262531,
    0x0038230C, 0x000B3238, 0x002F2714, 0x001F283B, 0x0015302B, 0x001C301C, 0x00312135, 0x00292827, 0x00322428, 0x002A2E09,
    0x00103C0C, 0x000F3B1D, 0x002D263D, 0x00232F2E, 0x003E2327, 0x0017333B, 0x000B3D2A, 0x003D271A, 0x00193B16, 0x003E213D,
    0x00302E1F, 0x00263422, 0x002B3416, 0x002F2D30, 0x00223C0C, 0x003B2832, 0x003D2E0E, 0x0028303C, 0x001D3C22, 0x001A3B30,
    0x00123E3D, 0x0030332A, 0x003E2A3E, 0x003A331B, 0x0033303D, 0x00273834, 0x003D3026, 0x00293D1B, 0x00303C0F, 0x00263D29,
    0x00213D3D, 0x003A3433, 0x00323D24, 0x003D3D0E, 0x00323D31, 0x002B3F3E, 0x003F363F, 0x003D3F1A, 0x003F3E23, 0x00353E3F,
    0x003E3E2D, 0x003E3E36, 0x003F3F3F
};

typedef union ABGRColor {
    struct {
        uint8_t r;
        uint8_t g;
        uint8_t b;
        uint8_t a;
    };
    uint32_t c;
} ABGRColor; // FF FF FF FF - AA BB GG RR 

typedef struct Star {
    int x;
    int y;
    int z;
    float r;
    float g;
    float b;
    float px;
    float py;
    float radius;
    float radsquared;
    float szfactor;
} Star;

static inline float rand_float()
{
    return 1.0f / (float)RAND_MAX * (float)rand();
}

static inline float clamp(float v)
{
    if (v < 0.0f) return 0.0f;
    else if (v > 1.0f) return 1.0f;
    else return v;
}

static inline float min(float a, float b)
{
    if (a < b) return a;
    return b;
}

int main()
{
    uint32_t color_table[PAL_SIZE];

    for (int p = 0; p < PAL_SIZE; ++p) {
        uint8_t *c = (uint8_t*)&quanti_colors[p];
        uint8_t R = *(c+2);
        uint8_t G = *(c+1);
        uint8_t B = *(c);

        ABGRColor col;
        col.b = B;
        col.g = G;
        col.r = R;
        col.a = 0xFF;

        color_table[p] = col.c;
        
        pal[p][0] = ((float)R / 63.0f); palG[p][0] = powf(pal[p][0], gamma);
        pal[p][1] = ((float)G / 63.0f); palG[p][1] = powf(pal[p][1], gamma);
        pal[p][2] = ((float)B / 63.0f); palG[p][2] = powf(pal[p][2], gamma);
    }

    const int cand_count = 64;     // color candidate count
    
    // Set up a 8x8 bayer-dithering matrix
    int dither8x8[8][8];
    int q, p;
    for (int y = 0; y < 8; ++y) {
        for (int x = 0; x < 8; ++x) {
            q = x ^ y;
            p = (x & 4) / 4 + (x & 2) * 2 + (x & 1) * 16;
            q = (q & 4) / 2 + (q & 2) * 4 + (q & 1) * 32;
            dither8x8[y][x] = (p + q) * cand_count / 64;
        }
    }

    Star *stars = (Star*)malloc(sizeof(Star) * N);
    Star **vstars = (Star**)malloc(sizeof(Star*) * N);
    int vstars_size = N;

    const int W = 320, H = 200;
    int hW = W/2, hH = H/2;
    //float fieldx = 820.0f, fieldy = 820.0f, fieldz = 600.0f;
    float fieldx = 1000.0f, fieldy = 1000.0f, fieldz = 400.0f;
    int zclip = 1;
    int velocity = 2;
    for (int c = 0; c < N; ++c) {
        // Random RGB color              Random 3D position
        stars[c].r = rand_float() + 0.01f; stars[c].x = lrintf(rand_float() * fieldx) - fieldx/2;
        stars[c].g = rand_float() + 0.01f; stars[c].y = lrintf(rand_float() * fieldy) - fieldy/2;
        stars[c].b = rand_float() + 0.01f; stars[c].z = lrintf(rand_float() * fieldz) + zclip;
                                           //stars[c].z = (int)fieldz;

        // normalize hue (maximize brightness)
        float maxhue = stars[c].r;
        if (stars[c].g > maxhue) maxhue = stars[c].g;
        if (stars[c].b > maxhue) maxhue = stars[c].b;
        stars[c].r *= 1.0f / maxhue;
        stars[c].g *= 1.0f / maxhue;
        stars[c].b *= 1.0f / maxhue;
    }
    // Main loop
    // Blur buffer
    const int blur_sz = W * H;
    float *blur = (float*)malloc(sizeof(float) * csz * blur_sz);
    const float ambient = 0.05f / sqrtf((float)N);
    float cop = 0.0f;   // Center of Projection
    int vp = 400;       // Is Vanishing Point the right nomenclature?
    
    float *distsquaredA = (float*)malloc(sizeof(float) * N);
    float *distsqrt_gt_radsqrt = (float*)malloc(sizeof(float) * N);

    printf("Begins...\n");
    int frames = 1440; // 1 sec at 24 frames
    const unsigned channels = 4;
    for (int f = 0; f < frames; ++f) {
        double stime = omp_get_wtime();

        size_t imsz = W*H*channels;
        uint8_t *im = (uint8_t*)malloc(imsz);
        uint32_t *pixel = (uint32_t*)im;

        // Move each star
        int vsi = 0;     // valida stars index
        for (int s = 0; s < N; ++s) {
            int newz = stars[s].z - velocity;
            if (newz < zclip) {
                newz = (int)fieldz + lrintf(rand_float() * 10.0f);
                stars[s].x = lrintf(rand_float() * fieldx) - fieldx/2;
                stars[s].y = lrintf(rand_float() * fieldy) - fieldy/2;
            }
            stars[s].z = newz;
            // Do perspective transformation
            stars[s].px = (float)stars[s].x * 200.0f / (float)stars[s].z;
            stars[s].py = (float)stars[s].y * 180.0f / (float)stars[s].z;
            stars[s].radius = 900.0f / (float)(stars[s].z);

            stars[s].radsquared = stars[s].radius * stars[s].radius;
            stars[s].szfactor = (1.0f - (float)stars[s].z / fieldz);
            if ((stars[s].szfactor) < 0.0f) stars[s].szfactor = 0;
            else stars[s].szfactor *= stars[s].szfactor;

            if ((fabsf(stars[s].px) + stars[s].radius) > 230) continue;
            if ((fabsf(stars[s].py) + stars[s].radius) > 180) continue;
            //if (stars[s].z > vp) continue;

            vstars[vsi] = stars + s;
            vsi++;
        }
        vstars_size = vsi;

        // Render each pixel
        omp_set_num_threads(4);
#pragma omp parallel
        {
            for (int y = 0; y < H; ++y) {
#pragma omp for schedule(static)
                for (int x = 0; x < W; ++x) {
                    int bi = (y*W)+x;       // Buffer index
                    float R = blur[bi*csz];
                    float G = blur[bi*csz + 1];
                    float B = blur[bi*csz + 2];
//#pragma omp simd
                    float xf = (float)(x-hW), yf = (float)(y-hH);
#define SIMD_ 1
#if SIMD_ == 0
                    for (int s = 0; s < vstars_size; ++s) {
                        float distx = xf - vstars[s]->px;
                        float disty = yf - vstars[s]->py;
                        float distsquared = distx * distx + disty * disty;
                        if (distsquared < vstars[s]->radsquared) {
                            float distance = sqrtf(distsquared);
                            float scaleddist = distance / vstars[s]->radius;
                            float sz_sum = (1.0f - sqrtf(scaleddist)) * vstars[s]->szfactor + ambient;
                            R += vstars[s]->r * sz_sum;
                            G += vstars[s]->g * sz_sum;
                            B += vstars[s]->b * sz_sum;
                        }
                    }
#elif SIMD_ == 1
                    // NOTE(daniel): No speed gain in any of these SIMD versions.
    #define SIMD_VER 2
                    __m256 xf256 = _mm256_set1_ps(xf);
                    __m256 yf256 = _mm256_set1_ps(yf);

                    __m256 distx256, disty256, px256, py256;
                    __m256 distsq256, radsq256, cmp_sq256;

                    float r_result[8];
                    float g_result[8];
                    float b_result[8];
                    int test_ps[8];
                    for (int s = 7; s < vstars_size; s+=8) {
                        px256 = _mm256_set_ps(
                            vstars[s]->px,   vstars[s-1]->px, vstars[s-2]->px, vstars[s-3]->px,
                            vstars[s-4]->px, vstars[s-5]->px, vstars[s-6]->px, vstars[s-7]->px);
                        py256 = _mm256_set_ps(
                            vstars[s]->py,   vstars[s-1]->py, vstars[s-2]->py, vstars[s-3]->py,
                            vstars[s-4]->py, vstars[s-5]->py, vstars[s-6]->py, vstars[s-7]->py);
                        distx256 = _mm256_sub_ps(xf256, px256); // distx
                        disty256 = _mm256_sub_ps(yf256, py256); // disty

                        distsq256 = _mm256_add_ps(_mm256_mul_ps(distx256, distx256), _mm256_mul_ps(disty256, disty256));
                        float distances[8];
                        _mm256_storeu_ps(distances, distsq256);

                        radsq256 = _mm256_set_ps(
                            vstars[s]->radsquared,   vstars[s-1]->radsquared, vstars[s-2]->radsquared, vstars[s-3]->radsquared,
                            vstars[s-4]->radsquared, vstars[s-5]->radsquared, vstars[s-6]->radsquared, vstars[s-7]->radsquared);
                        // mask
                        cmp_sq256 = _mm256_cmp_ps(distsq256, radsq256, _CMP_LT_OQ);

    #if SIMD_VER == 0
                        // NOTE(daniel): This calculation causes loss of quality, I'm not sure
                        // where or why exactly, however, I think, it can  be the accumulation of
                        // the float imprecision, the branchless nature of this piece of code.
                        // This code branch deviates a few times to do the calculation, so a bunch
                        // of unessecery operations, and zero sums are maded.
                        __m256 scaleddist256, radius256, szsum256, const256, dist256;
                        __m256 szfactor256, r256, g256, b256, ambient256;

                        radius256 = _mm256_set_ps(
                            vstars[s]->radius,   vstars[s-1]->radius, vstars[s-2]->radius, vstars[s-3]->radius,
                            vstars[s-4]->radius, vstars[s-5]->radius, vstars[s-6]->radius, vstars[s-7]->radius);
                        scaleddist256 = _mm256_div_ps(/* distance */_mm256_sqrt_ps(distsq256), radius256);

                        szfactor256 = _mm256_set_ps(
                            vstars[s]->szfactor,   vstars[s-1]->szfactor, vstars[s-2]->szfactor, vstars[s-3]->szfactor,
                            vstars[s-4]->szfactor, vstars[s-5]->szfactor, vstars[s-6]->szfactor, vstars[s-7]->szfactor);
                        const256 = _mm256_set1_ps(1.0f);
                        ambient256 = _mm256_set1_ps(ambient);
                        szsum256 = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(const256, _mm256_sqrt_ps(scaleddist256)), szfactor256), ambient256);

                        const256 = _mm256_setzero_ps();
                        r256 = _mm256_set_ps(
                            vstars[s]->r,   vstars[s-1]->r, vstars[s-2]->r, vstars[s-3]->r,
                            vstars[s-4]->r, vstars[s-5]->r, vstars[s-6]->r, vstars[s-7]->r);
                        g256 = _mm256_set_ps(
                            vstars[s]->g,   vstars[s-1]->g, vstars[s-2]->g, vstars[s-3]->g,
                            vstars[s-4]->g, vstars[s-5]->g, vstars[s-6]->g, vstars[s-7]->g);
                        b256 = _mm256_set_ps(
                            vstars[s]->b,   vstars[s-1]->b, vstars[s-2]->b, vstars[s-3]->b,
                            vstars[s-4]->b, vstars[s-5]->b, vstars[s-6]->b, vstars[s-7]->b);
                        r256 = _mm256_blendv_ps(const256, _mm256_mul_ps(r256, szsum256), cmp_sq256);
                        g256 = _mm256_blendv_ps(const256, _mm256_mul_ps(g256, szsum256), cmp_sq256);
                        b256 = _mm256_blendv_ps(const256, _mm256_mul_ps(b256, szsum256), cmp_sq256);

                        const256 = _mm256_setzero_ps();
                        R += r256.m256_f32[0] + r256.m256_f32[1] + r256.m256_f32[2] + r256.m256_f32[3] + 
                            r256.m256_f32[4] + r256.m256_f32[5] + r256.m256_f32[6] + r256.m256_f32[7];
                        G += g256.m256_f32[0] + g256.m256_f32[1] + g256.m256_f32[2] + g256.m256_f32[3] + 
                            g256.m256_f32[4] + g256.m256_f32[5] + g256.m256_f32[6] + g256.m256_f32[7];
                        B += b256.m256_f32[0] + b256.m256_f32[1] + b256.m256_f32[2] + b256.m256_f32[3] + 
                            b256.m256_f32[4] + b256.m256_f32[5] + b256.m256_f32[6] + b256.m256_f32[7];
    #elif SIMD_VER == 1
                        // NOTE(daniel): This branch is to mitgate the problem of the code above,
                        // SIMD_VER 0, this makes more branch operation and less simd operation,
                        // The code is a mix of branch and branchless code. This has the same
                        // problem of SIMD_VER 0, however isolate alot of unecessary batches of 8.
                        // This code may be less efficient than SIMD_VER 0 if many stars need to
                        // calculate and if the distribution of stars is good
                        __m256i cmp_sq256i;
                        cmp_sq256i = _mm256_castps_si256(cmp_sq256);
                        if (cmp_sq256i.m256i_u32[0] | cmp_sq256i.m256i_u32[1] | cmp_sq256i.m256i_u32[2] | cmp_sq256i.m256i_u32[3] |
                            cmp_sq256i.m256i_u32[4] | cmp_sq256i.m256i_u32[5] | cmp_sq256i.m256i_u32[6] | cmp_sq256i.m256i_u32[7]) {
                            __m256 scaleddist256, radius256, szsum256, const256, dist256;
                            __m256 szfactor256, r256, g256, b256, ambient256;

                            radius256 = _mm256_set_ps(
                                vstars[s]->radius,   vstars[s-1]->radius, vstars[s-2]->radius, vstars[s-3]->radius,
                                vstars[s-4]->radius, vstars[s-5]->radius, vstars[s-6]->radius, vstars[s-7]->radius);
                            scaleddist256 = _mm256_div_ps(/* distance */_mm256_sqrt_ps(distsq256), radius256);

                            szfactor256 = _mm256_set_ps(
                                vstars[s]->szfactor,   vstars[s-1]->szfactor, vstars[s-2]->szfactor, vstars[s-3]->szfactor,
                                vstars[s-4]->szfactor, vstars[s-5]->szfactor, vstars[s-6]->szfactor, vstars[s-7]->szfactor);
                            const256 = _mm256_set1_ps(1.0f);
                            ambient256 = _mm256_set1_ps(ambient);
                            szsum256 = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(const256, _mm256_sqrt_ps(scaleddist256)), szfactor256), ambient256);

                            const256 = _mm256_setzero_ps();
                            r256 = _mm256_set_ps(
                                vstars[s]->r,   vstars[s-1]->r, vstars[s-2]->r, vstars[s-3]->r,
                                vstars[s-4]->r, vstars[s-5]->r, vstars[s-6]->r, vstars[s-7]->r);
                            g256 = _mm256_set_ps(
                                vstars[s]->g,   vstars[s-1]->g, vstars[s-2]->g, vstars[s-3]->g,
                                vstars[s-4]->g, vstars[s-5]->g, vstars[s-6]->g, vstars[s-7]->g);
                            b256 = _mm256_set_ps(
                                vstars[s]->b,   vstars[s-1]->b, vstars[s-2]->b, vstars[s-3]->b,
                                vstars[s-4]->b, vstars[s-5]->b, vstars[s-6]->b, vstars[s-7]->b);
                            r256 = _mm256_blendv_ps(const256, _mm256_mul_ps(r256, szsum256), cmp_sq256);
                            g256 = _mm256_blendv_ps(const256, _mm256_mul_ps(g256, szsum256), cmp_sq256);
                            b256 = _mm256_blendv_ps(const256, _mm256_mul_ps(b256, szsum256), cmp_sq256);

                            const256 = _mm256_setzero_ps();
                            R += r256.m256_f32[0] + r256.m256_f32[1] + r256.m256_f32[2] + r256.m256_f32[3] + 
                                r256.m256_f32[4] + r256.m256_f32[5] + r256.m256_f32[6] + r256.m256_f32[7];
                            G += g256.m256_f32[0] + g256.m256_f32[1] + g256.m256_f32[2] + g256.m256_f32[3] + 
                                g256.m256_f32[4] + g256.m256_f32[5] + g256.m256_f32[6] + g256.m256_f32[7];
                            B += b256.m256_f32[0] + b256.m256_f32[1] + b256.m256_f32[2] + b256.m256_f32[3] + 
                                b256.m256_f32[4] + b256.m256_f32[5] + b256.m256_f32[6] + b256.m256_f32[7];
                        }
    #elif SIMD_VER == 2
                        // NOTE(daniel): This code is to suppress the lost of quality caused by
                        // making this branchless.
                        // TODO(daniel): Make the index more readable.
                        __m256i cmp_sq256i = _mm256_castps_si256(cmp_sq256);
                        uint32_t vals[8];
                        _mm256_storeu_si256((__m256i*)vals, cmp_sq256i);
                        for (int i = 0; i < 8; ++i) {
                            if (vals[i]) {
                                float distance = sqrtf(distances[i]);
                                float scaleddist = distance / vstars[s-7+i]->radius;
                                float sz_sum = (1.0f - sqrtf(scaleddist)) * vstars[s-7+i]->szfactor + ambient;
                                
                                R += vstars[s-7+i]->r * sz_sum;
                                G += vstars[s-7+i]->g * sz_sum;
                                B += vstars[s-7+i]->b * sz_sum;
                            }
                        }
    #endif
                    }
#endif
                    blur[bi*csz]     = R * 0.83f;
                    blur[bi*csz + 1] = G * 0.83f;
                    blur[bi*csz + 2] = B * 0.83f;
                    // Leak (some of) possible excess brightness to other color channels
                    // NOTE: This algorithm was fixed and improved after the Youtube video
                    float luma = R * 0.299f + G * 0.298f + B * 0.114f;
                    if (luma >= 1.0f) {
                        R = 1.0f; G = 1.0f; B = 1.0f;
                    }
                    else if (luma <= 0.0f) {
                        R = 0.0f; G = 0.0f; B = 0.0f;
                    }
                    else {
                        float sat = 1.0f;
                        if (R > 1.0f) {
                            sat = min(sat, (luma - 1.0f)/ (luma - R));
                        }
                        else if (R < 0.0f) {
                            sat = min(sat, luma / (luma - R));
                        }
                        if (G > 1.0f) {
                            sat = min(sat, (luma - 1.0f)/ (luma - G));
                        }
                        else if (G < 0.0f) {
                            sat = min(sat, luma / (luma - G));
                        }
                        if (B > 1.0f) {
                            sat = min(sat, (luma - 1.0f)/ (luma - B));
                        }
                        else if (B < 0.0f) {
                            sat = min(sat, luma / (luma - B));
                        }
                        if (sat < 1.0f) {
                            R = (R - luma) * sat + luma;
                            G = (G - luma) * sat + luma;
                            B = (B - luma) * sat + luma;
                        }
                    }
                    // Quantize (use gamma-aware Knoll-Yliluoma positional dithering)
                    float errorR = 0, gammaR = powf(R, gamma);
                    float errorG = 0, gammaG = powf(G, gamma);
                    float errorB = 0, gammaB = powf(B, gamma);
                    // Create color candidate table
                    int cand_list[cand_count];
                    for (int c = 0; c < cand_count; ++c) {
                        float tryR = powf(clamp(gammaR + errorR), ungamma);
                        float tryG = powf(clamp(gammaG + errorG), ungamma);
                        float tryB = powf(clamp(gammaB + errorB), ungamma);

                        // Find out which palette color is the best match
                        int chosen = 0;
                        float best = 0.0f;
                        for (int p = 0; p < 253; ++p) {
                            float eR = pal[p][0] - tryR;
                            float eG = pal[p][1] - tryG;
                            float eB = pal[p][2] - tryB;
                            float test = eR * eR + eG * eG + eB * eB;
                            if ((p == 0) || (test < best)) {
                                best = test; chosen = p;
                            }
                        }
                        cand_list[c] = chosen;
                        // Find out how much it differs from the desired value
                        errorR = gammaR - palG[chosen][0];
                        errorG = gammaG - palG[chosen][1];
                        errorB = gammaB - palG[chosen][2];
                    }
                    // Sort the color candidate table by luma.
                    // Since palette colors are already sorted by luma, we can
                    // simply sort by palette indices.
                    // Use insertion sort. (A bug was fixed here after publication.)
                    for (int j = 1; j < cand_count; ++j) {
                        unsigned k = cand_list[j];
                        int i = j;
                        for (; i >= 1; --i) {
                            if (cand_list[i] <= k) break;
                            cand_list[i] = cand_list[i - 1];
                        }
                        cand_list[i] = k;
                    }
                    // Plot the pixel to the screen
                    // NOTE: Double-buffering was removed for QB64 because it does
                    // not support the assembler function.
                    ABGRColor color;
                    color.c = color_table[cand_list[dither8x8[x & 7][y & 7]]];
                    float brightness_boost = 2.0f;
                    float r_bb = (float)color.r * brightness_boost;
                    float g_bb = (float)color.g * brightness_boost;
                    float b_bb = (float)color.b * brightness_boost;
                    // Clamp > 255.0f
                    if (r_bb > 255.0f) r_bb = 255.0f;
                    if (g_bb > 255.0f) g_bb = 255.0f;
                    if (b_bb > 255.0f) b_bb = 255.0f;
                    
                    color.r = (uint8_t)r_bb;
                    color.g = (uint8_t)g_bb;
                    color.b = (uint8_t)b_bb;

                    *(pixel+bi) = color.c;
                }
            }
        }

        double etime = omp_get_wtime();
        mkdir("output");
        char buf[64]; sprintf(buf, "output/field%04d.jpg\n", f);
        fprintf(stderr, "\"%s\" - %d in %f seconds.", buf, f, etime-stime);
        stbi_write_jpg(buf, W, H, channels, im, 100);
        free(im);
    }
}
