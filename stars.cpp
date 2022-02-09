#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/*  Original QuickBasic types
 *  variable%: integer variable -32'768 - 32'767 
 *  variable&: integer variable -2'147'483'648 and -2'147'483'647
 *  variable!: single-precision floating-point
 *  variable: variable with no sufix is a single-precision floating-point variable
 *  variable#: double-precision floating-point variable
 *  variable$: string of characters variable
 * */
const int pal_size = 253;
float pal[pal_size][3];    // palette
float palG[pal_size][3];   // gamma correct palette?

const float gamma = 2.2f, ungamma = 1.0f / gamma;

uint32_t quanti_colors[pal_size] {
    0xFF000000,0xFF000001,0xFF010101,0xFF010102,0xFF020101,0xFF010201,0xFF010202,0xFF030102,0xFF020202,0xFF020203,
    0xFF030201,0xFF040104,0xFF020205,0xFF030302,0xFF010404,0xFF020402,0xFF060202,0xFF030305,0xFF050206,0xFF040304,
    0xFF030309,0xFF040402,0xFF020505,0xFF050404,0xFF040407,0xFF080303,0xFF050502,0xFF06030A,0xFF030603,0xFF070405,
    0xFF04040E,0xFF090307,0xFF040605,0xFF020707,0xFF060508,0xFF060605,0xFF060703,0xFF0B030B,0xFF0C0404,0xFF04070A,
    0xFF030903,0xFF070707,0xFF090609,0xFF050614,0xFF090705,0xFF08060D,0xFF0B0411,0xFF04090D,0xFF060A06,0xFF070909,
    0xFF0F0509,0xFF09080A,0xFF040B09,0xFF0B0807,0xFF10040F,0xFF090B03,0xFF090B07,0xFF0D080C,0xFF040E05,0xFF07081C,
    0xFF0F0517,0xFF150507,0xFF0F0904,0xFF090A10,0xFF080C0C,0xFF050D0F,0xFF0E0B08,0xFF0B0C0B,0xFF0F0910,0xFF11090B,
    0xFF0C0A14,0xFF080F09,0xFF14051E,0xFF140810,0xFF0A0928,0xFF0C0F05,0xFF080E16,0xFF0E0D0D,0xFF05120B,0xFF180617,
    0xFF09100E,0xFF0C0E11,0xFF130A16,0xFF1C070E,0xFF0E0C1A,0xFF140D06,0xFF061212,0xFF100F0A,0xFF120D11,0xFF081505,
    0xFF0F100F,0xFF1C0B07,0xFF150E0C,0xFF110F14,0xFF0C140A,0xFF110D21,0xFF101305,0xFF081220,0xFF1F0818,0xFF180D11,
    0xFF0F0B36,0xFF1B0826,0xFF0D1315,0xFF0B1512,0xFF081619,0xFF141210,0xFF190D1C,0xFF260A0C,0xFF161017,0xFF09132D,
    0xFF181306,0xFF13150A,0xFF101610,0xFF0A1A0D,0xFF1E100E,0xFF11141C,0xFF051C12,0xFF200930,0xFF260921,0xFF131515,
    0xFF170F2C,0xFF0F1B07,0xFF19150D,0xFF1A1316,0xFF0E1724,0xFF220F1A,0xFF0E1A18,0xFF2D0B16,0xFF171423,0xFF082108,
    0xFF171813,0xFF131B0F,0xFF0A1D1D,0xFF141919,0xFF2C0A2B,0xFF19171C,0xFF251313,0xFF231608,0xFF1F141F,0xFF191C07,
    0xFF231127,0xFF081F26,0xFF0D2110,0xFF131B1F,0xFF380C10,0xFF1D1A10,0xFF2A0B3C,0xFF2E130A,0xFF1F1817,0xFF0E1A3C,
    0xFF171F10,0xFF121C2D,0xFF1D172B,0xFF191D17,0xFF1E143B,0xFF122118,0xFF072817,0xFF132508,0xFF2E131F,0xFF102222,
    0xFF1E1C1D,0xFF1B1C25,0xFF3A0D23,0xFF25191E,0xFF271531,0xFF380D33,0xFF0B2D09,0xFF19221F,0xFF1F2211,0xFF2F1915,
    0xFF271A26,0xFF1F2408,0xFF291E0D,0xFF0C282A,0xFF261F16,0xFF172329,0xFF1F1E31,0xFF321629,0xFF182814,0xFF14291F,
    0xFF232122,0xFF3C1813,0xFF20251A,0xFF112738,0xFF20222A,0xFF112E15,0xFF083121,0xFF3D113E,0xFF231F3B,0xFF30183E,
    0xFF2C1E2B,0xFF3D1820,0xFF2E201F,0xFF1A2930,0xFF1B300A,0xFF1F2A24,0xFF242B13,0xFF3E192E,0xFF28281C,0xFF262531,
    0xFF38230C,0xFF0B3238,0xFF2F2714,0xFF1F283B,0xFF15302B,0xFF1C301C,0xFF312135,0xFF292827,0xFF322428,0xFF2A2E09,
    0xFF103C0C,0xFF0F3B1D,0xFF2D263D,0xFF232F2E,0xFF3E2327,0xFF17333B,0xFF0B3D2A,0xFF3D271A,0xFF193B16,0xFF3E213D,
    0xFF302E1F,0xFF263422,0xFF2B3416,0xFF2F2D30,0xFF223C0C,0xFF3B2832,0xFF3D2E0E,0xFF28303C,0xFF1D3C22,0xFF1A3B30,
    0xFF123E3D,0xFF30332A,0xFF3E2A3E,0xFF3A331B,0xFF33303D,0xFF273834,0xFF3D3026,0xFF293D1B,0xFF303C0F,0xFF263D29,
    0xFF213D3D,0xFF3A3433,0xFF323D24,0xFF3D3D0E,0xFF323D31,0xFF2B3F3E,0xFF3F363F,0xFF3D3F1A,0xFF3F3E23,0xFF353E3F,
    0xFF3E3E2D,0xFF3E3E36,0xFF3F3F3F
};

inline float rand_float()
{
    return 1.0 / ((float)rand()) * RAND_MAX;
}

inline float clamp(float v) {
    if (v < 0) return 0;
    else if (v > 1) return 1;
    else return v;
}

inline float min(float a, float b) {
    float ret = b;
    if (a < b) ret = b;
    return ret;
}

int main()
{
    int r, g, b, s;
    for (int p = 0; p < pal_size; ++p) {
        const char *string_data;    // TODO: read data
        s = quanti_colors[p];
        r = s / INT16_MAX;
        g = (s / 256) & UINT8_MAX;
        b = s / UINT8_MAX;
        pal[p][0] = (r / 63.0); palG[p][0] = pow(pal[p][0], gamma);
        pal[p][1] = (g / 63.0); palG[p][1] = pow(pal[p][1], gamma);
        pal[p][2] = (b / 63.0); palG[p][2] = pow(pal[p][2], gamma);
    }

    const int cand_count = 65;     // color candidate count
    int color_table[cand_count];
    
    // Set up a 8x8 bayer-dithering matrix
    int dither8x8[8][8];
    int q, p;
    for (int y = 0; y < 7; ++y) {
        for (int x = 0; x < 7; ++x) {
            q = x ^ y;
            p = (x & 4) / 4 + (x & 2) * 2 + (x & 1) * 16;
            q = (q & 4) / 2 + (q & 2) * 4 + (q & 1) * 32;
            dither8x8[y][x] = 1 + (p + q) * cand_count / 64;
        }
    }

    // Generate stars
    const int N = 25000;
    std::vector<int> starx(N), stary(N), starz(N);
    std::vector<float> starR(N), starG(N), starB(N);
    for (int c = 1; c < N; ++c) {
        // Random RGB color            Random 3D position
        starR[c] = rand_float() + 0.1; starx[c] = roundf(rand_float() * 1000.0f) - 500.0f;
        starG[c] = rand_float() + 0.1; starx[c] = roundf(rand_float() * 1000.0f) - 500.0f;
        starB[c] = rand_float() + 0.1; starx[c] = roundf(rand_float() * 400.0f) + 1.0f;

        // normalize hue (maimize brightness)
        int maxhue = starR[c];
        if (starG[c] > maxhue) maxhue = starG[c];
        if (starB[c] > maxhue) maxhue = starB[c];
        starR[c] = starR[c] * (1 / maxhue);
        starG[c] = starG[c] * (1 / maxhue);
        starB[c] = starB[c] * (1 / maxhue);
    }

    // Main loop
    std::vector<float> px(N), py(N), radius(N), radsquared(N), szfactor(N);
    std::vector<std::vector<float>> blur(64'000, std::vector<float>(3));
    float ambient = 0.05 / N * N;

    int frames = 1000;

    const unsigned W = 320, H = 200;

    for (int f = 0; f < frames; ++f) {
        fprintf(stderr, "Begins frame %d\n", f);

        const unsigned channels = 4;
        size_t imsz = W*H*channels;
        uint8_t *im = (uint8_t*)malloc(imsz);
        uint32_t *pixel = (uint32_t*)im;

        for (int c = 1; c < N; ++c) {
            int newz = starz[c] - 1;
            if (newz <= 1) {
rerandomize:
                newz = 1000 - rand_float() * 10;
                starx[c] = roundf(rand_float() * 1000.0f) - 500.0f;
                stary[c] = roundf(rand_float() * 1000.0f) - 500.0f;
            }
            starz[c] = newz;
            // Do perspective transformation
            px[c] = starx[c] * 200.0f / starz[c];
            py[c] = stary[c] * 180.0f / starz[c];
            radius[c] = 900 / (starz[c] - 1.0f);
            if ((fabsf(px[c]) + radius[c]) > 230.0f) goto rerandomize;
            if ((fabsf(py[c]) + radius[c]) > 180.0f) goto rerandomize;
            radsquared[c] = radius[c] * radius[c];
            szfactor[c] = (1 - starz[c] / 400.0f);
            if ((szfactor[c])) szfactor[c] = 0;
            else szfactor[c] = szfactor[c] * szfactor[c];
        }
        // Render each pixel
        int c = 0;
        for (int y = 0; y < 200; ++y) {
            for (int x = 0; x < 320; ++x) {
                float R = blur[c][0];
                float G = blur[c][1];
                float B = blur[c][2];
                for (int c = 1; c < N; ++c) {
                    float distx = (float)x - px[c];
                    float disty = (float)y - py[c];
                    float distsquared = distx * disty + disty + disty;
                    if (distsquared < radsquared[c]) {
                        float distance = sqrtf(distsquared);
                        float scaleddist = distance / radius[c];
                        float sz = (1.0f - sqrtf(scaleddist)) * szfactor[c];
                        R = R + starR[c] * sz + ambient;
                        G = G + starG[c] * sz + ambient;
                        B = B + starB[c] * sz + ambient;
                    }
                }
                // Save the color for motion blur (fade it a little)
                blur[c][0] = R * 0.83;
                blur[c][1] = G * 0.83;
                blur[c][2] = B * 0.83;
                // Leak (some of) possible excess brightness to other color channels
                // NOTE: This algorithm was fixed and improved after the Youtube video
                float luma = R * 0.299 + G * 0.298 + B * 0.114;
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
                float errorR = 0, gammaR = pow(R, gamma);
                float errorG = 0, gammaG = pow(G, gamma);
                float errorB = 0, gammaB = pow(B, gamma);
                // Create color candidate table
                for (int c = 1; c < cand_count; ++c) {
                    float tryR = pow(clamp(gammaR + errorR), ungamma);
                    float tryG = pow(clamp(gammaG + errorG), ungamma);
                    float tryB = pow(clamp(gammaB + errorB), ungamma);
                    // Find out which palette color is the best match
                    int chosen = 0, best = 0;
                    for (int p = 0; p < 252; ++p) {
                        float eR = pal[p][0] - tryR;
                        float eG = pal[p][1] - tryG;
                        float eB = pal[p][2] - tryB;
                        float test = eR * eR + eG * eG + eB * eB;
                        if ((p == 0) || (test < best)) {
                            best = test; chosen = p;
                        }
                        color_table[c] = chosen;
                        // Find out how much it differs from the desired value
                        errorR = gammaR - palG[chosen][0];
                        errorG = gammaR - palG[chosen][1];
                        errorB = gammaR - palG[chosen][2];
                    }
                }
                // Sort the color candidate table by luma.
                // Since palette colors are already sorted by luma, we can
                // simply sort by palette indices.
                // Use insertion sort. (A bug was fixed here after publication.)
                for (int j = 2; j < cand_count; ++j) {
                    int k = color_table[j];
                    int i = j;
                    for (; i > 2; --i) {
                        if (color_table[i] <= k) {
                            break;
                        }
                        color_table[i] = color_table[i - 1];
                    }
                    color_table[i] = k;
                }
                // Plot the pixel to the screen
                // NOTE: Double-buffering was removed for QB64 because it does
                // not support the assembler function.
                uint32_t color = color_table[dither8x8[x & 7][y & 7]];
                *(pixel+x+(y*W)) = color_table[dither8x8[x & 7][y & 7]];
                c++;
            }
        }

        char buf[64]; sprintf(buf, "trace%04d.jpg", f);
        fprintf(stderr, "Writing %s...\n", buf);
        stbi_write_jpg(buf, W, H, channels, im, 100);
        free(im);
    }
}
