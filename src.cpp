#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <glm.hpp>
#include <glad/glad.h>
#include <glfw3.h>

#include <Windows.h>
#include <WinUser.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"

const int pixels_x = 1920;
const int pixels_y = 1080;

float mouseWheeloffset = 0;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }
}

namespace helpers {

    struct ShaderReader {

        ShaderReader(const char* path) : _path(path) {}

        const char* read() {
            std::ifstream shaderFile(_path);
            std::stringstream source;
            source << shaderFile.rdbuf();
            shaderFile.close();
            return source.str().c_str();
        }

        void update(const char* path) { _path = path; }

    private:

        std::string _path;

    };

    template<typename T>
    T calculate_scaling_factor(T p1, T p2, T q1, T q2) {
        return ((p2 - p1)) / ((q2 - q1));
    }

    template<typename T>
    T calculate_scaled(T val, T scaling_factor, T median) {
        return (val * scaling_factor) + median;
    }

    void takeScreenshot() {
        // get the device context of the screen
        HDC hScreenDC = CreateDC(L"DISPLAY", NULL, NULL, NULL);
        // and a device context to put it in
        HDC hMemoryDC = CreateCompatibleDC(hScreenDC);

        int width = GetDeviceCaps(hScreenDC, HORZRES);
        int height = GetDeviceCaps(hScreenDC, VERTRES);

        // maybe worth checking these are positive values
        HBITMAP hBitmap = CreateCompatibleBitmap(hScreenDC, width, height);

        // get a new bitmap
        HBITMAP hOldBitmap = (HBITMAP)SelectObject(hMemoryDC, hBitmap);

        BitBlt(hMemoryDC, 0, 0, width, height, hScreenDC, 0, 0, SRCCOPY);
        hBitmap = (HBITMAP)SelectObject(hMemoryDC, hOldBitmap);

        // clean up
        DeleteDC(hMemoryDC);
        DeleteDC(hScreenDC);
    }
}

#define color_MAX 256

unsigned int colorPalette[color_MAX];

void fillcolorp() {
    int colors[] = {
        0xFFFFFF,
        0xEBF5FF,
        0xD6EBFF,
        0xC2E1FF,
        0xADD8FF,
        0x96CDFF,
        0x85c4ff,
        0x70BAFF,
        0x5CB0ff,
        0x47a6ff,
        0x339cff,
        0x1f93ff,
        0x0a89ff,
        0x007ef5
    };
    colorPalette[0] = colors[0];
    colorPalette[1] = colors[1];
    colorPalette[2] = colors[2];
    colorPalette[3] = colors[3];
    colorPalette[4] = colors[4];
    colorPalette[5] = colors[5];
    colorPalette[6] = colors[6];
    colorPalette[7] = colors[7];
    colorPalette[8] = colors[8];
    colorPalette[9] = colors[9];
    colorPalette[10] = colors[10];
    colorPalette[11] = colors[11];
    colorPalette[12] = colors[12];
    colorPalette[13] = colors[13];
    for (unsigned int i = 13; i < 256; i++) {
        colorPalette[i] = 0;
    }
}

void fillcolorp2() {
    unsigned int color = 0xFFFFFF;
    unsigned int i = 0;
    unsigned int depr = 0x030200;
    while (i < color_MAX) {
        color -= depr;
        if (color - depr <= 0x00FFFF) {
            depr = 0x000201;
            if (color - depr <= 0x000100) {
                depr = 1;
            }
        }
        colorPalette[i] = color;
        i++;
    }
}

void fillcolorp3() {
    unsigned int i = 0;
    unsigned int c = 0x000000;
    while (i < color_MAX) {
        colorPalette[i] = c;
        i++;
        c += 0x000100;
    }
}

void fillcolorp4() {
    for (unsigned i = 0; i < color_MAX; i++) {
        colorPalette[i] = (i * 550) + (std::rand() % 550);
    }
}

void fillcolorp5() {
    unsigned int i = 0;
    unsigned int c = 0xFFFFFF;
    while (i < color_MAX) {
        colorPalette[i] = c;
        i++;
        if (i % 25 == 0) {
            c -= 0x010101;
        }
    }
}

void fillcolorp_grey() {
    unsigned c = 0xFFFFFF;
    unsigned i = 255;
    while (true) {
        colorPalette[i] = c;
        c -= 0x010101;
        i--;
        if (i == 0) break;
    }
}

void printContentsOfCP() {
    for (unsigned int i = 0; i < 256; i++) {
        std::cout << std::hex << "0x" << colorPalette[i] << '\n';
    }
}

void fillTexture(unsigned char* data, unsigned int i, int w, int h, int p) {
    for (unsigned int c = 0; c < w * h * p;) {
        unsigned int color = colorPalette[i];
        data[c] = (color >> 0x10) & 0xFF;
        data[c + 1] = (color >> 0x8) & 0xFF;
        data[c + 2] = color & 0xFF;
        c += 3;
    }
}

#include <chrono>
#include <ctime>

void generatePaletteTexture(unsigned int     *data, int size) {
    float red_depr = (float)0xFF / (float)size * 1.5f;
    float green_depr = (float)0xff / (float)size;
    float blue_depr = (float)0xFF / (float)size / (float)2;
    float red_c = 0xFF;
    float green_c = 0xFF;
    float blue_c = 0xFF;
    unsigned int c = 0;
    for (unsigned int px = 0; px < size; px++) {
        if (red_c <= 0.f) {
            red_c = 0;
            red_depr = 0;
        } else if (green_c <= 0.f) {
            green_c = 0;
            green_depr = 0;
        } else if (blue_c <= 0.f) {
            blue_c = 0;
            blue_depr = 0;
        }
        unsigned int color = (unsigned int)red_c << 0x10;
        color = color | ((unsigned int)green_c << 0x8);
        color = color | 0xff;
        red_c -= red_depr;
        green_c -= green_depr;
        blue_c -= blue_depr;
        data[c] = color;
        c++;
    }
}

void colorsquares(unsigned char* data, int w, int h) {
    unsigned int color = 0x000000;
    unsigned int sdim = w / 8;
    unsigned int c = 0;
    unsigned int bt = 0;
    for (unsigned int i = 0; i < h; i++) {
        unsigned int old_color = color;
        for (unsigned int j = 0; j < w; j += sdim) {
            for (unsigned int k = 0; k < sdim; k++) {
                data[c]     = (color >> 0x10) & 0xFF;
                data[c + 1] = (color >> 0x8) & 0xFF;
                data[c + 2] = color & 0xFF;
                c += 3;
            }
            color += 0x010204;
            printf("color: %d\n", color);
        }
        bt++;
        if (bt % sdim != 0) {
            color = old_color;
        } else {
            bt = 0;
        }
    }
}

void mandelbrotSet(unsigned char* data, int w, int h) {

    unsigned int size = w * h * 3;
    const double mx2 = 1;        // Mandelbrot X scale (-2.5, 1)
    const double mx1 = -2.5;
    const double my2 = 1;        // Mandelbrot Y scale (-1, 1)
    const double my1 = -1;
    double scalingFactorx = helpers::calculate_scaling_factor(mx1, mx2, 0.0, static_cast<double>(w));
    double scalingFactory = helpers::calculate_scaling_factor(my1, my2, 0.0, static_cast<double>(h));

    unsigned int c = 0;
    for (unsigned int ph = 0; ph < pixels_y; ph++) {
        double py = helpers::calculate_scaled(static_cast<double>(ph + 0.5), scalingFactory, my1);
        for (unsigned int pw = 0; pw < pixels_x; pw++) {
            double px = helpers::calculate_scaled(static_cast<double>(pw + 0.5), scalingFactorx, mx1);
            double x = 0.f;
            double y = 0.f;
            double x2 = 0.f;
            double y2 = 0.f;
            int i = 0;
            int max_i = color_MAX;
            while (x2 + y2 <= 4 && i < max_i) {
                y = ((x + x) * y) + py;
                x = x2 - y2 + px;
                x2 = x * x;
                y2 = y * y;
                i++;
            }   
            int color = colorPalette[i];
            data[c]     =   color >> 0x10;
            data[c + 1] =   color >> 0x8;
            data[c + 2] =   color;
            c += 3;
        }
    }
}


template<typename T>
T linear_interpolation(T a, T b, T c) {
    return a + ((1/3) * (b - a - c));
}

#define SIZE_PX (pixels_x)

struct pixel_array {
    float px[SIZE_PX];
};
namespace {
    void calculate_pxo(unsigned int w, pixel_array & pa) {
        float scalingFactorx = helpers::calculate_scaling_factor(-2.5f, 1.f, 0.f, static_cast<float>(w));
        for (unsigned int i = 0; i < w; i++) {
            pa.px[i] = helpers::calculate_scaled(static_cast<float>(i), scalingFactorx, -2.5f);
        }
    }

    bool calculate_line_color(unsigned int px, unsigned int py, unsigned int len, unsigned char direction) {
        // Mandelbrot algorithm to calculate pixel color using continous coloring
        if (direction == 2); // go down
        if (direction == 4); // go left
        if (direction == 6); // go right
        if (direction == 8); // go up
        if (px % 2 == 0) return true;       // placeholder, ofc replace this
        return false;
    }

    void cc_mandelbrot_rect(unsigned char* data, int w, int h, pixel_array & pa) {
        float log2 = std::logf(2);
        unsigned int c = 0;
        float cond1 = 1 << 16;
        float scalingFactory = helpers::calculate_scaling_factor(-1.f, 1.f, 0.f, static_cast<float>(h));
        auto start = std::chrono::steady_clock::now();
        unsigned int color;
        int de_iter = 0;

        bool squareBoxes = true;

        unsigned int r_x = 1;
        unsigned int r_y = 1;

        unsigned int p_x = 0;
        unsigned int p_y = 0;

        // Boxes are 25x25 pixels large
        while (squareBoxes) {

            bool box_finished = false;
            while (!box_finished) {
                auto equalcolor = calculate_line_color(r_x * 25, r_y * 25, 25, 4);
                if (!equalcolor) {

                }
                equalcolor = calculate_line_color((r_x - 1) * 25, r_y * 25, 25, 2);
                equalcolor = calculate_line_color((r_x - 1) * 25, (r_y - 1) * 25, 25, 6);
                equalcolor = calculate_line_color(r_x * 25, (r_y - 1) * 25, 25, 8);
            }

        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\t" << de_iter << " log calculations performed!\n";
    }

    void cc_mandelbrot(unsigned char* data, int w, int h, pixel_array & pa) {
        float log2 = std::logf(2);
        unsigned int c = 0;
        float cond1 = 1 << 16;
        float scalingFactory = helpers::calculate_scaling_factor(-1.f, 1.f, 0.f, static_cast<float>(h));
        auto start = std::chrono::steady_clock::now();
        unsigned int color;
        int de_iter = 0;
        for (unsigned int ph = 0; ph < h; ph++) {
            float py = helpers::calculate_scaled(static_cast<float>(ph), scalingFactory, -1.f);
            float _y = py * py;
            for (unsigned int pw = 0; pw < w; pw++) {
                // Does the point lie in the main cardoid
                float _x = pa.px[pw];
                float __x = _x - 0.25;
                float p = sqrt((__x * __x) + _y);
                float _p = p * p;
                if (_x <= p - (_p + _p) + 0.25) {
                    color = colorPalette[0];
                }
                else {
                    float x = 0, y = 0, xold = 0, yold = 0, period = 0, iter = 0;
                    unsigned int max_iter = color_MAX;
                    float x_ = 0, y_ = 0;
                    while (x_ + y_ <= cond1 && (unsigned)iter < max_iter) {
                        float xtemp = x_ - y_ + _x;
                        y = (2 * x * y) + py;
                        x = xtemp;
                        iter += 1;
                        x_ = x * x;
                        y_ = y * y;

                        if (x == xold && y == yold) {
                            iter = max_iter;
                            break;
                        }

                        period += 1.f;
                        if (period >= 20.f) {
                            period = 0;
                            xold = x;
                            yold = y;
                        }
                    }

                    if ((unsigned int)iter < max_iter) {
                        auto log_zn = std::log(x_ + y_) / 2;
                        auto nu = std::log(log_zn / log2) / log2;
                        iter += 1 - nu;
                        de_iter++;
                    }

                    double fpart;
                    color = linear_interpolation((float)colorPalette[(unsigned)iter],
                        (float)colorPalette[(unsigned)iter + 1],
                        (float)std::modf(iter, &fpart));
                }
                data[c] = color >> 0x10;
                data[c + 1] = color >> 0x8;
                data[c + 2] = color;
                c += 3;
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\t" << de_iter << " log calculations performed!\n";
    }

    bool equal_color(unsigned char* data, unsigned int c, unsigned int color) {
        if (data[c] == (unsigned char)(color >> 0x10)
            && data[c + 1] == (unsigned char)(color >> 0x08)
            && data[c + 2] == (unsigned char)color) return true;
        return false;
    }

    void bt_contour(unsigned char* data, int w, int h) {

        unsigned int lcolor = 0;
        int iter = 0;
        unsigned int c = 0;

        int off = 0;

        while (iter != 4) {

            for (unsigned int i = c; i < w * h * 3; i += 3) {
                if ((data[c] != (lcolor >> 0x10) || data[c + 1] != (lcolor >> 0x08) || data[c + 2] != lcolor)
                    && (data[c] != 0 || data[c + 1] != 0 || data[c + 2] != 0)) {
                    lcolor = data[c] << 0x10;
                    lcolor += data[c + 1] << 0x08;
                    lcolor += data[c + 2];
                    break;
                }
                c += 3;
            }

            if (lcolor == 0) return;

            int sp = c;
            int np = 0;
            bool end = false;

            for (int y = c / w; y < h && end != true;) {
                while (true) {
                    if (sp + (w * 3) + 3 >= w * h * 3) {
                        end = true;
                        break;
                    }
                    if (equal_color(data, sp + (w * 3), lcolor)) {
                        np = sp + (w * 3);
                        y++;
                        break;
                    }
                    if (equal_color(data, sp + (w * 3) + 3, lcolor)) {
                        np = sp + (w * 3) + 3;
                        y++;
                        break;
                    }
                    if (equal_color(data, sp + 3, lcolor)) {
                        np = sp + 3;
                        break;
                    }
                    else {
                        end = true;
                        break;
                    }
                }
                if (end == true) break;

                for (int x = 3; x < (sp + (w * 3)) % (w * 3); x += 3) {
                    if (data[sp - x] == 0x00 && data[sp - x + 1] == 0x00 && data[sp - x + 2] == 0x00) {
                        data[sp - x] = lcolor >> 0x10;
                        data[sp - x + 1] = lcolor >> 0x08;
                        data[sp - x + 2] = lcolor;
                    }
                    else {
                        break;
                    }
                }
                sp = np;
            }
            iter++;

            // Find the next neighboring pixels and increment c accordingly
            for (unsigned int i = c; i < w * h * 3; i += 3) {
                unsigned char r = lcolor >> 0x10;
                unsigned char g = lcolor >> 0x08;
                unsigned char b = lcolor;
                if (data[c] == r && data[c + 1] == g && data[c + 2] == b);
                else break;
                c += 3;
            }

        }
    }

    unsigned pixelcolor(unsigned maxiter, float px, float py) {
        int i = 0;
        float x = 0.f;
        float y = 0.f;
        float x2 = 0.f;
        float y2 = 0.f;
        while (x2 + y2 <= 4 && i < maxiter) {
            y = ((x + x) * y) + py;
            x = x2 - y2 + px;
            x2 = x * x;
            y2 = y * y;
            i++;
        }
        if (i >= maxiter) return 0x000000;
        return colorPalette[i];
    }
    unsigned pixelcolor(unsigned maxiter, float px, float py, unsigned& c_iter) {
        int i = 0;
        float x = 0.f;
        float y = 0.f;
        float x2 = 0.f;
        float y2 = 0.f;
        while (x2 + y2 <= 4 && i < maxiter) {
            y = ((x + x) * y) + py;
            x = x2 - y2 + px;
            x2 = x * x;
            y2 = y * y;
            i++;
        }
        c_iter = i;
        if (i >= maxiter) {
            c_iter = maxiter - 1;
            return 0x000000;
        }
        return colorPalette[i];
    }

    unsigned p8trace(unsigned& sx, unsigned& sy, unsigned char start, unsigned cc) {
        float kx = helpers::calculate_scaling_factor(-2.5f, 1.f, 0.f, static_cast<float>(pixels_x));
        float ky = helpers::calculate_scaling_factor(-1.f, 1.f, 0.f, static_cast<float>(pixels_y));
        unsigned kt = start + 1;
        char match_key = -1;
        char mismatch_key = -1;

        for (int i = 0; i < 8; i++) {
            if (kt == 1 && sy > 0 && sx > 0) {
                float py = helpers::calculate_scaled(static_cast<float>(sy - 1), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 1;
                }
                else {
                    match_key = 1;
                }
            }
            else if (kt == 2 && sy > 0) {
                float py = helpers::calculate_scaled(static_cast<float>(sy - 1), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx - 1), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 2;
                }
                else {
                    match_key = 2;
                }
            }
            else if (kt == 3 && sx > 0) {
                float py = helpers::calculate_scaled(static_cast<float>(sy), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx - 1), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 3;
                }
                else {
                    match_key = 3;
                }
            }
            else if (kt == 4 && sy < pixels_y && sx > 0) {
                float py = helpers::calculate_scaled(static_cast<float>(sy + 1), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx - 1), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 4;
                }
                else {
                    match_key = 4;
                }
            }
            else if (kt == 5 && sy < pixels_y) {
                float py = helpers::calculate_scaled(static_cast<float>(sy + 1), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 5;
                }
                else {
                    match_key = 5;
                }
            }
            else if (kt == 6 && sy < pixels_y && sx < pixels_x) {
                float py = helpers::calculate_scaled(static_cast<float>(sy + 1), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx + 1), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 6;
                }
                else {
                    match_key = 6;
                }
            }
            else if (kt == 7 && sx < pixels_x) {
                float py = helpers::calculate_scaled(static_cast<float>(sy), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx + 1), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 7;
                }
                else {
                    match_key = 7;
                }
            }
            else if (kt == 8 && sx < pixels_x && sy > 0) {
                float py = helpers::calculate_scaled(static_cast<float>(sy - 1), ky, -1.f);
                float px = helpers::calculate_scaled(static_cast<float>(sx + 1), kx, -2.5f);
                if (pixelcolor(color_MAX, px, py) != cc) {
                    mismatch_key = 8;
                }
                else {
                    match_key = 8;
                }
            }
            if (match_key != -1 && mismatch_key != -1) {
                break;
            }
            kt++;
            if (kt >= 8) kt = 1;
        }

        if (match_key == 1) {
            if (start == 1) return 0;
            sy--;
            return 5;
        }
        else if (match_key == 2) {
            if (start == 2) return 0;
            sx--;
            sy--;
            return 6;
        }
        else if (match_key == 3) {
            if (start == 3) return 0;
            sx--;
            return 7;
        }
        else if (match_key == 4) {
            if (start == 4) return 0;
            sx--;
            sy++;
            return 8;
        }
        else if (match_key == 5) {
            if (start == 5) return 0;
            sy++;
            return 1;
        }
        else if (match_key == 6) {
            if (start == 6) return 0;
            sx++;
            sy++;
            return 2;
        }
        else if (match_key == 7) {
            if (start == 7) return 0;
            sx++;
            return 3;
        }
        else if (match_key == 8) {
            if (start == 8) return 0;
            sx++;
            sy--;
            return 4;
        }

        return 0;
    }

    unsigned bytef_pixel_loc(unsigned sx, unsigned sy, unsigned key) {
        unsigned val = sy * pixels_x * 3;
        val += sx * 3;
        unsigned tmp = 0;
        if (key == 1) {
            if (sy - 1 == UINT32_MAX) return UINT32_MAX;
            return val - (pixels_x * 3);
        }
        else if (key == 2) {
            if (sy - 1 == UINT32_MAX) return UINT32_MAX;
            if (sx - 1 == UINT32_MAX) return UINT32_MAX;
            return val - (pixels_x * 3) - 3;
        }
        else if (key == 3) {
            if (sx - 1 == UINT32_MAX) return UINT32_MAX;
            return val - 3;
        }
        else if (key == 4) {
            if (sx - 1 == UINT32_MAX) return UINT32_MAX;
            if (sy + 1 >= pixels_y) return UINT32_MAX;
            return val + (pixels_x * 3) - 3;
        }
        else if (key == 5) {
            if (sy + 1 >= pixels_y) return UINT32_MAX;
            return val + (pixels_x * 3);
        }
        else if (key == 6) {
            if (sy + 1 >= pixels_y) return UINT32_MAX;
            if (sx + 1 >= pixels_x) return UINT32_MAX;
            return val + (pixels_x * 3) + 3;
        }
        else if (key == 7) {
            if (sx + 1 >= pixels_x) return UINT32_MAX;
            return val + 3;
        }
        else if (key == 8) {
            if (sx + 1 >= pixels_x) return UINT32_MAX;
            if (sy - 1 == UINT32_MAX) return UINT32_MAX;
            return val - (pixels_x * 3) + 3;
        }
    }

    unsigned colorBorder(unsigned char* data, unsigned sx, unsigned sy, unsigned color) {
        unsigned npx = sx;
        unsigned npy = sy;
        unsigned start = 0;
        unsigned xk = 0;
        unsigned c = 0;

        unsigned char r = color >> 0x10;
        unsigned char g = color >> 0x08;
        unsigned char b = color;

        unsigned iter = 0;

        while (true) {
            start = p8trace(npx, npy, start, color);
            if (start == 0) break;
            if (start <= 4) xk = start + 4;
            else xk = start - 4;

            c = bytef_pixel_loc(npx, npy, xk);
            if (c != UINT32_MAX) {
                if (data[c] == r && data[c + 1] == g && data[c + 2] == b) {
                    break;
                }
                else {
                    data[c] = r;
                    data[c + 1] = g;
                    data[c + 2] = b;
                }
            }
            else {
                break;
            }
            iter++;
        }
        return iter;
    }

    void bt_mandelbrot(unsigned char* data, int w, int h, unsigned sx, unsigned sy) {
        float kx = helpers::calculate_scaling_factor(-2.5f, 1.f, 0.f, static_cast<float>(w));
        float ky = helpers::calculate_scaling_factor(-1.f, 1.f, 0.f, static_cast<float>(h));

        unsigned ix = 0;
        unsigned nc = 0;

        unsigned xxc = 0;
        auto start = std::chrono::steady_clock::now();
        for (unsigned y = 0; y < h; y++) {
            float py = helpers::calculate_scaled(static_cast<float>(y), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(0), kx, -2.5f);
            unsigned cc = pixelcolor(color_MAX, px, py);
            for (unsigned x = 1; x < w; x++) {
                px = helpers::calculate_scaled(static_cast<float>(x), kx, -2.5f);
                nc = pixelcolor(color_MAX, px, py);
                if (nc != cc) {
                    unsigned scl = (y * w * 3) + (x * 3);
                    if ((data[scl] == 0 && data[scl + 1] == 0 && data[scl + 2] == 0)) {
                        colorBorder(data, x - 1, y, cc);
                        xxc++;
                    }
                    cc = nc;
                }
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

        std::cout << xxc << " borders traced!\n";

    }
}

namespace fr_experimental {

    int cx = 0;
    int cy = 0;

    int oox = 0;
    int ooy = 0;

    double dx = 0;
    double dy = 0;

    void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
    {
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
            oox = cx;
            ooy = cy;
        }
        else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
            dx = cx - oox;
            dy = cy - ooy;
        }
    }

    double zoomfactor = 0;

    void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
        zoomfactor = yoffset;
    }

    bool cardoid_check(float px, float py) {
        float __x = px - 0.25;
        float p = sqrt((__x * __x) + py);
        float _p = p * p;
        if (px <= p - (_p + _p) + 0.25) {
            return true;
        }
        return false;
    }

    unsigned compute_pixel_color_oet(float px, float py) {
        float x = 0, y = 0, xold = 0, yold = 0, period = 0;
        unsigned int iter = 0;
        unsigned int max_iter = color_MAX;
        float x_ = 0, y_ = 0;
        unsigned cond1 = 1 << 16;
        while (x_ + y_ <= cond1 && iter < max_iter) {
            float xtemp = x_ - y_ + px;
            y = (2 * x * y) + py;
            x = xtemp;
            iter += 1;
            x_ = x * x;
            y_ = y * y;

            if (x == xold && y == yold) {
                iter = max_iter;    // ???
                break;
            }

            period += 1.f;
            if (period >= 20.f) {
                period = 0;
                xold = x;
                yold = y;
            }
        }
        return iter;
    }

    unsigned compute_color(float px, float py) {       
        unsigned color;
        if (cardoid_check(px, py * py)) {
            color = 0x000000;
        }
        else {
            unsigned iter = compute_pixel_color_oet(px, py);
            color = colorPalette[iter];
        }
        return color;
    }
    
    class Mset {

        double lowx, lowy, highx, highy;

    public:

        Mset() : lowx(-2.5), lowy(-1.0), highx(1.0), highy(1.0) {}
        void oet_mandelbrot(unsigned char* data, int w, int h) {
            unsigned int c = 0;

            double scalingFactorx = helpers::calculate_scaling_factor(lowx, highx, 0.0, static_cast<double>(w));
            double scalingFactory = helpers::calculate_scaling_factor(lowy, highy, 0.0, static_cast<double>(h));

            if (dx != 0 && dy != 0) {
                double sdx = -dx / static_cast<double>(w);
                double sdy = dy / static_cast<double>(h);
                dx = 0;
                dy = 0;
                lowx += sdx;
                highx += sdx;
                lowy += sdy;
                highy += sdy;
            }

            if (zoomfactor > 0) {
                double zx = std::fabs(highx - lowx);
                double zy = std::fabs(highy - lowy);
                double tx = helpers::calculate_scaled(static_cast<double>(cx), scalingFactorx, lowx);
                double ty = helpers::calculate_scaled(static_cast<double>(cy), scalingFactory, lowy);
                lowx    = tx - (zx / 4);
                highx   = tx + (zx / 4);
                lowy    = ty - (zy / 4);
                highy   = ty + (zy / 4);
                zoomfactor = 0;
                printf("zx: %f\tzy: %f\n", zx, zy);
                printf("tx: %f\tty: %f\n", tx, ty);
                printf("cx: %d\tcy: %d\n", cx, cy);
                printf("(%f, %f) <-> (%f, %f)\n", lowx, lowy, highx, highy);
            }

            unsigned int color;
            for (unsigned int ph = 0; ph < h; ph++) {
                double py = helpers::calculate_scaled(static_cast<double>(ph), scalingFactory, lowy);
                for (unsigned int pw = 0; pw < w; pw++) {
                    double px = helpers::calculate_scaled(static_cast<double>(pw), scalingFactorx, lowx);
                    color = compute_color(px, py);
                    data[c] = color >> 0x10;
                    data[c + 1] = color >> 0x8;
                    data[c + 2] = color;
                    c += 3;
                }
            }
        }
        void oet_mandelbrot_il(unsigned char* data, int w, int h) {
            unsigned int c = 0;

            double scalingFactorx = helpers::calculate_scaling_factor(lowx, highx, 0.0, static_cast<double>(w));
            double scalingFactory = helpers::calculate_scaling_factor(lowy, highy, 0.0, static_cast<double>(h));

            if (dx != 0 && dy != 0) {
                double tx = helpers::calculate_scaled(static_cast<double>(-dx), scalingFactorx, 0.0);
                double ty = helpers::calculate_scaled(static_cast<double>(dy), scalingFactory, 0.0);
                dx = 0;
                dy = 0;
                printf("Positional change before dragging: L(%f, %f) H(%f, %f)\n", lowx, lowy, highx, highy);
                lowx += tx;
                highx += tx;
                lowy += ty;
                highy += ty;
                printf("After dragging: L(%f, %f) H(%f, %f)\n", lowx, lowy, highx, highy);
            }

            if (zoomfactor > 0) {
                double zx = std::fabs(highx - lowx);
                double zy = std::fabs(highy - lowy);
                double tx = helpers::calculate_scaled(static_cast<double>(cx), scalingFactorx, lowx);
                double ty = helpers::calculate_scaled(static_cast<double>(cy), scalingFactory, lowy);
                lowx = tx - (zx / 4);
                highx = tx + (zx / 4);
                lowy = ty - (zy / 4);
                highy = ty + (zy / 4);
                zoomfactor = 0;
                printf("zx: %f\tzy: %f\n", zx, zy);
                printf("tx: %f\tty: %f\n", tx, ty);
                printf("cx: %d\tcy: %d\n", cx, cy);
                printf("(%f, %f) <-> (%f, %f)\n", lowx, lowy, highx, highy);
            }

            unsigned int color, bc;
            auto start = std::chrono::steady_clock::now();
            unsigned stride = 60;
            unsigned maxstride = stride;
            unsigned size = stride / 3;

            unsigned adj_c = 0;

            for (unsigned int y = 0; y < h; y++) {
                double py = helpers::calculate_scaled(static_cast<double>(y), scalingFactory, lowy);
                for (unsigned int x = 0; x < w; x += size) {

                    unsigned decr = 0;

                    double px = helpers::calculate_scaled(static_cast<double>(x), scalingFactorx, lowx);
                    bc = compute_color(px, py);
                    data[c] = bc >> 0x10;
                    data[c + 1] = bc >> 0x8;
                    data[c + 2] = bc;
                    c += stride;

                    do {
                        px = helpers::calculate_scaled(static_cast<double>(x + size - decr), scalingFactorx, lowx);
                        color = compute_color(px, py);
                        data[c - (decr * 3)] = color >> 0x10;
                        data[c + 1 - (decr * 3)] = color >> 0x8;
                        data[c + 2 - (decr * 3)] = color;
                        decr++;
                    } while (color != bc && decr < size);

                    unsigned mismatch_points = decr - 1;

                    for (; decr < size; decr++) {
                        data[c - (decr * 3)] = bc >> 0x10;
                        data[c + 1 - (decr * 3)] = bc >> 0x8;
                        data[c + 2 - (decr * 3)] = bc;
                    }
                }
            }
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
            std::cout << adj_c << " times adjusted!\n";
        }
    };

    // This should primarily be used for zooming
    void oet_mandelbrot_il(unsigned char* data, int w, int h) {
        unsigned int c = 0;
        float scalingFactory = helpers::calculate_scaling_factor(0.0f, 1.0f, 0.0f, static_cast<float>(h));
        float scalingFactorx = helpers::calculate_scaling_factor(0.0f, 1.0f, 0.0f, static_cast<float>(w));
        unsigned int color, bc;
        auto start = std::chrono::steady_clock::now();
        unsigned stride = 60;
        unsigned maxstride = stride;
        unsigned size = stride / 3;

        unsigned adj_c = 0;

        for (unsigned int y = 0; y < h; y++) {
            float py = helpers::calculate_scaled(static_cast<float>(y), scalingFactory, 0.0f);
            for (unsigned int x = 0; x < w; x += size) {

                unsigned decr = 0;

                float px = helpers::calculate_scaled(static_cast<float>(x), scalingFactorx, 0.f);
                bc = compute_color(px, py);
                data[c] = bc >> 0x10;
                data[c + 1] = bc >> 0x8;
                data[c + 2] = bc;
                c += stride;

                do {
                    px = helpers::calculate_scaled(static_cast<float>(x + size - decr), scalingFactorx, -2.5f);
                    color = compute_color(px, py);
                    data[c - (decr * 3)] = color >> 0x10;
                    data[c + 1 - (decr * 3)] = color >> 0x8;
                    data[c + 2 - (decr * 3)] = color;
                    decr++;
                } while (color != bc && decr < size);

                unsigned mismatch_points = decr - 1;

                for (; decr < size; decr++) {
                    data[c - (decr * 3)] = bc >> 0x10;
                    data[c + 1 - (decr * 3)] = bc >> 0x8;
                    data[c + 2 - (decr * 3)] = bc;
                }
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
        std::cout << adj_c << " times adjusted!\n";
    }

}   // namespace fr_experimental

void generateRandomContour(unsigned char *data, int w, int h, int px, int py, int color) {

    unsigned int c = 0;
    int i = 0;

    while (true) {

        c = (py * (w * 3)) + (px * 3);
        data[c] = color >> 0x10;
        data[c + 1] = color >> 0x08;
        data[c + 2] = color;

        int cell = std::rand() % 3;
        if(cell == 0) {
            py++;
        } else if (cell == 1) {
            py++; px++;
        } else if (cell == 2) {
            px++;
        }

        if (px >= w || py >= h) {
            break;
        }
    }
}

unsigned int nonblackpixels(unsigned char* data, unsigned int size) {
    unsigned int c = 0;
    for (unsigned int i = 0; i < size; i += 3) {
        unsigned char r = data[i] >> 0x10;
        unsigned char g = data[i + 1] >> 0x08;
        unsigned char b = data[i + 2];
        if (r != 0 || g != 0 || b != 0) c++;
    }
    return c;
}

unsigned int ddxyo(unsigned char* data, unsigned int radius, unsigned int x, unsigned int y, unsigned int ox) {
    unsigned int count = radius * 2;
    for (unsigned int py = 0; py < y; py++) {
        for (unsigned int px = 0; px < ox; px++) {
           
        }
    }
    return count;
}

void generateCircle(unsigned char* data, unsigned int color, unsigned int pradius) { // pradius = radius in pixels
    double pi = 3.14159265358979323846;
    unsigned int ox = pixels_x / 2;
    unsigned int oy = pixels_y / 2;

    float gradient = (2 * pi / (float)pradius);
    float t = 0;

    for (unsigned int px = 0; px < pradius * 200; px++) {
        double x = ox + ((float)pradius * std::cos(t));
        double y = oy + ((float)pradius * std::sin(t));

        unsigned int c = (std::lround(y) * pixels_x * 3) + (std::lround(x) * 3);
        if (c + 3 <= pixels_x * pixels_y * 3) {
            data[c] = color >> 0x10;
            data[c + 1] = color >> 0x08;
            data[c + 2] = color;
        }

        t += gradient;
    }
}

void generatepl(unsigned char* data, unsigned int size) {

    for (unsigned int y = 200; y < 700; y++) {
        data[(y * pixels_x * 3) + (700 * 3)] = 0xFF;
        data[(y * pixels_x * 3) + ((700 * 3) + 1)] = 0xFF;
        data[(y * pixels_x * 3) + ((700 * 3) + 2)] = 0xFF;
        data[(y * pixels_x * 3) + (1100 * 3)] = 0xFF;
        data[(y * pixels_x * 3) + ((1100 * 3) + 1)] = 0xFF;
        data[(y * pixels_x * 3) + ((1100 * 3) + 2)] = 0xFF;
    }

}

void print_pxo(pixel_array& pa) {
    for (unsigned int x = 0; x < pixels_x; x++) {
        printf("%d: %f\n", x, pa.px[x]);
    }
}

void dividersOfX(unsigned x) {
    for (unsigned int i = 2; i < x / 2; i++) {
        if (x % i == 0) printf("%d, ", i);
    }
}

unsigned colorAtXY(unsigned x, unsigned y) {
    float kx = helpers::calculate_scaling_factor(-2.5f, 1.f, 0.f, static_cast<float>(pixels_x));
    float ky = helpers::calculate_scaling_factor(-1.f, 1.f, 0.f, static_cast<float>(pixels_y));
    float px = helpers::calculate_scaled(static_cast<float>(x), kx, -2.5f);
    float py = helpers::calculate_scaled(static_cast<float>(y), ky, -1.f);
    return pixelcolor(color_MAX, px, py);
}

unsigned dataAtXY(unsigned char *data, unsigned x, unsigned y) {
    unsigned xy = (x * 3) + (y * pixels_x * 3);
    unsigned color = data[xy] << 0x10;
    color += data[xy + 1] << 0x08;
    color += data[xy + 2];
    return color;
}

void pixelsWithColor(unsigned char* data, unsigned x, unsigned y, unsigned color) {
    unsigned matches = 0;
    for (unsigned py = 0; py < y; py++) {
        for (unsigned px = 0; px < x; px++) {
            unsigned c = dataAtXY(data, px, py);
            if (c == color) {
                std::cout << "(" << px << ", " << py << ")\n";
                matches++;
            }
        }
    }
    std::cout << std::dec << matches << " pixels found matching color 0x" << std::hex << color << '\n';
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    fr_experimental::cx = xpos;
    fr_experimental::cy = ypos;
}

int main(int argc, char **argv) {

    if (argc == 3) {
        // Double-precision floating point arithmetic
    }

    std::srand(std::time(nullptr));

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    if (!glfwInit()) {
        glfwTerminate();
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(pixels_x, pixels_y, "Mandelbrot Renderer", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, GLFW_TRUE);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetScrollCallback(window, fr_experimental::scroll_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, fr_experimental::mouse_button_callback);

    float vertices[] = {
         1.f,  1.f, 0.f,    1.f, 0.f, 0.f,  1.f, 1.f,
         1.f, -1.f, 0.f,    0.f, 1.f, 0.f,  1.0f, 0.0f,
        -1.f, -1.f, 0.f,    0.f, 0.f, 1.f,  0.0f, 0.0f,
        -1.f,  1.f, 0.f,    1.f, 1.f, 0.f,  0.0f, 1.f
    };

    unsigned int indices[] = {
        0, 1, 3, // first triangle
        1, 2, 3  // second triangle
    };

    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO); // we can also generate multiple VAOs or buffers at the same time
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    helpers::ShaderReader sr("C:\\generic\\shader.vs");
    const char* vs_src = sr.read();
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vs_src, NULL);
    glCompileShader(vertexShader);

    sr.update("C:\\generic\\fragment.fs");
    const char* fs_src = sr.read();
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fs_src, NULL);
    glCompileShader(fragmentShader);

    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    unsigned int mandelTexture;
    glGenTextures(1, &mandelTexture);
    glBindTexture(GL_TEXTURE_2D, mandelTexture);

    glUseProgram(shaderProgram);
    glUniform1i(glGetUniformLocation(shaderProgram, "mandelTexture"), 0);

    unsigned char* pixels = (unsigned char*)malloc(pixels_x * pixels_y * 3);
    memset(pixels, 0, pixels_x * pixels_y * 3);

    fillcolorp_grey();
    fr_experimental::Mset mset;
    mset.oet_mandelbrot_il(pixels, pixels_x, pixels_y);

    if (argc == 1) {

        while (!glfwWindowShouldClose(window)) {

            processInput(window);

            glClearColor(0.5f, 0.5f, 0.5f, 1.f);
            glClear(GL_COLOR_BUFFER_BIT);
 
            mset.oet_mandelbrot_il(pixels, pixels_x, pixels_y);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, mandelTexture);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pixels_x, pixels_y, 0, GL_RGB, GL_UNSIGNED_BYTE, (void*)pixels);
            glGenerateMipmap(GL_TEXTURE_2D);

            glUseProgram(shaderProgram);
            glBindVertexArray(VAO);
            glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

            glfwSwapBuffers(window);
            glfwPollEvents();
        }
    }

    free(pixels);

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);

}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}