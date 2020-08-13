#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include <glad/glad.h>
#include <glfw3.h>

#include "stb_image.h"

int sHeight = 1080;
int sWidth = 1920;

const int pixels_x = 1920;
const int pixels_y = 1080;

float mouseWheeloffset = 0;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        //medianx += 0.01f;
        printf("D pressed!\n");
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        //medianx -= 0.01f;
        printf("A pressed!\n");
    }
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        //mediany += 0.01f;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        //mediany -= 0.01f;
    }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
        //mx1 -= 0.01f;
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        //my2 -= 0.01f;
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

}

unsigned int colorPalette[256];

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
    while (i < 256) {
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
    while (i < 256) {
        colorPalette[i] = c;
        i++;
        c += 0x010101;
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
    const float mx2 = 1;        // Mandelbrot X scale (-2.5, 1)
    const float mx1 = -2.5;
    const float my2 = 1;        // Mandelbrot Y scale (-1, 1)
    const float my1 = -1;
    float scalingFactorx = helpers::calculate_scaling_factor(mx1, mx2, 0.f, static_cast<float>(w));
    float scalingFactory = helpers::calculate_scaling_factor(my1, my2, 0.f, static_cast<float>(h));

    unsigned int c = 0;
    auto start = std::chrono::steady_clock::now();
    for (unsigned int ph = 0; ph < pixels_y; ph++) {
        float py = helpers::calculate_scaled(static_cast<float>(ph), scalingFactory, my1);
        for (unsigned int pw = 0; pw < pixels_x; pw++) {
            float px = helpers::calculate_scaled(static_cast<float>(pw), scalingFactorx, mx1);
            float x = 0.f;
            float y = 0.f;
            float x2 = 0.f;
            float y2 = 0.f;
            int i = 0;
            int max_i = 256;
            while (x2 + y2 <= 4 && i < max_i) {
                y = ((x + x) * y) + py;
                x = x2 - y2 + px;
                x2 = x * x;
                y2 = y * y;
                i++;
            }   
            int color = colorPalette[i];
            data[c]     =   (color >> 0x10) & 0xFF;
            data[c + 1] =   (color >> 0x8) & 0xFF;
            data[c + 2] =   color & 0xFF;
            c += 3;
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}


template<typename T>
T linear_interpolation(T a, T b, T c) {
    return a + ((1/3) * (b - a - c));
}

#define SIZE_PX (1920)

struct pixel_array {
    float px[SIZE_PX];
};

void calculate_pxo(unsigned int w, pixel_array& pa) {
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

void cc_mandelbrot_rect(unsigned char* data, int w, int h, pixel_array& pa) {
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

void cc_mandelbrot(unsigned char* data, int w, int h, pixel_array& pa) {
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
            } else {
                float x = 0, y = 0, xold = 0, yold = 0, period = 0, iter = 0;
                unsigned int max_iter = 256;
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
            data[c]     = color >> 0x10;
            data[c + 1] = color >> 0x8;
            data[c + 2] = color;
            c += 3;
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\t" << de_iter << " log calculations performed!\n";
}

// (x * x) + (y * y) <= (1 << 16) && iter < max_iter

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


void bt_circle(unsigned char* data, int w, int h) {

}
/*
void bt_mandelbrot(unsigned char *data, int w, int h) {
    int max_i = 256;
    for(unsigned int y = 0; y < h; y++) {
        float py = helpers::calculate_scaled(static_cast<float>(ph), scalingFactory, my1);
        for(unsigned int x = 0; x < w; x++) {
            float px = helpers::calculate_scaled(static_cast<float>(pw), scalingFactorx, mx1);
            int i = 0;
            float x = 0.f;
            float y = 0.f;
            float x2 = 0.f;
            float y2 = 0.f;
            while (x2 + y2 <= 4 && i < max_i) {
                y = ((x + x) * y) + py;
                x = x2 - y2 + px;
                x2 = x * x;
                y2 = y * y;
                i++;
            }
            if(i < max_i) break;
        }
    }
}
*/

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

unsigned p8trace(unsigned& sx, unsigned& sy, unsigned char start, unsigned cc) {
    float kx = helpers::calculate_scaling_factor(-2.5f, 1.f, 0.f, static_cast<float>(pixels_x));
    float ky = helpers::calculate_scaling_factor(-1.f, 1.f, 0.f, static_cast<float>(pixels_y));
    unsigned kt = start;
    char match_key = -1;
    char mismatch_key = -1;

    unsigned b8color[8];

    for (int i = 0; i < 8; i++) {
        if (kt == 1 && sy > 0 && sx > 0) {
            float py = helpers::calculate_scaled(static_cast<float>(sy - 1), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 1;
            } else {
                match_key = 1;
            }
        }
        else if (kt == 2 && sy > 0) {
            float py = helpers::calculate_scaled(static_cast<float>(sy - 1), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx - 1), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 2;
            } else {
                match_key = 2;
            }
        }
        else if (kt == 3 && sx > 0) {
            float py = helpers::calculate_scaled(static_cast<float>(sy), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx - 1), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 3;
            } else {
                match_key = 3;
            }
        }
        else if (kt == 4 && sy < pixels_y && sx > 0) {
            float py = helpers::calculate_scaled(static_cast<float>(sy + 1), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx - 1), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 4;
            } else {
                match_key = 4;
            }
        }
        else if (kt == 5 && sy < pixels_y) {
            float py = helpers::calculate_scaled(static_cast<float>(sy + 1), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 5;
            } else {
                match_key = 5;
            }
        }
        else if (kt == 6 && sy < pixels_y && sx < pixels_x) {
            float py = helpers::calculate_scaled(static_cast<float>(sy + 1), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx + 1), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 6;
            } else {
                match_key = 6;
            }
        }
        else if (kt == 7 && sx < pixels_x) {
            float py = helpers::calculate_scaled(static_cast<float>(sy), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx + 1), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 7;
            } else {
                match_key = 7;
            }
        }
        else if (kt == 8 && sx < pixels_x && sy > 0) {
            float py = helpers::calculate_scaled(static_cast<float>(sy - 1), ky, -1.f);
            float px = helpers::calculate_scaled(static_cast<float>(sx + 1), kx, -2.5f);
            if (pixelcolor(256, px, py) != cc) {
                mismatch_key = 8;
            } else {
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
        sy--;
        return 5;
    }
    else if (match_key == 2) {
        sx--;
        sy--;
        return 6;
    }
    else if (match_key == 3) {
        sx--;
        return 7;
    }
    else if (match_key == 4) {
        sx--;
        sy++;
        return 8;
    }
    else if (match_key == 5) {
        sy++;
        return 1;
    }
    else if (match_key == 6) {
        sx++;
        sy++;
        return 2;
    }
    else if (match_key == 7) {
        sx++;
        return 3;
    }
    else if (match_key == 8) {
        sx++;
        sy--;
        return 4;
    }

    return 0;
}

unsigned bytef_pixel_loc(unsigned sx, unsigned sy, unsigned key) {
    auto val = sy * pixels_x * 3;
    val += sx * 3;
    if (key == 1) {
        return val - (pixels_x * 3);
    } else if (key == 2) {
        return val - (pixels_x * 3) - 3;
    } else if (key == 3) {
        return val - 3;
    } else if (key == 4) {
        return val + (pixels_x * 3) - 3;
    } else if (key == 5) {
        return val + (pixels_x * 3);
    } else if (key == 6) {
        return val + (pixels_x * 3) + 3;
    } else if (key == 7) {
        return val + 3;
    } else if (key == 8) {
        return val - (pixels_x * 3) + 3;
    }
}

void bt_mandelbrot(unsigned char* data, int w, int h) {
    float kx = helpers::calculate_scaling_factor(-2.5f, 1.f, 0.f, static_cast<float>(w));
    float ky = helpers::calculate_scaling_factor(-1.f, 1.f, 0.f, static_cast<float>(h));
    unsigned cc = pixelcolor(256, -1.f, -2.5f);
    unsigned ctable[50];
    unsigned ctable_iter = 0;
    unsigned spx = 0;
    unsigned spy = 0;
    bool switchAlgoToBT = false;
    for (unsigned y = 0; y < h; y++) {
        float py = helpers::calculate_scaled(static_cast<float>(y), ky, -1.f);
        for (unsigned x = 1; x < w; x++) {
            float px = helpers::calculate_scaled(static_cast<float>(x), kx, -2.5f);
            float nc = pixelcolor(256, px, py);
            if (nc != cc) {
                std::cout << "nc: " << nc << '\n';
                switchAlgoToBT = true;
                ctable[ctable_iter] = cc;
                ctable_iter++;
                spx = x;
                spy = y;
                break;
            }
        }
        if (switchAlgoToBT) break;
    }

    std::cout << "Color: " << cc << '\n';
    // Border-tracing
    unsigned npx = spx;
    unsigned npy = spy;
    bool cont = false;

    unsigned start = 0;
    unsigned xk = 0;
    unsigned c = 0;

    for (unsigned int i = 0; i < pixels_y - 1; i++) {
        start = p8trace(npx, npy, start + 1, cc);
        if (start <= 4) xk = start + 4;
        else xk = start - 4;
        
        c = bytef_pixel_loc(npx, npy, xk);
        data[c] = 0xFF;
        data[c + 1] = 0x00;
        data[c + 2] = 0x00;

        std::cout << start << '\n';
    }
}

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

int main(int argc, char **argv) {

    //calculate_pxo(pixels_x);
    std::srand(std::time(nullptr));
    pixel_array pa;
    calculate_pxo(pixels_x, pa);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    if (!glfwInit()) {
        glfwTerminate();
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(sWidth, sHeight, "Mandelbrot Renderer", NULL, NULL);
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

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetScrollCallback(window, scroll_callback);

    const char* img = "C:\\Users\\abkour\\source\\repos\\TextRenderer\\img\\wall.jpg";
    const char* img2 = "C:\\Users\\abkour\\source\\repos\\TextRenderer\\img\\awesomeface.png";

    float vertices[] = {
        // positions          // colors           // texture coords
         1.f,  1.f, 0.f,   1.f, 0.f, 0.f,   1.f, 1.f, // top right
         1.f, -1.f, 0.f,   0.f, 1.f, 0.f,   1.0f, 0.0f, // bottom right
        -1.f, -1.f, 0.f,   0.f, 0.f, 1.f,   0.0f, 0.0f, // bottom left
        -1.f,  1.f, 0.f,   1.f, 1.f, 0.f,   0.0f, 1.f  // top left 
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

    helpers::ShaderReader sr("C:\\Users\\abkour\\source\\repos\\TextRenderer\\shader.vs");
    const char* vs_src = sr.read();
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vs_src, NULL);
    glCompileShader(vertexShader);

    sr.update("C:\\Users\\abkour\\source\\repos\\TextRenderer\\fragment.fs");
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

    generatePaletteTexture(colorPalette, 256);

    unsigned int mandelTexture;
    glGenTextures(1, &mandelTexture);
    glBindTexture(GL_TEXTURE_2D, mandelTexture);

    glUseProgram(shaderProgram);
    glUniform1i(glGetUniformLocation(shaderProgram, "mandelTexture"), 0);

    unsigned char* pixels = (unsigned char*)malloc(pixels_x * pixels_y * 3);
    memset(pixels, 0, pixels_x * pixels_y * 3);
    fillcolorp3();
    //generateRandomContour(pixels, pixels_x, pixels_y, 0, 0, 0xabcdef);
    //generateRandomContour(pixels, pixels_x, pixels_y, 200, 0, 0x123456);
    //generateRandomContour(pixels, pixels_x, pixels_y, 400, 0, 0x789abc);
    //generateRandomContour(pixels, pixels_x, pixels_y, 600, 0, 0x00FFFF);
    //generateCircle(pixels, 0xFFFFFF, 700);
    //generateCircle(pixels, 0xFF00FF, 300);
    //generatepl(pixels, 2);
    //auto count = nonblackpixels(pixels, pixels_x * pixels_y * 3);
    //std::cout << "Non-black pixels: " << count << '\n';
    //generateCircle(pixels, 0xFFFFFF, 250);
    //generateCircle(pixels, 0x00FFFF, 200);
    //bt_circle(pixels, pixels_x, pixels_y);
    //mandelbrotSet(pixels, pixels_x, pixels_y);
    bt_mandelbrot(pixels, pixels_x, pixels_y);
    //bt_contour(pixels, pixels_x, pixels_y);
    //colorsquares(pixels, pixels_x, pixels_y);
    //printContentsOfCP();

    if (argc == 1) {

        while (!glfwWindowShouldClose(window)) {

            processInput(window);

            glClearColor(0.5f, 0.5f, 0.5f, 1.f);
            glClear(GL_COLOR_BUFFER_BIT);

            //cc_mandelbrot(pixels, pixels_x, pixels_y, pa);
            //colorsquares(pixels, pixels_x, pixels_y);     
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

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    mouseWheeloffset = yoffset;
}