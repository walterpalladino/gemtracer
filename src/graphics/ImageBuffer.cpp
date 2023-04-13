
#include "graphics/ImageBuffer.h"

void ImageBuffer::Create(int width, int height)
{
    this->width = width;
    this->height = height;
    data.resize(width * height * 4, 0);
}

void ImageBuffer::Clear()
{
    data.resize(width * height * 4, 0);
}

void ImageBuffer::Plot(long x, long y, unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
    data[4 * width * y + 4 * x + 0] = r;
    data[4 * width * y + 4 * x + 1] = g;
    data[4 * width * y + 4 * x + 2] = b;
    data[4 * width * y + 4 * x + 3] = a;
}
