
#include <vector>

class ImageBuffer
{
public:
    int width;
    int height;
    std::vector<unsigned char> data;

    void Create(int width, int height);
    void Clear();
    void Plot(long x, long y, unsigned char r, unsigned char g, unsigned char b, unsigned char a);
};
