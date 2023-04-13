#include <iostream>
// #include "App.h"
// #include "lodepng.h"

#include "core/renderer/View.h"

#include "lodepng.h"

// using namespace std;

int main(int argc, char *argv[])
{
    const char *filename = argc > 1 ? argv[1] : "test.png";

    View::GetInstance()->Init();

    float fAngle = 0;
    float fXPos = 128 << 2;
    float fYPos = 128 << 2;
    float fZPos = 50;
    float fStep = 32;
    float fDelta = 15; // 3.14159 / 12.0 ;

    View::GetInstance()->RenderScene(
        fAngle,
        fZPos,
        fXPos,
        fYPos);
    /*
        encodeOneStep(filename,
                      View::GetInstance()->GetImageBuffer().data,
                      View::GetInstance()->GetImageBuffer().width,
                      View::GetInstance()->GetImageBuffer().height);
    */
    std::vector<unsigned char> image = View::GetInstance()->GetImageBuffer().data;

    unsigned error = lodepng::encode(filename,
                                     image,
                                     View::GetInstance()->GetImageBuffer().width,
                                     View::GetInstance()->GetImageBuffer().height);

    return 0;
}
