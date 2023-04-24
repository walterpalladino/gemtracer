#include <iostream>

#include "core/renderer/View.h"
#include "core/scene/SceneLoader.h"

#include "lodepng.h"

// using namespace std;

int main(int argc, char *argv[])
{
    const char *filename = argc > 1 ? argv[1] : "test.png";

    SceneLoader::GetInstance()->Load("./resources/scenes/scene_1.txt");

    //    View::GetInstance()->Init();
    View::GetInstance()->Init(1280, 720);

    View::GetInstance()->SetCamera(Vector3d(0.0, 100.0, -200.0),
                                   Vector3d(0.0, 0.0, 1.0),
                                   Vector3d(1280.0 / 720.0, 0.0, 0.0),
                                   Vector3d(0.0, 1.0, 0.0),
                                   90.0);

    //  Render the scene
    View::GetInstance()->RenderScene(false);

    //  Save image as PNG
    std::vector<unsigned char> image = View::GetInstance()->GetImageBuffer().data;

    unsigned error = lodepng::encode(filename,
                                     image,
                                     View::GetInstance()->GetImageBuffer().width,
                                     View::GetInstance()->GetImageBuffer().height);

    return 0;
}
