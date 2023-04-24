#pragma once

class SceneLoader
{
public:
    static SceneLoader *GetInstance()
    {
        static SceneLoader instance;
        return &instance;
    };

    void Load(std::string fileName);
};