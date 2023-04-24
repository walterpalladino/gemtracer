#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <any>

#include "core/scene/SceneLoader.h"

using namespace std;

void SceneLoader::Load(string fileName)
{
    //    std::unordered_map<std::string, std::any> map = {};
    std::unordered_map<std::string, std::string> map;

    ifstream configFile;
    configFile.open(fileName);

    if (!configFile.is_open())
    {
        //  TODO : Rise an exception
        cout << "File not found" << endl;
        return;
    }

    string line;
    size_t lineNbr = 1;

    while (getline(configFile, line))
    {
        //  Remove spaces
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

        //  Skip if empty line
        if (line.empty())
        {
            lineNbr++;
            continue;
        }

        //  Skip comments
        if (line.at(0) == '#')
        {
            lineNbr++;
            continue;
        }

        vector<string> tokens;
        stringstream ss(line);
        while (ss.good())
        {
            string tmp;
            getline(ss, tmp, '=');
            tokens.push_back(tmp);
        }
        //  TODO : validate amount of elements
        map[tokens[0]] = tokens[1];

        cout << "Line " << lineNbr << " contains " << tokens[0] << " with value " << tokens[1] << endl;

        lineNbr++;
    }

    configFile.close();
}
