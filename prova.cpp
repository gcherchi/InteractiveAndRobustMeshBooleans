//
// Created by Michele on 13/02/25.
//

#include "booleans.h"
#include "cinolib/meshes/trimesh.h"
#include "cinolib/meshes/meshes.h"
#include <iostream>
#include <filesystem>
#include <unordered_set>
#include <string>
#include <fstream>

namespace fs = std::filesystem;

// Function to execute the command, capture output, and write to log
void executeAndLog(const std::string& command, const std::string& modelName, std::ofstream& logFile) {
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error executing command" << std::endl;
        return;
    }

    std::string manifold, watertight, localOrient, globalOrient, intersection;
    char buffer[128];

    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::cout << buffer;  // Print output to console in real-time

        std::string line(buffer);
        if (line.find("Manifold check:") != std::string::npos)
            manifold = (line.find("passed") != std::string::npos) ? "yes" : "no";
        else if (line.find("Watertight check:") != std::string::npos)
            watertight = (line.find("passed") != std::string::npos) ? "yes" : "no";
        else if (line.find("Local  Orientation check:") != std::string::npos)
            localOrient = (line.find("passed") != std::string::npos) ? "yes" : "no";
        else if (line.find("Global Orientation check:") != std::string::npos)
            globalOrient = (line.find("passed") != std::string::npos) ? "yes" : "no";
        else if (line.find("Intersection check:") != std::string::npos)
            intersection = (line.find("passed") != std::string::npos) ? "yes" : "no";
    }

    pclose(pipe);

    // Print parsed results before writing to file
    std::cout << "Parsed Output: " << modelName << " " << manifold << " " << watertight << " "
              << localOrient << " " << globalOrient << " " << intersection << std::endl;

    // Write results to statistics file
    logFile << modelName << " " << manifold << " " << watertight << " "
            << localOrient << " " << globalOrient << " " << intersection << std::endl;
}

int main(int argc, char **argv) {
    std::string path_folder_test = "../folder_test";
    std::string path_folder_origin = "/Tinghi10K";
    std::string name_folder_rotated = "/mesh_rotated";
    std::string name_folder_output = "/mesh_output";

    fs::path folderPath = path_folder_test + path_folder_origin;
    fs::path folderPathRotated = path_folder_test + name_folder_rotated;
    fs::path folderPathOutput = path_folder_test + name_folder_output;
    fs::path logFilePath = path_folder_test + "/statistics.txt";

    std::ofstream logFile(logFilePath, std::ios::app);
    if (!logFile.is_open()) {
        std::cerr << "Error opening statistics file!" << std::endl;
        return 1;
    }

    std::unordered_set<std::string> filesInFolder2;
    for (const auto &entry : fs::directory_iterator(folderPathRotated)) {
        if (fs::is_regular_file(entry.path())) {
            filesInFolder2.insert(entry.path().filename().string());
        }
    }

    for (const auto &entry: fs::directory_iterator(folderPath)) {
        if (fs::is_regular_file(entry.path())) {
            std::string fileName = entry.path().filename().string();
            fs::path fileOutput = folderPathOutput / fileName;

            if (filesInFolder2.find(fileName) != filesInFolder2.end()) {
                fs::path fileRotated = folderPathRotated / fileName;
                std::cout << "Processing: " << fileName << " with " << fileRotated << std::endl;

                std::string operation = "union";
                fs::path name_bool_output = folderPathOutput / fileName;
                const char* exe = "../cmake-build-debug/mesh_booleans";

                // Construct the command WITHOUT redirection (popen captures output directly)
                std::string command = std::string(exe) + " " + operation + " "
                                      + entry.path().string() + " "
                                      + fileRotated.string() + " "
                                      + name_bool_output.string();

                // Execute command and log output
                executeAndLog(command, name_bool_output.filename().string(), logFile);
            } else {
                std::cout << "   -> No corresponding file found in " << folderPathRotated << std::endl;
            }
        }
    }

    logFile.close();
    return 0;
}
