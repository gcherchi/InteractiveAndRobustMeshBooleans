#include "booleans.h"
#include "cinolib/meshes/trimesh.h"
#include "cinolib/meshes/meshes.h"
#include <iostream>
#include <filesystem>
#include <unordered_set>
#include <string>
#include <fstream>

namespace fs = std::filesystem;
bool rotation_enabled = false;

// Function to execute the command, capture output, and log errors
void executeAndLog(const std::string& command, const std::string& modelName, std::ofstream& logFile, std::ofstream& exceptionLog) {
    FILE* pipe = popen((command + " 2>&1").c_str(), "r"); // Capture both stdout and stderr
    if (!pipe) {
        std::cerr << "Error executing command" << std::endl;
        exceptionLog << modelName << " Error executing command" << std::endl;
        return;
    }

    std::string manifold, watertight, localOrient, globalOrient, intersection;
    std::string errors;
    char buffer[256];

    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::string line(buffer);
        std::cout << buffer;  // Print output live

        // Extract results
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

        // Capture errors or assertion failures
        if (line.find("Assertion failed") != std::string::npos || line.find("Error") != std::string::npos) {
            errors += line + " | ";  // Store errors in one line
        }
    }
    pclose(pipe);

    // Print parsed results before writing to file
    std::cout << "Parsed Output: " << modelName << " " << manifold << " " << watertight << " "
              << localOrient << " " << globalOrient << " " << intersection << std::endl;

    // Write results to statistics file
    logFile << modelName << " " << manifold << " " << watertight << " "
            << localOrient << " " << globalOrient << " " << intersection << std::endl;

    // Log errors if any
    if (!errors.empty()) {
        exceptionLog << modelName << " " << errors << std::endl;
    }
}

int main(int argc, char **argv) {

    fs::path script_dir = fs::absolute(fs::path(argv[0])).parent_path();
    std::string path_folder_test = script_dir.string();
    std::string path_folder_origin = "/Tinghi10K";
    std::string name_folder_rotated = "/mesh_rotated";
    std::string name_folder_output = "/mesh_bool_output";

    std::vector<std::string> operations = {"union", "intersection", "subtraction"};

    fs::path folderPath = path_folder_test + path_folder_origin;
    fs::path folderPathRotated = path_folder_test + name_folder_rotated;
    fs::path folderPathOutput = path_folder_test + name_folder_output;


    if(rotation_enabled) {
        for (const auto& mesh : fs::directory_iterator(folderPath)) {
            if (fs::is_regular_file(mesh.path())) {
                cinolib::Trimesh<> m(mesh.path().c_str());
                m.rotate(cinolib::vec3d(1, 0, 0), 90);
                fs::path rotatedMesh = folderPathRotated / mesh.path().filename();
                m.save(rotatedMesh.c_str());
            }
        }
    }

    for (const auto& operation : operations) {
        fs::create_directories(folderPathOutput / operation);
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

            if (filesInFolder2.find(fileName) != filesInFolder2.end()) {
                fs::path fileRotated = folderPathRotated / fileName;
                std::cout << "Processing: " << fileName << " with " << fileRotated << std::endl;

                for (const auto& operation : operations) {
                    fs::path name_bool_output = folderPathOutput / operation / fileName;
                    fs::path logFilePath = path_folder_test + "/statistics_" + operation + ".txt";
                    fs::path exceptionFilePath = path_folder_test + "/exceptions_" + operation + ".txt";
                    fs::path exePath = fs::absolute(script_dir / "mesh_booleans");
                    const std::string exe = exePath.string();


                    std::ofstream logFile(logFilePath, std::ios::app);
                    if (!logFile.is_open()) {
                        std::cerr << "Error opening statistics file for " << operation << "!" << std::endl;
                        continue;
                    }

                    std::ofstream exceptionLog(exceptionFilePath, std::ios::app);
                    if (!exceptionLog.is_open()) {
                        std::cerr << "Error opening exception log file for " << operation << "!" << std::endl;
                        continue;
                    }

                    std::string command = std::string(exe) + " " + operation + " "
                                          + entry.path().string() + " "
                                          + fileRotated.string() + " "
                                          + name_bool_output.string();

                    executeAndLog(command, name_bool_output.filename().string(), logFile, exceptionLog);
                }
            } else {
                std::cout << "   -> No corresponding file found in " << folderPathRotated << std::endl;
            }
        }
    }
    return 0;
}
