#include "booleans.h"
#include "cinolib/meshes/trimesh.h"
#include "cinolib/meshes/meshes.h"
#include <iostream>
#include <filesystem>
#include <unordered_set>
#include <string>
#include <fstream>
#include <vector>
#include <mutex>
#include <queue>
#include <thread>
#include <condition_variable>
#include <functional>
#include <map>
#include <sstream>

namespace fs = std::filesystem;
bool rotation_enabled = false;
bool relative = false;
bool parallel = true;



std::mutex file_mutex; // mutex globale per sincronizzare accesso al file

// ================== Thread Pool ====================
class ThreadPool {
public:
    ThreadPool(size_t threads);
    ~ThreadPool();
    void enqueue(std::function<void()> task);

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

ThreadPool::ThreadPool(size_t threads) : stop(false) {
    for(size_t i = 0; i < threads; ++i)
        workers.emplace_back([this] {
            while(true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                    if(this->stop && this->tasks.empty()) return;
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }
                task();
            }
        });
}

void ThreadPool::enqueue(std::function<void()> task) {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        tasks.emplace(std::move(task));
    }
    condition.notify_one();
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers) worker.join();
}

// ============= Execute and Log Function ==============
void executeAndLog(const std::string& command, const std::string& modelName,
                   const std::string& logFilePath, std::ofstream& exceptionLog) {
    FILE* pipe = popen((command + " 2>&1").c_str(), "r");
    if (!pipe) {
        std::cerr << "Error executing command" << std::endl;
        exceptionLog << modelName << " Error executing command" << std::endl;
        return;
    }

    std::string manifold = "no", watertight = "no" , localOrient = "no", globalOrient = "no", intersection = "no";
    std::string errors;
    char buffer[256];

    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::string line(buffer);
        std::cout << buffer;

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

        if (line.find("Assertion failed") != std::string::npos || line.find("Error") != std::string::npos) {
            errors += line + " | ";
        }
    }
    pclose(pipe);

    std::ostringstream newEntryStream;
    newEntryStream << modelName << " " << manifold << " " << watertight << " "
                   << localOrient << " " << globalOrient << " " << intersection;
    std::string newEntry = newEntryStream.str();

    // Scrittura sicura con mutex
    {
        std::lock_guard<std::mutex> lock(file_mutex);

        std::map<std::string, std::string> log_entries;
        std::ifstream inFile(logFilePath);
        std::string line;
        while (std::getline(inFile, line)) {
            std::istringstream iss(line);
            std::string name;
            iss >> name;
            log_entries[name] = line;
        }
        inFile.close();

        log_entries[modelName] = newEntry;

        std::ofstream outFile(logFilePath, std::ios::trunc);
        for (const auto& [key, value] : log_entries) {
            outFile << value << "\n";
        }
    }

    if (!errors.empty()) {
        exceptionLog << modelName << " " << errors << std::endl;
    }
}

// ====================== MAIN =========================
int main(int argc, char **argv) {
    //fs::path script_dir = fs::absolute(fs::path(argv[0])).parent_path();
    fs::path script_dir = "/home/michele/Documents/GitHub/InteractiveAndRobustMeshBooleans/cmake-build-release";
    std::string path_folder_test = relative ? script_dir.string() : "../folder_test";
    std::string path_folder_origin = "/Tinghi10K";
    std::string name_folder_rotated = "/mesh_rotated";
    std::string name_folder_output = "/mesh_bool_output";

    std::vector<std::string> operations = {"union"};

    fs::path folderPath = path_folder_test + path_folder_origin;
    fs::path folderPathRotated = path_folder_test + name_folder_rotated;
    fs::path folderPathOutput = path_folder_test + name_folder_output;

    const char* meshNames[] = {
        "67608_sf_a.obj",
        "103824_sf_a.obj",
        "1582379_sf_a.obj",
        "59756_sf_a.obj",
        "42194_sf_a.obj",
        "153968_sf_a.obj",
        "472003_sf_a.obj",
        "627211_sf_c.obj",
        "103355_sf_a.obj",
        "471984_sf_a.obj",
        "472105_sf_a.obj",
        "205449_sf_a.obj",
        "68509_sf_a.obj",
        "103825_sf_a.obj",
        "49333_sf_a.obj",
        "1356644_sf_a.obj",
        "472071_sf_a.obj",
        "1458693_sf_a.obj",
        "135216_sf_a.obj",
        "118037_sf_a.obj",
        "51808_sf_a.obj",
        "348502_sf_a.obj",
        "565652_sf_a.obj",
        "74150_sf_a.obj",
        "472161_sf_c.obj",
        "112421_sf_a.obj",
        "87430_sf_a.obj",
        "129919_sf_a.obj",
        "1387346_sf_a.obj",
        "47568_sf_a.obj",
        "1063851_sf_a.obj",
        "1505034_sf_a.obj",
        "1095207_sf_a.obj",
        "137977_sf_a.obj",
        "87431_sf_a.obj",
        "187284_sf_a.obj",
        "1458694_sf_a.obj",
        "471995_sf_a.obj",
        "42652_sf_a.obj",
        "59771_sf_a.obj",
        "1322464_sf_a.obj",
        "1275117_sf_a.obj",
        "65778_sf_a.obj",
        "77938_sf_a.obj",
        "472136_sf_a.obj",
        "472098_sf_a.obj",
        "42192_sf_a.obj",
        "49888_sf_a.obj",
        "472052_sf_a.obj",
        "87429_sf_a.obj",
        "111013_sf_a.obj",
        "39929_sf_a.obj",
        "59764_sf_a.obj",
        "1384321_sf_a.obj",
        "257599_sf_a.obj",
        "80942_sf_a.obj",
        "129895_sf_a.obj",
        "87432_sf_a.obj",
        "104559_sf_c.obj",
        "84781_sf_a.obj",
        "113890_sf_a.obj",
        "118302_sf_a.obj",
        "43393_sf_a.obj",
        "1582423_sf_a.obj",
        ".DS_Store",
        "1582435_sf_a.obj",
        "314000_sf_a.obj",
        "69057_sf_a.obj",
        "1439537_sf_a.obj",
        "37323_sf_a.obj",
        "109925_sf_a.obj",
        "43769_sf_a.obj",
        "49423_sf_a.obj",
        "471994_sf_a.obj",
        "111006_sf_c.obj",
        "51511_sf_a.obj",
        "472008_sf_a.obj",
        "326895_sf_a.obj"
    };


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

    ThreadPool pool(24);
    bool startProcessing = false;
    std::string startFrom = "280281_sf_a.obj";

    /*for (const auto &entry: fs::directory_iterator(folderPath)) {
        if (fs::is_regular_file(entry.path())) {
            std::string fileName = entry.path().filename().string();

            if (!startProcessing) {
                if (fileName == startFrom) startProcessing = true;
                else continue;
            }

            if (filesInFolder2.find(fileName) != filesInFolder2.end()) {
                fs::path fileRotated = folderPathRotated / fileName;

                for (const auto& operation : operations) {
                    pool.enqueue([entry, fileRotated, operation, script_dir, path_folder_test]() {
                        fs::path name_bool_output = fs::path(path_folder_test) / "mesh_bool_output" / operation / entry.path().filename();
                        fs::path logFilePath = fs::path(path_folder_test) / ("statistics_" + operation + ".txt");
                        fs::path exceptionFilePath = fs::path(path_folder_test) / ("exceptions_" + operation + ".txt");
                        fs::path exePath = fs::absolute(script_dir / "mesh_booleans");

                        std::ofstream exceptionLog(exceptionFilePath, std::ios::app);

                        if (!exceptionLog.is_open()) {
                            std::cerr << "Error opening logs for " << entry.path() << std::endl;
                            return;
                        }

                        std::string command = exePath.string() + " " + operation + " "
                                              + entry.path().string() + " "
                                              + fileRotated.string() + " "
                                              + name_bool_output.string();

                        executeAndLog(command, name_bool_output.filename().string(), logFilePath.string(), exceptionLog);
                    });
                }
            } else {
                std::cout << "   -> No corresponding file found in " << folderPathRotated << std::endl;
            }
        }
    }*/
    for (const std::string& fileName : meshNames) {


        fs::path entryPath = folderPath / fileName;

        if (!fs::exists(entryPath)) {
            std::cerr << "File not found: " << entryPath << std::endl;
            continue;
        }

        if (filesInFolder2.find(fileName) != filesInFolder2.end()) {
            fs::path fileRotated = folderPathRotated / fileName;

            for (const auto& operation : operations) {
                auto task = [entryPath, fileRotated, operation, script_dir, path_folder_test]() {
                    fs::path name_bool_output = fs::path(path_folder_test) / "mesh_bool_output" / operation / entryPath.filename();
                    fs::path logFilePath = fs::path(path_folder_test) / ("statistics_" + operation + ".txt");
                    fs::path exceptionFilePath = fs::path(path_folder_test) / ("exceptions_" + operation + ".txt");
                    fs::path exePath = fs::absolute(script_dir / "mesh_booleans");
                    std::ofstream exceptionLog(exceptionFilePath, std::ios::app);

                    std::string command = exePath.string() + " " + operation + " "
                                          + entryPath.string() + " "
                                          + fileRotated.string() + " "
                                          + name_bool_output.string();

                    std::cout << "Processing :" << entryPath.filename()  <<  std::endl;
                    executeAndLog(command, name_bool_output.filename().string(), logFilePath, exceptionLog);
                };

                if (parallel) {
                    pool.enqueue(task);
                } else {
                    task();
                }
            }

        } else {
            std::cout << "   -> No corresponding rotated file found for " << fileName << std::endl;
        }
    }


    return 0;
}
