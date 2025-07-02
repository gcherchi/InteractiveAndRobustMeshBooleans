#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
#include <thread>
#include <mutex>
#include <functional>
#include <filesystem>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <map>

namespace fs = std::filesystem;

std::mutex file_mutex;

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
    for (size_t i = 0; i < threads; ++i)
        workers.emplace_back([this] {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock, [this] {
                        return this->stop || !this->tasks.empty();
                    });
                    if (this->stop && this->tasks.empty()) return;
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
    for (std::thread &worker : workers)
        worker.join();
}

struct BooleanResult {
    std::string filename;
    std::string operation;
    std::string manifold = "no";
    std::string watertight = "no";
    std::string local_orient = "no";
    std::string global_orient = "no";
    std::string intersection = "no";
    int poly_input = -1;
    int poly_arrangement = -1;
    double pipeline_time = -1.0;
};

BooleanResult executeAndParse(const std::string &command, const std::string &filename) {
    BooleanResult result;
    result.filename = filename;

    FILE *pipe = popen((command + " 2>&1").c_str(), "r");
    if (!pipe) {
        std::cerr << "Failed to run command\n";
        return result;
    }

    char buffer[512];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::string line(buffer);
        std::cout << buffer;


        if (line.find("Operation :") != std::string::npos) {
            std::istringstream iss(line);
            std::string tmp;
            iss >> tmp >> tmp >> result.operation;
        } else if (line.find("Boolean Pipeline time") != std::string::npos) {
            auto start = line.find('[') + 1;
            auto end = line.find('s');
            result.pipeline_time = std::stod(line.substr(start, end - start));
        } else if (line.find("Dimension input:") != std::string::npos) {
            result.poly_input = std::stoi(line.substr(line.find_last_of('-') + 2));
        } else if (line.find("Dimension arrangement:") != std::string::npos) {
            result.poly_arrangement = std::stoi(line.substr(line.find_last_of('-') + 2));
        } else if (line.find("Manifold check") != std::string::npos) {
            result.manifold = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Watertight check") != std::string::npos) {
            result.watertight = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Local  Orientation check") != std::string::npos) {
            result.local_orient = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Global Orientation check") != std::string::npos) {
            result.global_orient = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Intersection check") != std::string::npos) {
            result.intersection = line.find("passed") != std::string::npos ? "yes" : "no";
        }
    }

    pclose(pipe);
    return result;
}

void writeToLog(const std::string &logFile, const BooleanResult &res, std::mutex &file_mutex) {
    std::lock_guard<std::mutex> lock(file_mutex);
    std::ofstream out(logFile, std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Could not open log file: " << logFile << std::endl;
        return;
    }
    out << res.filename << " "
        << res.manifold << " "
        << res.watertight << " "
        << res.local_orient << " "
        << res.global_orient << " "
        << res.intersection << " "
        << res.poly_input << " "
        << res.poly_arrangement << " "
        << res.pipeline_time << "\n";
}

// Parse multipli risultati da output che contiene tutte e 3 le operazioni
std::vector<BooleanResult> executeAndParseMultiple(const std::string &command, const std::string &filename) {
    std::vector<BooleanResult> results;
    BooleanResult currentResult;
    currentResult.filename = filename;

    int poly_input = -1;
    int poly_arrangement = -1;
    double pipeline_time = -1.0;

    FILE *pipe = popen((command + " 2>&1").c_str(), "r");
    if (!pipe) {
        std::cerr << "Failed to run command\n";
        return results;
    }

    char buffer[512];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::string line(buffer);
        std::cout << buffer;

        // Estraggo valori globali che appaiono solo una volta
        if (line.find("Boolean Pipeline time") != std::string::npos) {
            auto start = line.find('[') + 1;
            auto end = line.find('s');
            pipeline_time = std::stod(line.substr(start, end - start));
        } else if (line.find("Dimension input:") != std::string::npos) {
            poly_input = std::stoi(line.substr(line.find_last_of('-') + 2));
        } else if (line.find("Dimension arrangement:") != std::string::npos) {
            poly_arrangement = std::stoi(line.substr(line.find_last_of('-') + 2));
        }

        // Quando inizio un nuovo risultato
        if (line.find("Operation :") != std::string::npos) {
            if (!currentResult.operation.empty()) {
                // Assegno i valori globali a ogni risultato
                currentResult.poly_input = poly_input;
                currentResult.poly_arrangement = poly_arrangement;
                currentResult.pipeline_time = pipeline_time;

                results.push_back(currentResult);
                currentResult = BooleanResult();
                currentResult.filename = filename;
            }
            std::istringstream iss(line);
            std::string tmp;
            iss >> tmp >> tmp >> currentResult.operation;
        } else if (line.find("Manifold check") != std::string::npos) {
            currentResult.manifold = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Watertight check") != std::string::npos) {
            currentResult.watertight = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Local  Orientation check") != std::string::npos) {
            currentResult.local_orient = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Global Orientation check") != std::string::npos) {
            currentResult.global_orient = line.find("passed") != std::string::npos ? "yes" : "no";
        } else if (line.find("Intersection check") != std::string::npos) {
            currentResult.intersection = line.find("passed") != std::string::npos ? "yes" : "no";
        }
    }

    if (!currentResult.operation.empty()) {
        currentResult.poly_input = poly_input;
        currentResult.poly_arrangement = poly_arrangement;
        currentResult.pipeline_time = pipeline_time;
        results.push_back(currentResult);
    }

    pclose(pipe);
    return results;
}


void printHelp(const char* progName) {
    std::cout << "Usage: " << progName << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help                Show this help message and exit\n"
              << "  --serial              Run tests in serial mode (no parallelism)\n"
              << "  --operation [op]      Select operation to perform:\n"
              << "                        union          Perform union operation only\n"
              << "                        intersection   Perform intersection operation only\n"
              << "                        subtraction    Perform subtraction operation only\n"
              << "                        all            Perform all operations in a single run\n"
              << "\n"
              << "If --operation is not specified, default is 'all'.\n"
              << "Example:\n"
              << "  " << progName << " --operation union\n"
              << "  " << progName << " --serial --operation all\n";
}


int main(int argc, char **argv) {
    bool parallel = true;
    std::string selected_op = "all";
    std::string multiple_test = "none";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            printHelp(argv[0]);
            return 0;
            }else if (arg == "--serial") {
            parallel = false;
        } else if (arg == "--operation" && i + 1 < argc) {
            selected_op = argv[++i];
            if(selected_op == "all") multiple_test = "all";

            if (selected_op != "union" && selected_op != "intersection" && selected_op != "subtraction" && selected_op != "all") {
                std::cerr << "Invalid operation: " << selected_op << std::endl;
                return 1;
            }
        }

    }

    std::vector<std::string> all_operations = {"union", "intersection", "subtraction"};

    std::vector<std::string> operations;
    if (selected_op == "all") {
        operations = {"all"}; // singolo comando che produce tutte
    } else {
        operations = {selected_op}; // singolo comando per operazione
    }

    std::string base_path = "../folder_test";
    std::string input_folder = base_path + "/Tinghi10K";
    std::string rotated_folder = base_path + "/mesh_rotated";
    std::string output_folder = base_path + "/mesh_bool_output";
    std::string exe_path = fs::absolute(fs::path(argv[0])).parent_path().string() + "/mesh_booleans";

    std::vector<std::string> meshFiles;
    for (const auto &entry : fs::directory_iterator(input_folder)) {
        if (fs::is_regular_file(entry.path()))
            meshFiles.push_back(entry.path().filename().string());
    }

    ThreadPool pool(std::thread::hardware_concurrency());

    for (const auto &meshName : meshFiles) {
        std::string input_mesh = input_folder + "/" + meshName;
        std::string rotated_mesh = rotated_folder + "/" + meshName;

        if (!fs::exists(rotated_mesh)) continue;

        for (const auto &op : operations) {
            auto task = [=]() {
                std::string output_mesh = output_folder + "/";
                if (op == "all") {
                    output_mesh += "operation/" + meshName; // eventuale cartella per all
                } else {
                    output_mesh += op + "/" + meshName;
                }

                std::string command = exe_path + " ";
                if (op == "all") {
                    command += "all ";
                } else {
                    command += op + " ";
                }
                command += input_mesh + " " + rotated_mesh + " " + output_mesh ;


                std::cout << "Running test on: " << input_mesh << " - " << rotated_mesh << " WITH: "<<  selected_op << std::endl;

                if (op == "all") {
                    auto results = executeAndParseMultiple(command, meshName);
                    for (const auto &res : results) {
                        std::string logFile = base_path + "/statistics_" + res.operation + ".txt";
                        std::cout << "Writing to log: " << logFile << std::endl;
                        writeToLog(logFile, res, file_mutex);
                    }
                } else {
                    BooleanResult result = executeAndParse(command, meshName);
                    std::string logFile = base_path + "/statistics_" + op + ".txt";
                    std::cout << "Writing to log: " << logFile << std::endl;
                    writeToLog(logFile, result, file_mutex);
                }
            };

            if (parallel) {
                pool.enqueue(task);
            } else {
                task();
            }
        }
    }

    return 0;
}