import smtplib
from email.message import EmailMessage
import os
import subprocess
import time
import re


passed = 0 
failed = 0 

def completion_bar(iterations, total):
    # Calculate the completion ratio
    ratio = iterations / total

    # Determine the length of the progress bar
    bar_length = 50
    progress = int(bar_length * ratio)

    # Create the completion bar
    bar = "[" + "#" * progress + "-" * (bar_length - progress) + "]"

    return bar

def compile_and_run_cpp(mesh_file, mode="release"):
    global passed,  failed
    executable = os.path.join("cmake-build-"+ mode, 'mesh_booleans_inputcheck')
    
    
    run_process = subprocess.run([executable, mesh_file], capture_output=True, text=True)

    if run_process.returncode == 0:
        #save the mesh file in mesh_ready folder 
        #mesh_ready_folder = "./mesh_ready_to_test"
        #os.system("cp " + mesh_file + " " + mesh_ready_folder)

        #remove the mesh file 
        
        passed += 1
    else: 
        os.system("rm " + mesh_file)
        failed += 1



if __name__ == "__main__":
        
    os.system('clear')

    print("Starting the test")
    #path_folder = "/Users/michele/Desktop/VMCG/Thingi10K"
    path_folder = "/Users/michele/Documents/GitHub/InteractiveAndRobustMeshBooleans/mesh_ready_to_test/mesh_rotated"

    results = []
    files = os.listdir(path_folder)
    total_files = len(files)
    for i, file in enumerate(files):
        if file.endswith(".off"):
            progress_bar = completion_bar(i, total_files-1)
            print("\rProgress: {}% {}".format(int(((i+1) / total_files) * 100), progress_bar), "Files Processed -> ", i+1,"/", total_files," Running on file: ", file, end="")
            compile_and_run_cpp(os.path.join(path_folder, file), "release")
    
    
    print ("\nPassed: ", passed, " Failed: ", failed)


    