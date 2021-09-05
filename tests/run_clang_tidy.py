import os
import glob
import subprocess

if len(glob.glob("compile_commands.json")) == 0:
    print("compile_commands.json doesn't exists")

files = glob.glob("../src/*.cpp")

for counter,file in enumerate(files):
    print("File #"+str(counter+1)+"/"+str(len(files))+": "+file + " --> " + file.split("/")[-1].split(".")[0] +"_clang_tidy.txt")
    os.system("clang-tidy -checks=\"*\" "+file+" > " + file.split("/")[-1].split(".")[0]+"_clang_tidy.txt 2> /dev/null")

