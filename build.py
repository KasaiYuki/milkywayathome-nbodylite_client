import sys
import subprocess
import os as ospkg

# Simple script to start builds using jenkins CI
# No dependecy checking, etc

# Usage
def usage():
    print("Usage: python build.py (win,linux,mac) (32|64) (ON|OFF)")

# Simple wrapper function
def execute(args):
    proc = subprocess.Popen(args, shell=False)
    proc.communicate()
    if proc.returncode != 0:
        print("ERROR: \"" + " ".join(args) + "\" exited with a non-zero return code")
        exit(1)

# Check correct number of args
if len(sys.argv) != 4:
    usage()
    exit(1)

# Make sure os and architecture are valid
os  = sys.argv[1]
arch = sys.argv[2]
sep_opencl = sys.argv[3]
assert os in ["win", "linux", "mac"], "ERROR: Unknown OS " + os
assert arch in ["32", "64"], "ERROR: Unknown arch " + bit 
assert sep_opencl in ["ON", "OFF"], "ERROR: Set SEP_OPENCL to ON or OFF" #No longer using multithreaded

# CMake flags used for all platforms
cmake_shared_flags = ["-DSEPARATION=ON", "-DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER", "-DNBODY=" + sep_opencl, "-DSEPARATION_OPENCL=" + sep_opencl]

# CMake flags used for windows
cmake_static_flag = ["-DNBODY_STATIC=ON",  "-DBOINC_APPLICATION=ON",  "-DCMAKE_BUILD_TYPE=Release"]

cmake_windows_toolchain = ["-DCMAKE_TOOLCHAIN_FILE=cmake_modules/MinGW32-Cross-Toolchain.cmake"]

try:
    execute(["rm", "-r", "build"])
except:
    print "No build file found, continuing"
execute(["mkdir", "build"])
ospkg.chdir("./build")


# Linux
if os == "linux":

    if arch == "64":
        execute(["cmake", "..", "-DOPENCL_LIBRARIES=/usr/lib/x86_64-linux-gnu/libOpenCL.so", "-DOPENCL_INCLUDE_DIRS=/srv/chroot/hardy_amd64/usr/local/cuda/include/", "-DBUILD_32=OFF"] + cmake_shared_flags + cmake_static_flag)
    if arch == "32":
        print("Warning: Linux 32 bit not supported for Milky Way N-body")
        execute(["cmake", "..", "-DBUILD_32=ON", "-DNBODY=OFF", "-DSEPARATION_STATIC=ON", "-DCMAKE_FIND_ROOT_PATH=/srv/chroot/hardy_i386", "-DOPENCL_LIBRARIES=/srv/chroot/hardy_i386/usr/lib/libOpenCL.so", "-DOPENCL_INCLUDE_DIRS=/srv/chroot/hardy_i386/usr/local/cuda/include/"] + cmake_shared_flags + cmake_static_flag)
    
    execute(["make", "clean"])
    execute(["make"])

# Windows
if os == "win":

    if arch == "64":
        execute(["cmake", "..", "-DOPENCL_INCLUDE_DIRS=/home/jenkins/cuda/include", "-DOPENCL_LIBRARIES=/home/jenkins/x86_64/OpenCL.lib"] + cmake_shared_flags + cmake_static_flag + cmake_windows_toolchain)

    if arch == "32":
        #print("ERROR: Windows 32 bit not supported")
        execute(["cmake", "..", "-DBUILD_32=ON", "-DOPENCL_INCLUDE_DIRS=/home/jenkins/cuda/include", "-DOPENCL_LIBRARIES=/home/jenkins/x86/OpenCL.lib"] + cmake_shared_flags + cmake_static_flag + cmake_windows_toolchain)
    
    execute(["make", "clean"])
    execute(["make"])

# Mac OS X
if os == "mac":

    if arch == "64":
        execute(["cmake", ".."] + cmake_shared_flags + cmake_static_flag + ["-DCMAKE_C_FLAGS='-O3'" ])

    if arch == "32":
        print("ERROR: Mac OSX 32 bit not supported")

    execute(["/usr/bin/make", "clean"])
    execute(["/usr/bin/make"])

