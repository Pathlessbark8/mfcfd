#!/usr/bin/env python
# Copyright 2020 Birla Institute of Technology and Science - Pilani, Hyderabad Campus
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import argparse, json, os, platform, shutil, subprocess, sys

_version = sys.version_info.major

if _version == 2: # Python 2.x:
    _input = raw_input
elif _version == 3: # Python 3.x:
    _input = input
else:
    raise Exception("Incompatible Python version")

mfcfd_build_dir = os.path.join(os.getcwd(), "build")
mfcfd_parent_dir = os.path.dirname(os.getcwd())
mfcfd_working_dir = os.getcwd()

def cleanup_build_directory():
    if os.path.isdir("build"):
        try:
            print("Cleaning build directory")
            shutil.rmtree("build")
            os.mkdir("build")
        except:
            raise Exception("Unable to clean build directory")
            sys.exit(1)
    else:
        os.mkdir("build")

def remove_pre_existing_executable():
    if os.path.isfile("execname"):
        try:
            os.remove("execname")
        except:
            raise Exception("Unable to remove pre existing executables")
            sys.exit(1)

def build(mfcfd_type=None, extra_flags=[], dest_path=None):

    if mfcfd_type == "serial":
        mfcfd_var = "SERIAL"
    elif mfcfd_type == "tangent":
        mfcfd_var = "TANGENT"
    elif mfcfd_type == "adjoint":
        mfcfd_var = "ADJOINT"
    elif mfcfd_type == "cuda":
        mfcfd_var = "CUDA"
    else:
        mfcfd_var = "SERIAL"
        
    try:
        subprocess.check_call(["cmake", mfcfd_parent_dir, "-DMFCFD={}".format(mfcfd_var)], cwd=mfcfd_build_dir)
    except subprocess.CalledProcessError:
        print("CMAKE failed")
        sys.exit(1)
    try:
        subprocess.check_call(["make"] + extra_flags, cwd=mfcfd_build_dir)
    except:
        print("MAKE failed")
        sys.exit(1)
    if dest_path == None:
        shutil.move(os.path.join(mfcfd_build_dir, "execname"), mfcfd_working_dir)
    else:
        if os.path.isdir(dest_path):
            shutil.move(os.path.join(mfcfd_build_dir, "execname"), dest_path)
        else:
            print("Invalid Destination Path")
            sys.exit(1)
    print("MFCFD has been built and a 'execname' file has been generated in the current working directory.")

def install(mfcfd_type=None, extra_flags=[], dest_path=None):
    print("Building {} Meshfree Solver.".format(mfcfd_type))
    cleanup_build_directory()
    remove_pre_existing_executable()

    if os.path.isfile(os.path.join(mfcfd_parent_dir, "CMakeLists.txt")):
        fline=open(os.path.join(mfcfd_parent_dir, "CMakeLists.txt")).readline().rstrip()
        if fline == "#MFCFD CMAKE SCRIPT":
            print("Verified cmake script.")
            try:
                subprocess.check_output(["cmake", "--version"])
            except OSError:
                print("Error: CMake is not installed or otherwise not executable. Please check")
                print("your CMake installation and try again.")
                print("Attempted to execute: {}".format("cmake"))
                sys.exit(1)
            build(mfcfd_type, extra_flags, dest_path)
        else:
            raise Exception("Invalid CMakeLists.txt found.")
            sys.exit(1)
    else:
        raise Exception("Unable to find CMakeLists.txt. Are you running the script from the correct directory?")
        sys.exit(1)
    

def driver():
    parser = argparse.ArgumentParser(
        description="Compile and Install Meshfree Solver.")
    parser.add_argument(
        "--mfcfd", dest="mfcfd_type", required=True,
        choices=["cuda", "serial", "tangent", "adjoint"], default="serial",
        help="Meshfree Solver to install.")
    parser.add_argument(
        "--extra", dest="extra_flags", action="append", required=False,
        default=[],
        help="Extra flags for make command.")
    parser.add_argument(
        "--destination", dest="dest_path", required=False,
        default=None,
        help="Destination Path for saving execname.")
    args = parser.parse_args()

    install(**vars(args))

if __name__ == "__main__":
    driver()