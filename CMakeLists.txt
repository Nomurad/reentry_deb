# CMake version 指定
cmake_minimum_required(VERSION 2.8)

# プロジェクト名を指定
project(test_cmake CXX)

add_subdirectory(src)

# # a.exe を作成
# add_executable(a 
#     src/main.cpp
#     )