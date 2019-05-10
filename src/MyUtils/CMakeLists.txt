# CMake version 指定
cmake_minimum_required(VERSION 2.8)

# 静的ライブラリを作成
add_library(MyUtils STATIC 
            MyUtils.hpp
            MyUtils.cpp)