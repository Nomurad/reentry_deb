# CMake version 指定
cmake_minimum_required(VERSION 2.8)

# module dir を追加
add_subdirectory(MyUtils)
add_subdirectory(ReentryTrj)

add_executable(a main.cpp)

# 追加したライブラリとリンク
target_link_libraries(a MyUtils)
target_link_libraries(a ReentryTrj)