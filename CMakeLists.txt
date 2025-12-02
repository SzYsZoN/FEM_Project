cmake_minimum_required(VERSION 3.10)
project(mes)

set(CMAKE_CXX_STANDARD 17)

add_executable(mes
    main.cpp
    GlobalData.cpp
    Grid.cpp
    Gauss.cpp
    ElemUniv.cpp
    Jacobian.cpp
    
)


configure_file(Test1.txt Test1.txt COPYONLY) 
configure_file(Test2.txt Test2.txt COPYONLY)
configure_file(Test3.txt Test3.txt COPYONLY)
