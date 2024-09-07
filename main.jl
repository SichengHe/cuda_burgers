# --- Julia 1.10---
"""
@File          :   main.jl
@Date created  :   2024/09/07
@Last modified :   2024/09/07
@Author        :   Galen Ng
@Desc          :   Calling C/C++ functions from Julia
                   From: https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/
"""

# @ccall is used to call C functions from Julia
@ccall "./libburgers.so".main()::Int32
