# Instructions

## Compile and run
### `python` code
`python burgers.py`

### `c++` code

Compile use `g++ -o burgers burgers.cpp`.
Run with `./burgers`. (Authorize by `chmod +x burgers`)

### `julia` wrapper
Compile the `c++` code with the extra `-fPIC -shared` flags.
Also make the `-o libburgers.so` to follow best practices.
Run `julia main.jl`

### `CUDA` code

Compile using 
`nvcc -arch=sm_50 burgers.cu -o test`.
Run with `./test`. 

## Performance
![Performance](https://github.com/SichengHe/cuda_burgers/blob/main/results.png)