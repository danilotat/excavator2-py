# Focused numerical core

`fastcall.cpp` implements the scalar E-step, fixed-mean M-step, and posterior
kernels. Their readable specifications are in `src/excavator2/reference/fastcall.py`.
The binding validates arrays and releases the GIL during computation. Parameters,
iteration/stopping, cellularity, assignment, file access, and output stay in Python.

See `docs/fastcall-reference.md` for shapes, compatibility behavior and evidence.
HSLM is not implemented. Additional kernels and SIMD remain later milestone work.
