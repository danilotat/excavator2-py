# Focused numerical core

Only the extension binding exists today. Add `hslm.cpp` for transitions,
emissions, and Viterbi; add `fastcall.cpp` for batched E/M and posterior kernels.
Keep parameters, algorithm control, filtering, file access, and output in Python.

Every kernel needs documented array shapes, a readable implementation in
`src/excavator2/reference/`, and comparison with legacy checkpoints. Additional
C++ requires measured performance or numerical compatibility justification.
Do not add an I/O framework or SIMD backend before there is evidence for it.
