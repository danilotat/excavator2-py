# Benchmarks

After installing the package, run from the repository root:

```sh
python benchmarks/fastcall.py --segments 10003 --repeats 5
```

This deterministic synthetic benchmark checks posterior and iteration agreement
before timing both backends. JSON output records platform, Python, NumPy, input
size, iteration count, median full-fit times and speedup. It includes Python
control but excludes label assignment, HSLM and I/O; it does not measure pipeline
speedup or peak memory. Each fit runs scalar kernels on the calling thread.
Record compiler/build details separately when comparing environments.
