# Native kernel tests

Run the standalone smoke test under AddressSanitizer and UndefinedBehaviorSanitizer
from the repository root with Clang:

```sh
clang++ -std=c++17 -O1 -g -fno-fast-math -ffp-contract=off \
  -fsanitize=address,undefined -fno-omit-frame-pointer -Icpp/include \
  cpp/fastcall.cpp cpp/tests/fastcall_smoke.cpp -o /tmp/excavator2-fastcall-sanitized
/tmp/excavator2-fastcall-sanitized
```

The test exercises E/M/posterior buffers at 1, 17, and 10,003 segments. Numerical
parity, boundary validation and concurrent calls are checked through pytest in
`tests/test_fastcall.py` and `tests/test_fastcall_native.py`, also on installed wheels.

The same command with `hslm` replacing `fastcall` runs the HSLM smoke test.
It verifies traced and rolling-buffer paths agree at the same three sizes.
Both sanitizer checks are included in the Python port CI workflow.
