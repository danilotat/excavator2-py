"""Literal MakeReadCount/EXOMECOUNT state machine for debugging small streams."""

import numpy as np

from excavator2.reads import CHUNK_SIZE


def count_positions(chunks, starts, ends):
    counts = np.zeros(len(starts), dtype=np.int64)
    window, residual = 0, 0
    for chunk in chunks:
        if not len(chunk):
            raise ValueError("legacy MakeReadCount cannot process an empty chunk")
        if chunk[-1] <= ends[window]:
            residual = len(chunk)
        else:
            count = residual
            for position in chunk:
                while (
                    window < len(starts) and position >= starts[window] and position > ends[window]
                ):
                    counts[window] = count
                    count = 0
                    window += 1
                if window == len(starts):
                    return counts
                if starts[window] <= position <= ends[window]:
                    count += 1
            residual = count
        if len(chunk) < CHUNK_SIZE:
            break
    return counts
