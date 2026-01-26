"""
Pytest configuration and fixtures for EXCAVATOR2 tests.
"""

import pytest
import sys
from pathlib import Path

# Add src to path for testing
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


@pytest.fixture
def sample_log2_ratios():
    """Generate sample log2-ratio data for testing."""
    import numpy as np
    np.random.seed(42)
    # Simulate some simple CNV segments
    data = np.concatenate([
        np.random.normal(0.0, 0.1, 100),    # Normal region
        np.random.normal(-0.5, 0.1, 50),    # Deletion
        np.random.normal(0.0, 0.1, 100),    # Normal region
        np.random.normal(0.58, 0.1, 50),    # Duplication
        np.random.normal(0.0, 0.1, 100),    # Normal region
    ])
    return data


@pytest.fixture
def sample_positions():
    """Generate sample genomic positions."""
    import numpy as np
    # Positions every 1000bp
    return np.arange(1000, 401000, 1000, dtype=np.int64)
