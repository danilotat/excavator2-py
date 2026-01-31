"""
Basic tests to verify the package structure and C++ module.
"""

import pytest


def test_import_excavator2():
    """Test that the main package can be imported."""
    import excavator2

    assert excavator2.__version__ == "3.0.0"


def test_cpp_module_imports():
    """Test that the C++ module can be imported."""
    try:
        from excavator2 import _excavator_core

        assert hasattr(_excavator_core, "__version__")
        assert hasattr(_excavator_core, "hello")
        assert hasattr(_excavator_core, "has_openmp")
    except ImportError as e:
        pytest.skip(f"C++ module not built yet: {e}")


def test_cpp_hello():
    """Test the hello function from C++ module."""
    try:
        from excavator2._excavator_core import hello

        message = hello()
        assert "EXCAVATOR2" in message
        assert "C++" in message
    except ImportError:
        pytest.skip("C++ module not built yet")


def test_cpp_version_match():
    """Test that Python and C++ versions match."""
    try:
        import excavator2
        from excavator2._excavator_core import __version__ as cpp_version

        assert excavator2.__version__ == cpp_version
    except ImportError:
        pytest.skip("C++ module not built yet")
