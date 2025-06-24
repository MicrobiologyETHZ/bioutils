import numpy as np
import pytest
# Replace 'your_module' with the actual module containing the rarefy function
from bioinf_utils.rnaseq_qc import rarefy


def test_basic_rarefaction():
    x = np.array([10, 20, 30, 40, 50])
    rarefied = rarefy(x, depth=100, iterations=1, seed=42)
    assert len(rarefied) == len(x)
    assert np.all(np.isnan(rarefied) == False)


def test_rarefaction_depth_greater_than_total_counts():
    x = np.array([10, 20, 30, 40, 50])
    rarefied = rarefy(x, depth=300, iterations=1, seed=42)
    assert len(rarefied) == len(x)
    assert np.all(np.isnan(rarefied))


def test_multiple_iterations():
    x = np.array([10, 20, 30, 40, 50])
    rarefied = rarefy(x, depth=50, iterations=5, seed=42)
    assert len(rarefied) == len(x)
    assert np.all(np.isnan(rarefied) == False)


def test_zero_counts():
    x = np.array([0, 0, 0, 0, 0])
    rarefied = rarefy(x, depth=10, iterations=1, seed=42)
    assert len(rarefied) == len(x)
    assert np.all(np.isnan(rarefied))


def test_seed_reproducibility():
    x = np.array([10, 20, 30, 40, 50])
    rarefied_1 = rarefy(x, depth=50, iterations=1, seed=42)
    rarefied_2 = rarefy(x, depth=50, iterations=1, seed=42)
    np.testing.assert_allclose(rarefied_1, rarefied_2, rtol=1e-6)
