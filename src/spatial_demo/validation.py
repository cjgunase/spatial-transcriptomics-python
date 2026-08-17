"""Validation helpers kept outside the notebook for reuse."""

from __future__ import annotations

from typing import Any


def describe_spatial_adata(adata: Any) -> dict[str, object]:
    """Return a compact description and fail early if coordinates are missing.

    Parameters
    ----------
    adata
        An AnnData-like object with observations, variables, and ``obsm``.
    """
    if "spatial" not in adata.obsm:
        raise ValueError("Expected spatial coordinates in adata.obsm['spatial']")

    coordinates = adata.obsm["spatial"]
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("Spatial coordinates must have shape (n_observations, 2)")

    return {
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "coordinate_shape": tuple(coordinates.shape),
        "observation_columns": list(adata.obs.columns),
    }
