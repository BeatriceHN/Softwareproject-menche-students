# src/geometry.py
# Minimal geometry generator stubs for organelle placement.
# Returns point clouds for simple parametric organelles.
import numpy as np

def sample_sphere(center, radius, n_points=1000, rng=None):
    rng = np.random.default_rng(rng)
    u = rng.random(n_points)
    v = rng.random(n_points)
    theta = 2 * np.pi * u
    phi = np.arccos(2 * v - 1)
    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]
    return np.vstack((x, y, z)).T

def sample_ellipsoid(center, axes, n_points=1000, rng=None):
    # axes = (a, b, c)
    rng = np.random.default_rng(rng)
    u = rng.random(n_points)
    v = rng.random(n_points)
    theta = 2 * np.pi * u
    phi = np.arccos(2 * v - 1)
    a, b, c = axes
    x = a * np.sin(phi) * np.cos(theta) + center[0]
    y = b * np.sin(phi) * np.sin(theta) + center[1]
    z = c * np.cos(phi) + center[2]
    return np.vstack((x, y, z)).T

def sample_shell_surface(center, axes, thickness, n_points=1000, rng=None):
    # For membrane-proximal placements: sample two concentric ellipsoids and pick shell points.
    outer = sample_ellipsoid(center, axes, n_points, rng=rng)
    inner_axes = (max(axes[0] - thickness, 0.001), max(axes[1] - thickness, 0.001), max(axes[2] - thickness, 0.001))
    inner = sample_ellipsoid(center, inner_axes, n_points, rng=rng)
    return outer
