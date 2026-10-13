"""
Spherical/Cartesian conversions and reference-frame changes used by the fitting scripts.

Taken from pylib (shuanglongliu/programs, python/pylib.py) so that this directory has no
external dependency beyond numpy.

Convention: theta is the polar angle (0-180 deg, from +z) and phi the azimuthal angle
(0-360 deg). Angles are in degrees unless rad=True. An emat is a (3, 3) array whose
rows are the basis vectors of a frame in global coordinates.
"""

import math

import numpy as np


def appendSpherical_np(xyz: np.ndarray) -> np.ndarray:
    """Append spherical coordinates (r, theta, phi) to an (N, 3) Cartesian array.

    Returns an (N, 6) array; columns 3-5 hold r, theta (rad), phi (rad).
    phi is normalised to [0, 2 pi); at the poles phi is set to 0.
    """
    xyz = np.asarray(xyz, dtype=float)
    result = np.hstack((xyz, np.zeros((xyz.shape[0], 3))))
    xy2 = xyz[:, 0] ** 2 + xyz[:, 1] ** 2

    result[:, 3] = np.sqrt(xy2 + xyz[:, 2] ** 2)
    result[:, 4] = np.arctan2(np.sqrt(xy2), xyz[:, 2])
    result[:, 5] = np.arctan2(xyz[:, 1], xyz[:, 0])

    # phi is degenerate at the poles
    at_pole = (np.abs(result[:, 4]) < 1e-6) | (np.abs(result[:, 4] - np.pi) < 1e-6)
    result[at_pole, 5] = 0.0
    result[result[:, 5] < 0.0, 5] += 2 * np.pi

    return result


def cart2sph(xyz, rad: bool = False) -> list:
    """Cartesian vector to [r, theta, phi]."""
    row = appendSpherical_np(np.array([xyz]))[0]
    r, theta, phi = row[3], row[4], row[5]
    if not rad:
        theta = math.degrees(theta)
        phi = math.degrees(phi)
    return [r, theta, phi]


def cart2sph_deg(xyz) -> list:
    return cart2sph(xyz, rad=False)


def cart2sph_rad(xyz) -> list:
    return cart2sph(xyz, rad=True)


def sph2cart(sph, rad: bool = False) -> list:
    """[r, theta, phi] to a Cartesian vector. The input is not modified."""
    r, theta, phi = float(sph[0]), float(sph[1]), float(sph[2])
    if not rad:
        theta = math.radians(theta)
        phi = math.radians(phi)
    return [r * math.sin(theta) * math.cos(phi), r * math.sin(theta) * math.sin(phi), r * math.cos(theta)]


def sph2cart_deg(sph) -> list:
    return sph2cart(sph, rad=False)


def sph2cart_rad(sph) -> list:
    return sph2cart(sph, rad=True)


def get_opposite_direction_sph(direction) -> tuple:
    """Antipodal direction of (theta, phi) in degrees; phi is set to 0 at the poles."""
    theta, phi = direction
    theta = 180 - theta
    phi = (phi + 180) % 360
    if theta in (0, 180):
        phi = 0
    return (theta, phi)


def get_emat_local(z_direction=(0, 0)) -> np.ndarray:
    """Local frame whose z axis points along z_direction = (theta, phi) in degrees.

    The local x axis is chosen so that the global x axis lies in the local xz plane.
    When z_direction is along the global x axis, the axes are permuted (x -> z, y -> x, z -> y).
    """
    theta, phi = z_direction

    if abs(theta - 90) < 1e-6 and abs(phi) < 1e-6:
        return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)

    ez = np.array(sph2cart_deg([1, theta, phi]))
    ey = np.cross(ez, [1, 0, 0])
    ey = ey / np.linalg.norm(ey)
    ex = np.cross(ey, ez)
    return np.array([ex, ey, ez])


def change_frame_sph(vi_sph, emati, ematf, rad: bool = False) -> list:
    """Express a vector [r, theta, phi], given in the frame emati, in the frame ematf."""
    vi_cart = sph2cart(vi_sph, rad=rad)
    vf_cart = np.linalg.inv(np.transpose(ematf)) @ (np.transpose(emati) @ vi_cart)
    return cart2sph(vf_cart, rad=rad)


def change_frame_sph_rad(vi_sph, emati, ematf) -> list:
    return change_frame_sph(vi_sph, emati, ematf, rad=True)


def change_frame_sph_deg(vi_sph, emati, ematf) -> list:
    return change_frame_sph(vi_sph, emati, ematf, rad=False)
