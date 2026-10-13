"""Geometry helpers shared by the magnetic anisotropy job directories and data_processing/.

Spherical/cartesian conversion, reference-frame changes, a discretized sphere, and the
per-site magnetic moments used to write MAGMOM.

Angle conventions (when ``rad=False``, degrees):
    theta : polar angle measured from the +z axis, in [0, 180]
    phi   : azimuthal angle in the xy-plane from +x toward +y, in [0, 360)
"""

import math
import numpy as np


def appendSpherical_np(xyz):
    """Append [r, theta, phi] columns (radians) to an (N, 3) cartesian array.

    theta is the polar angle from the +z axis.
    """
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:, 0]**2 + xyz[:, 1]**2
    ptsnew[:, 3] = np.sqrt(xy + xyz[:, 2]**2)
    ptsnew[:, 4] = np.arctan2(np.sqrt(xy), xyz[:, 2])
    ptsnew[:, 5] = np.arctan2(xyz[:, 1], xyz[:, 0])
    return ptsnew


def sph2cart(sph, rad=False):
    """Convert spherical [r, theta, phi] to cartesian [x, y, z].

    Note: ``sph`` is mutated in place when ``rad=False`` (angles converted
    to radians).
    """
    if not rad:
        sph[1:3] = np.deg2rad(sph[1:3])
    z = sph[0] * math.cos(sph[1])
    x = sph[0] * math.sin(sph[1]) * math.cos(sph[2])
    y = sph[0] * math.sin(sph[1]) * math.sin(sph[2])
    return [x, y, z]


def cart2sph(xyz, rad=False):
    """Convert cartesian [x, y, z] to spherical [r, theta, phi]."""
    ptsnew = appendSpherical_np(np.array([xyz]))
    sph = [ptsnew[0][3]]
    if rad:
        sph.append(ptsnew[0][4])
        sph.append(ptsnew[0][5])
    else:
        sph.append(ptsnew[0][4] / math.pi * 180)
        sph.append(ptsnew[0][5] / math.pi * 180)
    return sph


def get_opposite_direction_sph(direction):
    theta, phi = direction
    theta = 180 - theta
    phi = phi + 180
    phi = phi % 360
    if theta == 0 or theta == 180:
        phi = 0
    return (theta, phi)

def get_angle(a, b, rad=False):
    cosab = np.dot(a, b)/np.linalg.norm(a)/np.linalg.norm(b)
    if cosab > 1:
        cosab = 1
    elif cosab < -1:
        cosab = -1
    if rad:
        theta = math.acos(cosab)
    else:
        theta = 180*math.acos(cosab)/math.pi
    return theta

def get_angle_sph(a, b, rad=False):
    a = sph2cart(a, rad)
    b = sph2cart(b, rad)
    theta = get_angle(a, b, rad)
    return theta

def change_frame_sph(vi_sph, emati, ematf, rad=False):
    """Re-express a spherical vector from one orthonormal frame in another.

    Parameters
    ----------
    vi_sph : [r, theta, phi] in the initial frame.
    emati  : 3x3 array; rows are the initial basis vectors in global coords.
    ematf  : 3x3 array; rows are the final basis vectors in global coords.

    Returns
    -------
    vf_sph : the same vector expressed as [r, theta, phi] in the final frame.
    """
    emati_t = np.transpose(emati)
    ematf_t = np.transpose(ematf)
    ematf_t_inv = np.linalg.inv(ematf_t)

    vi_cart = sph2cart(vi_sph, rad=rad)
    vf_cart = np.matmul(emati_t, vi_cart)
    vf_cart = np.matmul(ematf_t_inv, vf_cart)

    vf_sph = cart2sph(vf_cart, rad=rad)

    if rad:
        pi = np.pi
        twopi = 2 * np.pi
    else:
        pi = 180
        twopi = 360

    # On the poles, phi is undefined; pick a canonical value.
    if abs(vf_sph[1]) < 1.e-9 or abs(vf_sph[1] - pi) < 1.e-9:
        vf_sph[2] = 0
    # Wrap phi into [0, 2*pi).
    if vf_sph[2] < 0:
        vf_sph[2] = vf_sph[2] + twopi
    elif vf_sph[2] > twopi:
        vf_sph[2] = vf_sph[2] - twopi

    return vf_sph


def change_frame_angles(angles, emati, ematf):
    """Re-express arrays of directions (thetas, phis), in radians, from frame emati in frame ematf.

    Returns
    -------
    (thetas_prime, phis_prime) : arrays of the same directions in the final frame.
    """
    thetas, phis = angles
    thetas_prime = []
    phis_prime = []
    for theta, phi in zip(thetas, phis):
        dumb, theta_prime, phi_prime = change_frame_sph([1, theta, phi], emati, ematf, rad=True)
        thetas_prime.append(theta_prime)
        phis_prime.append(phi_prime)
    return (np.array(thetas_prime), np.array(phis_prime))


def get_emat_local(z_direction=(0, 0)):
    """Build a local right-handed orthonormal frame whose +z axis points
    along ``z_direction = (theta, phi)`` (degrees) in the global frame.

    The local +x axis is chosen so that the global +x axis lies in the
    local xz-plane (with the special case that aligns with the global x
    axis handled explicitly to avoid a degenerate cross product).
    """
    if z_direction == (90, 0):
        # Local z axis along global x: swap axes (x -> z, y -> x, z -> y).
        emat_local = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    else:
        ez = sph2cart([1, z_direction[0], z_direction[1]])
        ey = np.cross(ez, [1, 0, 0])
        ey = ey / np.linalg.norm(ey)
        ex = np.cross(ey, ez)
        emat_local = [ex, ey, ez]
    return np.array(emat_local)


class magmoms:
    def __init__(self, n_site=1, magmoms=[1], directions=[[0,0]]):
        self.n_site = n_site
        if n_site != len(magmoms):
            magmoms = [1 for i in range(n_site)]
        self.magmoms = magmoms
        if n_site != len(directions):
            directions = [[0, 0] for i in range(n_site)]
        self.directions = directions
        self.get_directions_magmom_to_spin()
        self.get_m_mat_car()
        self.get_mstring()
    
    def get_directions_magmom_to_spin(self, reverse=False):
        if reverse:
            self.directions = []
            for i in range(self.n_site):
                self.directions.append( get_opposite_direction_sph(self.directions_of_spin[i]) )
        else:
            self.directions_of_spin = []
            for i in range(self.n_site):
                self.directions_of_spin.append( get_opposite_direction_sph(self.directions[i]) )

    def get_m_mat_car(self):
        m_mat_car = []
        for i in range(self.n_site):
            m_vec_car = sph2cart([self.magmoms[i], self.directions[i][0], self.directions[i][1]])
            m_mat_car.append(m_vec_car)
        self.m_mat_car = np.array(m_mat_car)

    def get_mstring(self):
        mstring = ""
        for i in range(self.n_site):
            m_vec_car = self.m_mat_car[i]
            m_vec_car = (3*"{:15.10f}").format(*m_vec_car)
            mstring = mstring + m_vec_car + "   \\\n"
        self.mstring = mstring

    def set_directions(self, directions):
        self.directions = directions
        self.get_directions_magmom_to_spin()
        self.get_m_mat_car()
        self.get_mstring()
        #print(self.directions)

    def set_directions_of_spin(self, directions_of_spin):
        self.directions_of_spin = directions_of_spin
        self.get_directions_magmom_to_spin(reverse=True)
        self.get_m_mat_car()
        self.get_mstring()
        #print(self.directions)

    def set_magmoms(self, magmoms):
        self.magmoms = magmoms
        self.get_m_mat_car()
        self.get_mstring()


class sphere:
    """A discretized unit sphere indexed by (theta, phi) bins.

    The two poles each map to a single direction (i_direction = 0 and
    n_direction - 1), so the total number of distinct directions is
    ``n_phi * (n_theta - 2) + 2``.
    """

    def __init__(self, n_theta=1801, n_phi=3600):
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.dtheta = 180 / (self.n_theta - 1)
        self.dphi = 360 / self.n_phi
        self.n_direction = self.n_phi * (self.n_theta - 2) + 2
        self.i_direction = 0
        self.theta = 0.0
        self.phi = 0.0

    def get_i_direction(self):
        """Set ``self.i_direction`` from the current ``self.theta``/``self.phi``."""
        i_theta = int(round(self.theta / self.dtheta))
        i_phi = int(round(self.phi / self.dphi))
        if i_phi == self.n_phi:
            i_phi = 0
        if i_theta == 0:
            i_direction = 0
        elif i_theta == self.n_theta - 1:
            i_direction = self.n_direction - 1
        else:
            i_direction = self.n_phi * (i_theta - 1) + i_phi + 1
        self.i_direction = i_direction

    def get_theta(self):
        """Set ``self.theta`` from the current ``self.i_direction``."""
        if self.i_direction == 0:
            self.theta = 0.0
        elif self.i_direction == self.n_direction - 1:
            self.theta = 180.0
        else:
            i_theta = (self.i_direction - 1) // self.n_phi + 1
            self.theta = i_theta * self.dtheta

    def get_phi(self):
        """Set ``self.phi`` from the current ``self.i_direction``."""
        if self.i_direction == 0:
            self.phi = 0.0
        elif self.i_direction == self.n_direction - 1:
            self.phi = 0.0
        else:
            i_phi = (self.i_direction - 1) % self.n_phi
            self.phi = i_phi * self.dphi

    def set_i_direction(self, i_direction=0):
        """Clamp and set the linear index, then refresh theta and phi."""
        if i_direction < 0:
            i_direction = 0
        elif i_direction >= self.n_direction:
            i_direction = self.n_direction - 1
        self.i_direction = i_direction
        self.get_theta()
        self.get_phi()

    def set_angle(self, direction):
        """Clamp and set (theta, phi), then refresh the linear index and snap
        theta/phi to the nearest grid point."""
        theta, phi = direction
        if theta < 0:
            theta = 0.0
        elif theta > 180:
            theta = 180.0
        if phi < 0:
            phi = 0.0
        elif phi > 360:
            phi = 360.0
        self.theta = theta
        self.phi = phi
        self.get_i_direction()
        self.get_theta()
        self.get_phi()

    def set_n_theta_n_phi(self, n_theta=181, n_phi=360):
        """Re-grid the sphere; preserves ``self.theta``/``self.phi``."""
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.dtheta = 180 / (self.n_theta - 1)
        self.dphi = 360 / self.n_phi
        self.n_direction = self.n_phi * (self.n_theta - 2) + 2
        self.get_i_direction()
