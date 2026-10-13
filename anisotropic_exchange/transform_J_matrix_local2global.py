import numpy as np
from scipy.spatial.transform import Rotation as R

def do(emat, J, rotmat):

    """
    The local reference frame is (rotmat*emat^T)^T.
    """

    emat_ = np.transpose(np.matmul(rotmat, np.transpose(emat)))

    J_ = np.matmul(np.transpose(emat_), np.matmul(J, emat_))

    print("J matrix in the global reference frame is:")
    print("{:12.6f} {:12.6f} {:12.6f}".format(*J_[0]))
    print("{:12.6f} {:12.6f} {:12.6f}".format(*J_[1]))
    print("{:12.6f} {:12.6f} {:12.6f}".format(*J_[2]))

if __name__ == "__main__":

    emat = np.array( [ \
        [ 0.173909,    0.984761,   -0.001465], \
        [-0.984726,    0.173915,    0.008295], \
        [ 0.008423,    0.000000,    0.999965] ] )

    J = np.array([ \
      [ 1.674408,  0.003371,  0.040520], \
      [-0.003049,  1.671754,  0.032918], \
      [ 0.001758, -0.001655,  1.651589] ])

    r = R.from_euler('z',   0, degrees=True)
    rotmat = r.as_matrix()
    do(emat, J, rotmat)

    r = R.from_euler('z', 120, degrees=True)
    rotmat = r.as_matrix()
    do(emat, J, rotmat)

    r = R.from_euler('z', 240, degrees=True)
    rotmat = r.as_matrix()
    do(emat, J, rotmat)
