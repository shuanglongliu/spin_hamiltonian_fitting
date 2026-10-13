import numpy as np

### Change pair from [2, 1] to [1, 2] when the latter is convenient.
ostring = """
  - pair: [2, 1]
    conf1:   [{:6.1f}, {:6.1f}, {:6.1f}, {:6.1f}, {:18.9f}]
    conf2:   [{:6.1f}, {:6.1f}, {:6.1f}, {:6.1f}, {:18.9f}]
"""

def do_one_pair(i, j):
    data = np.loadtxt("energies.dat")
    ii = min(i, j)
    jj = max(i, j)
    for i in range(9):
        j = i+9
        print(ostring.format(data[i, 2*(ii-1)], data[i, 2*ii-1], data[i, 2*(jj-1)], data[i, 2*jj-1], data[i, -1], \
                             data[j, 2*(ii-1)], data[j, 2*ii-1], data[j, 2*(jj-1)], data[j, 2*jj-1], data[j, -1] ))

if __name__ == "__main__":

    do_one_pair(2, 4)

