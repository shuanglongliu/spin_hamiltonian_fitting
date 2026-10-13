import numpy as np

def read_data():
    data = np.loadtxt("energies.dat")
    return data

def get_As(data):
    m = data.shape[0]
    n = data.shape[1]
    As = []
    for i in range(m):
        A = 0
        for j in range(1, 7):
            A = A + data[i,0] * data[i,j]
        As.append(A)
    return As

def get_sumSS(data, pairs):
    npair = len(pairs)
    As = []
    for i in range(data.shape[0]):
        A = 0
        for ipair in range(npair):
            A = A + data[i, pairs[ipair][0]] * data[i, pairs[ipair][1]]
        As.append(A)
    return As

def get_Js(data, sumSSs, factor):
    # Assumption: n Js and n equations after subtracting E0.
    # sumSSs = [sumSS1, sumSS2]
    # Set up and solve the equation factor Cs Js = dEs
    n = data.shape[0] - 1 # number of Js
    Cs = np.zeros((n, n))
    dEs = []
    for i in range(n):
        for j in range(n):
            Cs[i, j] = sumSSs[j][n] - sumSSs[j][i]
        dEs.append(data[n, -1] - data[i, -1])
    Js = 1000*np.matmul(np.linalg.inv(Cs), dEs)/factor
    Js_wavenumber = 8.06554393738 * Js
    print(("Js = [" + n*" {:6.3f} " + "] meV\n   = [" + n*" {:6.2f} " + "] cm^-1").format(*Js, *Js_wavenumber))
    return (Js, Js_wavenumber)

if __name__ == "__main__":

    data = read_data()
    sumSS1 = get_sumSS(data, [[0,1], [0,2], [0,3], [0,4], [0,5], [0,6]])
    sumSS2 = get_sumSS(data, [[1,2], [1,3], [2,3], [4,5], [4,6], [5,6]])
    Js = get_Js(data, [sumSS1, sumSS2], -2)

