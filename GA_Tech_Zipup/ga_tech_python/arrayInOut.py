import numpy as np

def writeArray(array, fileName = 'output.txt', long_fmt = False):
    dim = len(array.shape)
    fileHandle = open(fileName, 'w')

    if (long_fmt):
        fmt_str = '%25.14e %+.14ej'
    else:
        fmt_str = '%15.4f %+.4fj'
    if (dim == 1):
        Ncol = 1
    #    np.savetxt(fileHandle, array, fmt_str*Ncol)
        np.savetxt(fileHandle, array)
    elif(dim == 2):
        Ncol = array.shape[1]
        np.savetxt(fileHandle, array)
    #    np.savetxt(fileHandle, array, fmt_str*Ncol)
    else: print("Error, array not accepted")

    fileHandle.close()

def readArray(fileName = 'output.txt'):
    array = np.loadtxt(fileName, dtype = np.complex128)
    chk_cmplx = np.any(np.iscomplex(array))
    if(not chk_cmplx):
        array = np.copy(np.real(array))

    return array

