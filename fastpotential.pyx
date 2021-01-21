import numpy as np
cimport numpy as np
from libc.math cimport abs,pow

cpdef double get_norm(np.ndarray[complex, ndim=2] psi,int N, double dx):
    cdef int i,j
    cdef double norm
    norm = 0.0
    for i in range(N):
        for j in range(N):
            norm = norm + pow(abs(psi[i,j]),2)*dx*dx
    return norm



cpdef np.ndarray[complex,ndim=2] laplace(np.ndarray[complex,ndim=2] psi,int N, double dx):
    cdef np.ndarray[complex,ndim=2] laplace
    cdef int i,j
    laplace = np.zeros((N,N), dtype=complex)
    for i in range(1, N - 1):
        for j in range(1,N-1):
            laplace[i,j] = (psi[i+1,j] + psi[i-1,j] + psi[i,j-1] +\
                                psi[i,j+1] - 4. * psi[i,j]) / (dx * dx)
    laplace[0,:]         = 0
    laplace[N-1,:]       = 0  
    laplace[:,0]         = 0
    laplace[:,N-1]       = 0
    return laplace

cpdef np.ndarray[complex,ndim=2] create_gaussian(int N, np.ndarray x, np.ndarray y, double centerx, double centery, double factor,double sigma):
    cdef int i,j 
    cdef np.ndarray[complex,ndim=2] gauss
    gauss = np.zeros((N,N),dtype=complex)
    for i in range(N):
        for j in range(N):
            gauss[i,j] = factor*np.exp(-(pow((x[i]-centerx),2)+pow((y[j]-centery),2)) / 2 / sigma / sigma)

    return gauss