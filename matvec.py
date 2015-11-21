#******************************************************************************#
# matvec.py 
#
#    matrix-vector multiplication module, used by VIE solver
#    contains functions:
#       Au
#
# Neil Budko (c) 2015
#******************************************************************************#
import scipy as sci
from numpy import fft

def Au(U,GF,EpsArr,NX,NY,NZ):
    """Returns the result of matrix-vector multiplication
       by the system matrix A=I-GX
    """
    # reshaping input vector into 4-D array
    Uarr=sci.reshape(U,(NX,NY,NZ,3))
    # extended zero-padded arrays
    Uext=sci.zeros((2*NX,2*NY,2*NZ,3),complex)
    Vext=sci.zeros((2*NX,2*NY,2*NZ,3),complex)
    Jext=sci.zeros((2*NX,2*NY,2*NZ,3),complex)
    JFext=sci.zeros((2*NX,2*NY,2*NZ,3),complex)
    Uext[0:NX,0:NY,0:NZ,:]=Uarr
    # contrast current array
    s=0
    while s<=2:
        Jext[0:NX,0:NY,0:NZ,s]=Uext[0:NX,0:NY,0:NZ,s]*(EpsArr[0:NX,0:NY,0:NZ]-1.0)
        JFext[:,:,:,s]=fft.fftn(sci.squeeze(Jext[:,:,:,s]))
        s=s+1
    Vext[:,:,:,0]=Uext[:,:,:,0]-\
    fft.ifftn(sci.squeeze(sci.multiply(GF[:,:,:,0,0],JFext[:,:,:,0])+\
                          sci.multiply(GF[:,:,:,0,1],JFext[:,:,:,1])+\
                          sci.multiply(GF[:,:,:,0,2],JFext[:,:,:,2])))
    Vext[:,:,:,1]=Uext[:,:,:,1]-\
    fft.ifftn(sci.squeeze(sci.multiply(GF[:,:,:,1,0],JFext[:,:,:,0])+\
                          sci.multiply(GF[:,:,:,1,1],JFext[:,:,:,1])+\
                          sci.multiply(GF[:,:,:,1,2],JFext[:,:,:,2])))
    Vext[:,:,:,2]=Uext[:,:,:,2]-\
    fft.ifftn(sci.squeeze(sci.multiply(GF[:,:,:,2,0],JFext[:,:,:,0])+\
                          sci.multiply(GF[:,:,:,2,1],JFext[:,:,:,1])+\
                          sci.multiply(GF[:,:,:,2,2],JFext[:,:,:,2])))
    # reshaping output into column vector
    V=sci.reshape(Vext[0:NX,0:NY,0:NZ,:],(NX*NY*NZ*3,1))

    return V
