#******************************************************************************#
# matvecpar.py
#
#    matrix-vector multiplication module, used by VIE solver
#
#    multi-core version
#
#    contains functions:
#       Au
#
# Neil Budko (c) 2015
#******************************************************************************#
import scipy as sci
from numpy import fft
from multiprocessing import Pool

def FFT(ArrIn):
    ArrOut = fft.fftn(ArrIn)
    return ArrOut

def IFFT(ArrIn):
    ArrOut = fft.ifftn(ArrIn)
    return ArrOut

def Au(U,GF,EpsArr,NX,NY,NZ,p):
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
    for s in range(3):
        Jext[0:NX,0:NY,0:NZ,s]=Uext[0:NX,0:NY,0:NZ,s]*(EpsArr[0:NX,0:NY,0:NZ]-1.0)

    JFext[:,:,:,0],JFext[:,:,:,1],JFext[:,:,:,2]=p.map(FFT,\
    [sci.squeeze(Jext[:,:,:,0]),sci.squeeze(Jext[:,:,:,1]),\
    sci.squeeze(Jext[:,:,:,2])])

    V1,V2,V3 = p.map(IFFT,\
    [sci.squeeze(sci.multiply(GF[:,:,:,0,0],JFext[:,:,:,0])+\
                          sci.multiply(GF[:,:,:,0,1],JFext[:,:,:,1])+\
                          sci.multiply(GF[:,:,:,0,2],JFext[:,:,:,2])),\
     sci.squeeze(sci.multiply(GF[:,:,:,1,0],JFext[:,:,:,0])+\
                          sci.multiply(GF[:,:,:,1,1],JFext[:,:,:,1])+\
                          sci.multiply(GF[:,:,:,1,2],JFext[:,:,:,2])),\
     sci.squeeze(sci.multiply(GF[:,:,:,2,0],JFext[:,:,:,0])+\
                          sci.multiply(GF[:,:,:,2,1],JFext[:,:,:,1])+\
                          sci.multiply(GF[:,:,:,2,2],JFext[:,:,:,2]))])

    Vext[:,:,:,0]=Uext[:,:,:,0]-V1
    Vext[:,:,:,1]=Uext[:,:,:,1]-V2
    Vext[:,:,:,2]=Uext[:,:,:,2]-V3
    # reshaping output into column vector
    V=sci.reshape(Vext[0:NX,0:NY,0:NZ,:],(NX*NY*NZ*3,1))

    return V
