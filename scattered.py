#******************************************************************************#
# scattered.py 
#
#    scattered field module, used by VIE solver
#    contains functions:
#       SingleReceiver
#
# Neil Budko (c) 2015
#******************************************************************************#
import scipy as sci

def SingleReceiver(E,EpsArr,Freq,EpsB,Cell,NX,NY,NZ,Xr,Yr,Zr):
    """Returns the scattered field 3-vector Esc
       at the location of the receiver (Xr,Yr,Zr)
    """
    c0 = 299792458.0 # speed of light in vacuum
    Mu0 = 4.0*sci.pi*1.0e-7 # vacuum permeability
    Eps0 = 1.0/(Mu0*c0*c0) # vacuum permittivity
    Omega = 2.0*sci.pi*Freq
    EtaB = -1.0j*Omega*Eps0*EpsB
    ZetaB = -1.0j*Omega*Mu0
    KB = Omega*sci.sqrt(Eps0*EpsB*Mu0)
    Esc = sci.zeros(3,complex)

    J = sci.zeros((NX,NY,NZ,3),complex)
    for s in range(3):
        J[0:NX,0:NY,0:NZ,s]=E[0:NX,0:NY,0:NZ,s]*(EpsArr[0:NX,0:NY,0:NZ]-1.0)

    G = sci.zeros((NX,NY,NZ,3,3),complex)
    # 3D arrays of x,y,z coordinates
    xx,yy,zz = sci.mgrid[0:NX*Cell:Cell,0:NY*Cell:Cell,0:NZ*Cell:Cell]
    dd = sci.zeros((NX,NY,NZ),complex)
    alpha = sci.zeros((NX,NY,NZ),complex)
    beta = sci.zeros((NX,NY,NZ),complex)
    Q11 = sci.zeros((NX,NY,NZ),complex)
    Q12 = sci.zeros((NX,NY,NZ),complex)
    Q13 = sci.zeros((NX,NY,NZ),complex)
    Q21 = sci.zeros((NX,NY,NZ),complex)
    Q22 = sci.zeros((NX,NY,NZ),complex)
    Q23 = sci.zeros((NX,NY,NZ),complex)
    Q31 = sci.zeros((NX,NY,NZ),complex)
    Q32 = sci.zeros((NX,NY,NZ),complex)
    Q33 = sci.zeros((NX,NY,NZ),complex)
    # 3D arrays of distances
    dd = sci.sqrt((Xr-xx)**2+(Yr-yy)**2+(Zr-zz)**2)
    dd2 = dd*dd
    # 3D arrays of components of the Q-matrix
    Q11 = sci.divide((Xr-xx)*(Xr-xx),dd2)
    Q12 = sci.divide((Xr-xx)*(Yr-yy),dd2)
    Q13 = sci.divide((Xr-xx)*(Zr-zz),dd2)
    Q21 = Q12
    Q22 = sci.divide((Yr-yy)*(Yr-yy),dd2)
    Q23 = sci.divide((Yr-yy)*(Zr-zz),dd2)
    Q31 = Q13
    Q32 = Q23
    Q33 = sci.divide((Zr-zz)*(Zr-zz),dd2)
    # alpha and beta scalar multipliers
    alpha = sci.divide(sci.exp(1.0j*KB*dd),4.0*sci.pi*dd)*\
    (-KB**2.0 - sci.divide(1.0j*3.0*KB,dd) + sci.divide(3.0,dd2))
    beta = sci.divide(sci.exp(1.0j*KB*dd),4.0*sci.pi*dd)*\
    (KB**2.0 + sci.divide(1.0j*KB,dd) - sci.divide(1.0,dd2))
    # Green's tensor
    G[:,:,:,0,0] = Q11*alpha+beta
    G[:,:,:,0,1] = Q12*alpha
    G[:,:,:,0,2] = Q13*alpha
    G[:,:,:,1,0] = Q21*alpha
    G[:,:,:,1,1] = Q22*alpha+beta
    G[:,:,:,1,2] = Q23*alpha
    G[:,:,:,2,0] = Q31*alpha
    G[:,:,:,2,1] = Q32*alpha
    G[:,:,:,2,2] = Q33*alpha+beta
    G = G*(Cell**3) # multiplying by the elementary volume

    Esc[0] = sci.sum(sci.squeeze(sci.multiply(G[:,:,:,0,0],J[:,:,:,0])+\
                                 sci.multiply(G[:,:,:,0,1],J[:,:,:,1])+\
                                 sci.multiply(G[:,:,:,0,2],J[:,:,:,2])))
    Esc[1] = sci.sum(sci.squeeze(sci.multiply(G[:,:,:,1,0],J[:,:,:,0])+\
                                 sci.multiply(G[:,:,:,1,1],J[:,:,:,1])+\
                                 sci.multiply(G[:,:,:,1,2],J[:,:,:,2])))
    Esc[2] = sci.sum(sci.squeeze(sci.multiply(G[:,:,:,2,0],J[:,:,:,0])+\
                                 sci.multiply(G[:,:,:,2,1],J[:,:,:,1])+\
                                 sci.multiply(G[:,:,:,2,2],J[:,:,:,2])))

    return Esc
