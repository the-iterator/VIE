#******************************************************************************#
# green.py 
#
#    Green's tensor module, used by VIE solver
#    contains functions:
#       greentensor
#
# Neil Budko (c) 2015
#******************************************************************************#
import scipy as sci
from numpy import fft

def greentensor(Freq,EpsB,Cell,NX,NY,NZ):
    """Returns the Fourier transform GF of the
       circular extension of the Green's tensor array
    """
    c0=299792458.0 # speed of light in vacuum
    Mu0=4.0*sci.pi*1.0e-7 # vacuum permeability
    Eps0=1.0/(Mu0*c0*c0) # vacuum permittivity
    Omega=2.0*sci.pi*Freq
    EtaB=-1.0j*Omega*Eps0*EpsB
    ZetaB=-1.0j*Omega*Mu0
    KB=Omega*sci.sqrt(Eps0*EpsB*Mu0)
    G=sci.zeros((NX,NY,NZ,3,3),complex)
    GC=sci.zeros((NX*2,NY*2,NZ*2,3,3),complex)
    GF=sci.zeros((NX*2,NY*2,NZ*2,3,3),complex)
    # 3D arrays of x,y,z coordinates
    xx,yy,zz=sci.mgrid[0:NX*Cell:Cell,0:NY*Cell:Cell,0:NZ*Cell:Cell]
    dd=sci.zeros((NX,NY,NZ),complex)
    alpha=sci.zeros((NX,NY,NZ),complex)
    beta=sci.zeros((NX,NY,NZ),complex)
    Q11=sci.zeros((NX,NY,NZ),complex)
    Q12=sci.zeros((NX,NY,NZ),complex)
    Q13=sci.zeros((NX,NY,NZ),complex)
    Q21=sci.zeros((NX,NY,NZ),complex)
    Q22=sci.zeros((NX,NY,NZ),complex)
    Q23=sci.zeros((NX,NY,NZ),complex)
    Q31=sci.zeros((NX,NY,NZ),complex)
    Q32=sci.zeros((NX,NY,NZ),complex)
    Q33=sci.zeros((NX,NY,NZ),complex)
    # 3D arrays of distances
    dd=sci.sqrt((xx)**2+(yy)**2+(zz)**2)
    dd2=dd*dd
    # 3D arrays of components of the Q-matrix
    Q11=sci.divide(xx*xx,dd2)
    Q12=sci.divide(xx*yy,dd2)
    Q13=sci.divide(xx*zz,dd2)
    Q21=Q12
    Q22=sci.divide(yy*yy,dd2)
    Q23=sci.divide(yy*zz,dd2)
    Q31=Q13
    Q32=Q23
    Q33=sci.divide(zz*zz,dd2)
    # alpha and beta scalar multipliers
    alpha=sci.divide(sci.exp(1.0j*KB*dd),4.0*sci.pi*dd)*\
    (-KB**2.0 - sci.divide(1.0j*3.0*KB,dd) + sci.divide(3.0,dd2))
    beta=sci.divide(sci.exp(1.0j*KB*dd),4.0*sci.pi*dd)*\
    (KB**2.0 + sci.divide(1.0j*KB,dd) - sci.divide(1.0,dd2))
    print('All divisions by zero are done')
    # Green's tensor without self-patch
    G[:,:,:,0,0]=Q11*alpha+beta
    G[:,:,:,0,1]=Q12*alpha
    G[:,:,:,0,2]=Q13*alpha
    G[:,:,:,1,0]=Q21*alpha
    G[:,:,:,1,1]=Q22*alpha+beta
    G[:,:,:,1,2]=Q23*alpha
    G[:,:,:,2,0]=Q31*alpha
    G[:,:,:,2,1]=Q32*alpha
    G[:,:,:,2,2]=Q33*alpha+beta
    G=G*(Cell**3) # multiplying by the elementary volume
    # self-patch
    G[0,0,0,0,0]=(2./3.)*(1.-1.j*KB*Cell*((3./(4.*sci.pi))**(1./3.)))*\
    sci.exp(1.j*KB*Cell*((3./(4.*sci.pi))**(1./3.)))-1.0
    G[0,0,0,0,1]=0.
    G[0,0,0,0,2]=0.
    G[0,0,0,1,0]=0.
    G[0,0,0,1,1]=G[0,0,0,0,0]
    G[0,0,0,1,2]=0.
    G[0,0,0,2,0]=0.
    G[0,0,0,2,1]=0.
    G[0,0,0,2,2]=G[0,0,0,0,0]
    #Circular extension of G
    GC[0:NX,0:NY,0:NZ,:,:]=G
    DeltaOp=sci.eye(3,3)
    s=0
    while s<=2:
      ss=0
      while ss<=2:
        GC[NX+1:,0:NY,0:NZ,s,ss]=(1-2*DeltaOp[0,s])*(1-2*DeltaOp[0,ss])*G[:0:-1,:,:,s,ss]
        GC[NX+1:,NY+1:,0:NZ,s,ss]=(1-2*DeltaOp[0,s])*(1-2*DeltaOp[0,ss])*\
                             (1-2*DeltaOp[1,s])*(1-2*DeltaOp[1,ss])*G[:0:-1,:0:-1,:,s,ss]
        GC[NX+1:,NY+1:,NZ+1:,s,ss]=(1-2*DeltaOp[0,s])*(1-2*DeltaOp[0,ss])*\
                             (1-2*DeltaOp[1,s])*(1-2*DeltaOp[1,ss])*\
                             (1-2*DeltaOp[2,s])*(1-2*DeltaOp[2,ss])*G[:0:-1,:0:-1,:0:-1,s,ss]
        GC[0:NX,NY+1:,0:NZ,s,ss]=(1-2*DeltaOp[1,s])*(1-2*DeltaOp[1,ss])*G[:,:0:-1,:,s,ss]
        GC[0:NX,NY+1:,NZ+1:,s,ss]=(1-2*DeltaOp[1,s])*(1-2*DeltaOp[1,ss])*\
                             (1-2*DeltaOp[2,s])*(1-2*DeltaOp[2,ss])*G[:,:0:-1,:0:-1,s,ss]
        GC[0:NX,0:NY,NZ+1:,s,ss]=(1-2*DeltaOp[2,s])*(1-2*DeltaOp[2,ss])*G[:,:,:0:-1,s,ss]
        GC[NX+1:,0:NY,NZ+1:,s,ss]=(1-2*DeltaOp[0,s])*(1-2*DeltaOp[0,ss])*\
                             (1-2*DeltaOp[2,s])*(1-2*DeltaOp[2,ss])*G[:0:-1,:,:0:-1,s,ss]
        ss=ss+1
      s=s+1
    # FFT of the Green's tensor array
    s=0
    while s<=2:
      ss=0
      while ss<=2:
        GF[:,:,:,s,ss]=fft.fftn(sci.squeeze(GC[:,:,:,s,ss]))
        ss=ss+1
      s=s+1

    return GF
