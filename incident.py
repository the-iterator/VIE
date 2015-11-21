#******************************************************************************#
# incident.py 
#
#    incident field module, used by VIE solver
#    contains functions:
#       PointDipole
#       PlaneWave
#
# Neil Budko (c) 2015
#******************************************************************************#
import scipy as sci

def PointDipole(Freq,EpsB,Cell,NX,NY,NZ,Xs,Ys,Zs,XPol,YPol,ZPol):
    """Returns the incident field Ein over the computational domain,
       produced by a point dipole source polarized as (XPol,YPol,ZPol)
       and located at (Xs,Ys,Zs)
    """
    c0=299792458.0 # speed of light in vacuum
    Mu0=4.0*sci.pi*1.0e-7 # vacuum permeability
    Eps0=1.0/(Mu0*c0*c0) # vacuum permittivity
    Field=sci.zeros((NX,NY,NZ,3),complex)
    Omega=2.0*sci.pi*Freq
    EtaB=-1.0j*Omega*Eps0*EpsB
    ZetaB=-1.0j*Omega*Mu0
    KB=Omega*sci.sqrt(Eps0*EpsB*Mu0)
    # 3D arrays of x,y,z coordinates
    xx,yy,zz=sci.mgrid[0:NX*Cell:Cell,0:NY*Cell:Cell,0:NZ*Cell:Cell]
    # 3D arrays of distances
    dd=sci.sqrt((xx-Xs)**2+(yy-Ys)**2+(zz-Zs)**2)
    dd2=dd*dd
    Xd=xx-Xs
    Yd=yy-Ys
    Zd=zz-Zs
    # 3D arrays of components of the Q-matrix
    Q11=(Xd*Xd)/dd2
    Q12=(Xd*Yd)/dd2
    Q13=(Xd*Zd)/dd2
    Q21=Q12
    Q22=(Yd*Yd)/dd2
    Q23=(Yd*Zd)/dd2
    Q31=Q13
    Q32=Q23
    Q33=(Zd*Zd)/dd2
    QJ1=Q11*XPol+Q12*YPol+Q13*ZPol
    QJ2=Q21*XPol+Q22*YPol+Q23*ZPol
    QJ3=Q31*XPol+Q32*YPol+Q33*ZPol

    dd2=dd*dd
    dd3=dd2*dd

    Field[0:NX,0:NY,0:NZ,0]=sci.exp(1.0j*KB*dd)*\
                (sci.divide((3.0*QJ1-XPol),(4.0*EtaB*sci.pi*dd3))\
                +sci.divide((3.0*QJ1-XPol)*(-1.0j*KB),(4.0*EtaB*sci.pi*dd2))\
                +sci.divide((QJ1-XPol),(4.0*ZetaB*sci.pi*dd)))

    Field[0:NX,0:NY,0:NZ,1]=sci.exp(1.0j*KB*dd)*\
                (sci.divide((3.0*QJ2-YPol),(4.0*EtaB*sci.pi*dd3))\
                +sci.divide((3.0*QJ2-YPol)*(-1.0j*KB),(4.0*EtaB*sci.pi*dd2))\
                +sci.divide((QJ2-YPol),(4.0*ZetaB*sci.pi*dd)))

    Field[0:NX,0:NY,0:NZ,2]=sci.exp(1.0j*KB*dd)*\
                (sci.divide((3.0*QJ3-ZPol),(4.0*EtaB*sci.pi*dd3))\
                +sci.divide((3.0*QJ3-ZPol)*(-1.0j*KB),(4.0*EtaB*sci.pi*dd2))\
                +sci.divide((QJ3-ZPol),(4.0*ZetaB*sci.pi*dd)))

    return Field

def PlaneWave(Freq,EpsB,Cell,NX,NY,NZ,Kx,Ky,Kz,Xs,Ys,Zs,XPol,YPol,ZPol):
    """Returns the incident field Ein over the computational domain,
       produced by a plane wave polarized as (XPol,YPol,Zpol),
       propagating in the direction of the vector (Kx,Ky,Kz),
       and located (strating at) at (Xs,Ys,Zs)
    """
    c0=299792458.0 # speed of light in vacuum
    Mu0=4.0*sci.pi*1.0e-7 # vacuum permeability
    Eps0=1.0/(Mu0*c0*c0) # vacuum permittivity
    Field=sci.zeros((NX,NY,NZ,3),complex)
    Omega=2.0*sci.pi*Freq
    KB=Omega*sci.sqrt(Eps0*EpsB*Mu0)
    # 3D arrays of x,y,z coordinates
    xx,yy,zz=sci.mgrid[0:NX*Cell:Cell,0:NY*Cell:Cell,0:NZ*Cell:Cell]
    # 3D arrays of distances
    KdotX=Kx*(xx-Xs)+Ky*(yy-Ys)+Kz*(zz-Zs)

    Field[0:NX,0:NY,0:NZ,0]=XPol*sci.exp(1.0j*KB*KdotX)

    Field[0:NX,0:NY,0:NZ,1]=YPol*sci.exp(1.0j*KB*KdotX)

    Field[0:NX,0:NY,0:NZ,2]=ZPol*sci.exp(1.0j*KB*KdotX)

    return Field
