#******************************************************************************#
# addobject.py 
#
#    geometry module, modifies the permittivity array, used by VIE solver
#    contains functions:
#       brick
#       randomfill
#       sphere
#
# Neil Budko (c) 2015
#******************************************************************************#
import scipy as sci

def brick(EpsArr,Cell,CornerX,CornerY,CornerZ,Lx,Ly,Lz,EpsRel):
  """This function adds a homogeneous parallelepiped and
     returns the 3D array of relative permittivity"""
  NxBrickStart=int(round(CornerX/Cell))
  NyBrickStart=int(round(CornerY/Cell))
  NzBrickStart=int(round(CornerZ/Cell))
  NxBrickEnd=int(round((CornerX+Lx)/Cell))
  NyBrickEnd=int(round((CornerY+Ly)/Cell))
  NzBrickEnd=int(round((CornerZ+Lz)/Cell))

  EpsArr[NxBrickStart:NxBrickEnd,NyBrickStart:NyBrickEnd,NzBrickStart:NzBrickEnd]=EpsRel

  return

def randomfill(EpsArr,Cell,CornerX,CornerY,CornerZ,Lx,Ly,Lz,EpsRel,EpsDev):
  """This function adds a random perturbation to the permittivity
     inside a given parallelepiped and returns the modified 3D array
     of relative permittivity"""
  NxBrickStart=int(round(CornerX/Cell))
  NyBrickStart=int(round(CornerY/Cell))
  NzBrickStart=int(round(CornerZ/Cell))
  NxBrickEnd=int(round((CornerX+Lx)/Cell))
  NyBrickEnd=int(round((CornerY+Ly)/Cell))
  NzBrickEnd=int(round((CornerZ+Lz)/Cell))

  EpsArr[NxBrickStart:NxBrickEnd,NyBrickStart:NyBrickEnd,NzBrickStart:NzBrickEnd]=EpsRel-\
    EpsDev+2.*EpsDev*sci.random.rand(NxBrickEnd-NxBrickStart,
    NyBrickEnd-NyBrickStart,NzBrickEnd-NzBrickStart)

  return

def sphere(EpsArr,Cell,Xc,Yc,Zc,Radius,EpsRel,NX,NY,NZ):
  """This function adds a homogeneous sphere and
     returns the 3D array of relative permittivity"""
  kk=0
  while kk<NZ:
    jj=0
    Zp=kk*Cell
    while jj<NY:
      Yp=jj*Cell
      ii=0
      while ii<NX:
        Xp=ii*Cell
        DistC=sci.sqrt((Xp-Xc)**2+(Yp-Yc)**2+(Zp-Zc)**2)
        if DistC<=Radius:
          EpsArr[ii,jj,kk]=EpsRel
        ii=ii+1
      jj=jj+1
    kk=kk+1

  return
