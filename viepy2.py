#******************************************************************************#
#
#  viepy2.py
#
#  VIE solver
#
#  Neil Budko (c) 2012-2015
#  n.v.budko@gmail.com
#
#******************************************************************************#
#
# This file contains the main program
#
# Standard modules: scipy, matplotlib, time, sys, os, inspect
# Extra modules: green, incident, scattered, addobject, visualization, matvec
#
# Input: see basic_input.py (inside ../Input folder)
# Output: matlab files (see ../Output folder)
# Visualization: 2D slices, GMRES convergence plot
#
# To run in background on Linux:
# at now (press 'Enter')
#   python viepy2.py input_file (press 'Enter')
#   (press 'Ctrl+D')
#
# To monitor the output:
# tail -f stdout.log
#
# To see if the job is running, in terminal
# atq
#
#******************************************************************************#
import scipy as sci
import scipy.io as scio
from scipy import sparse
from scipy.sparse import linalg as ssla
from scipy import linalg as scila
import matplotlib.pyplot as plt
import time
import sys
import os
import inspect
import incident # module for computing icident field
import scattered # module for computing scattered field
import green # module for computing Green's tensor
import addobject # module for adding objects to the domain
import visualization as vis # visualization module
import matvec # matrix-vector multiplication module

# redirecting the stdout and stderr to files without buffering
#sys.stdout = open('stdout.log', 'w', 0)
#sys.stderr = open('stderr.log', 'w', 0)

print '\n*************'
print 'PROGRAM START'
print 'VIE solver \nNeil Budko (c) 2012-2015\n'

#inputfile = sys.argv[1]
inputfile = 'basic_input.inp'

c0 = 299792458.0 # speed of light in vacuum
Mu0 = 4.0*sci.pi*1e-7 # vacuum permeability
Eps0 = 1.0/(Mu0*c0*c0) # vacuum permittivity

#******************************************************************************#
# Loading input parameters
#******************************************************************************#
# determining the current script directory

myDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
fileName = myDir+'/Input/'+inputfile # name of the file with basic input
f = open(fileName,'r') # open data file for reading
print 'data file:', fileName,'\n'
print 'Problem parameters:\n'
for line in f.readlines(): # reading and executing lines from input file
    exec(line)
f.close()

print 'NX =',NX,'NY =',NY,'NZ =',NZ  #
print 'Freq =',Freq, '(frequency, [Hz])'  #
print 'Lambda0 =',Lambda0, '(vacuum wavelength, [m])'  #
print 'EpsRelMax =',EpsRelMax, 'maximum relative permittivity'  #
print 'LambdaMed =',LambdaMed, '(minimum wavelength, [m])'  #
print 'Kappa =',Kappa, '(points per vacuum wavelength)'  #
print 'Cell =',Cell, '(discretization cell, [m])'  #
print 'KappaMed =',LambdaMed/Cell, '(points per medium wavelength)'  #
print 'Source: ',source 
print '(PolX,PolY,PolZ)=(',PolX,PolY,PolZ,') polarization of source'  #
print '(Xs,Ys,Zs)=(',Xs,Ys,Zs,') coordinates of source, [m]'  #

#******************************************************************************#
# Computing Ein
#******************************************************************************#
t0 = time.clock()
Ein = sci.zeros([NX,NY,NZ,3],complex)
if source == 'point dipole':
    Ein = incident.PointDipole(Freq,EpsRelB,Cell,NX,NY,NZ,Xs,Ys,Zs,PolX,PolY,PolZ)
elif source == 'plane wave':
    Ein = incident.PlaneWave(Freq,EpsRelB,Cell,NX,NY,NZ,Kx,Ky,Kz,Xs,Ys,Zs,PolX,PolY,PolZ)
else:
    print('Error: unknown source type')
t1 = time.clock()
# saving Ein in a *.mat file
outFile = myDir+'/Output/IncidentField' # file name
scio.savemat(outFile, {'Ein':Ein})
Bvec = sci.reshape(Ein,(NX*NY*NZ*3,1))
deltaT = t1-t0
sys.stdout.write('Incident field computed in '+str(deltaT)+' seconds\n')

# Showing Ein
figNum = 1
fieldName = 'E^{in}'
SliceX = round(NX/2.)
SliceY = round(NY/2.)
SliceZ = round(NZ/2.)
vis.sliceFieldArr(sci.real(Ein),NX,NY,NZ,SliceX,SliceY,SliceZ,figNum,fieldName)

#******************************************************************************#
# Computing G tensor
#******************************************************************************#
t0 = time.clock()
GF = green.greentensor(Freq,EpsRelB,Cell,NX,NY,NZ)
t1 = time.clock()
# saving GF in a *.mat file
#outFile = myDir+'/Output/GreenTensor' # file name
#scio.savemat(outFile, {'GF':GF})
deltaT = t1-t0
print 'GF tensor computed in',deltaT, 'seconds'

#******************************************************************************#
# Constructing an object
#******************************************************************************#
# setting all permittivity to EpsRelB
EpsArr = sci.ones((NX,NY,NZ),complex) # array with all elements 1+0j
EpsArr = EpsArr*EpsRelB

# adding a homogeneous brick to the domain (repeat if necessary)
CornerX,CornerY,CornerZ = 0.1*Xdim, 0.1*Ydim, 0.1*Zdim # corner coordinates of the brick, in [m] (use Lambda0)
Lx,Ly,Lz = 0.8*Xdim, 0.8*Ydim, 0.8*Zdim # dimensions of the brick, in [m] (use Lambda0)
EpsRel = 2.0 # permittivity of the brick, can be complex
addobject.brick(EpsArr,Cell,CornerX,CornerY,CornerZ,Lx,Ly,Lz,EpsRel)

# adding a homogeneous ball to the domain (repeat if necessary)
Xc,Yc,Zc = 0.5*Xdim, 0.5*Ydim, 0.5*Zdim
Radius = 0.3*Zdim
EpsRel = 6.0
addobject.sphere(EpsArr,Cell,Xc,Yc,Zc,Radius,EpsRel,NX,NY,NZ)

# adding a homogeneous ball to the domain (repeat if necessary)
Xc,Yc,Zc = 0.5*Xdim, 0.5*Ydim, 0.5*Zdim
Radius = 0.2*Zdim
EpsRel = 4.0
addobject.sphere(EpsArr,Cell,Xc,Yc,Zc,Radius,EpsRel,NX,NY,NZ)

# Plotting EpsArrClean
figNum = 2
vis.sliceEpsArr(EpsArr,NX,NY,NZ,figNum)

# saving EpsArrClean in a *.mat file
outFile = myDir+'/Output/EpsRelArr' # file name
scio.savemat(outFile, {'EpsArr':EpsArr})

#******************************************************************************#
# Computing the total field
#******************************************************************************#
# parameters for the GMRES solver
tol = 1e-6 # relative residual tolerance
restart = 300 # inner iterations
maxiter = 300 # max total iterations

# declearing the operator
def AuClean(U):
    V = matvec.Au(U,GF,EpsArr,NX,NY,NZ)
    return V
Aop = ssla.LinearOperator((NX*NY*NZ*3,NX*NY*NZ*3),matvec = AuClean,dtype = complex)

# preparing to read the residual vector
relres = []
nb = scila.norm(Bvec)
def resnorm(residual):
    relres.append(residual)
    global it
    print 'it =', it, ' rel. res. =', residual
    it = it+1
    return

print 'Computing total field...'
print 'Starting GMRES with: tol =', tol,' restart =',restart, ' maxiter =', maxiter
t0 = time.clock()
it = 1
(Vsol,info) = ssla.gmres(Aop,Bvec, tol=tol, restart=restart, maxiter=maxiter, callback=resnorm)
t1 = time.clock()
V0 = Vsol # for use as an initial guess
Etot = sci.reshape(Vsol,(NX,NY,NZ,3))
# saving Etot in a *.mat file
outFile = myDir+'/Output/TotalField' # file name
scio.savemat(outFile, {'Etot':Etot})

deltaT = t1-t0
print 'GMRES finished in', deltaT, 'seconds'
print 'Info:', info

# Plotting residual
plt.ion()
fig = plt.figure(3)
plt.semilogy(relres)
plt.title('GMRES residual')
fig.canvas.draw()
plt.ioff()

# Slicing E
figNum = 4
fieldName = 'E'
vis.sliceFieldArr(sci.absolute(Etot),NX,NY,NZ,SliceX,SliceY,SliceZ,figNum,fieldName)


#******************************************************************************#
# Computing force density
#******************************************************************************#
print 'Computing force density...'
t0 = time.clock()


t1 = time.clock()
deltaT = t1-t0
print 'Force density computed in', deltaT, 'seconds'

#plt.show() # show and keep all the figures untill closed
print '\nPROGRAM END :-)'
print '***************\n'

#******************************************************************************#
