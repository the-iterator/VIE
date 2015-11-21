#******************************************************************************#
# visualization.py 
#
#    module for plots and images, used by VIE solver
#    contains functions:
#       plotEpsArr
#       sliceEpsArr
#       sliceFieldArr
#
# Neil Budko (c) 2015
#******************************************************************************#


import scipy as sci
import numpy as np
import matplotlib.pyplot as plt

def plotEpsArr(EpsArr,NX,NY,NZ,figNum,hold,label,clear):
    """ Visualizing object with cuts 
    """
    cutxR=sci.real(EpsArr[:,round(NY/2.0),round(NZ/2.0)])
    cutxI=sci.imag(EpsArr[:,round(NY/2.0),round(NZ/2.0)])
    cutyR=sci.real(EpsArr[round(NX/2.0),:,round(NZ/2.0)])
    cutyI=sci.imag(EpsArr[round(NX/2.0),:,round(NZ/2.0)])
    cutzR=sci.real(EpsArr[round(NX/2.0),round(NY/2.0),:])
    cutzI=sci.imag(EpsArr[round(NX/2.0),round(NY/2.0),:])
    
    plt.ion()
    fig=plt.figure(figNum)
    fig.clear()
    if clear=='yes':
        fig.clear()

    plt.hold(hold)
    
    plt.subplot(231)
    plt.plot(cutxR,label)
    plt.title('real($\\varepsilon$), x-axis cut')
    plt.subplot(232)
    plt.plot(cutyR,label)
    plt.title('real($\\varepsilon$), y-axis cut')
    plt.subplot(233)
    plt.plot(cutzR,label)
    plt.title('real($\\varepsilon$), z-axis cut')
    plt.subplot(234)
    plt.plot(cutxI,label)
    plt.title('imag($\\varepsilon$), x-axis cut')
    plt.subplot(235)
    plt.plot(cutyI,label)
    plt.title('imag($\\varepsilon$), y-axis cut')
    plt.subplot(236)
    plt.plot(cutzI,label)
    plt.title('imag($\\varepsilon$), z-axis cut')
    
    plt.hold('False')
    fig.canvas.draw()
    plt.ioff()
    
    return

def sliceEpsArr(EpsArr,NX,NY,NZ,figNum):
    """ Visualizing object with slices 
    """
    vxyR=sci.real(EpsArr[:,:,round(NZ/2.0)])
    vxzR=sci.real(EpsArr[:,round(NY/2.0),:])
    vyzR=sci.real(EpsArr[round(NX/2.0),:,:])
    vxyI=sci.imag(EpsArr[:,:,round(NZ/2.0)])
    vxzI=sci.imag(EpsArr[:,round(NY/2.0),:])
    vyzI=sci.imag(EpsArr[round(NX/2.0),:,:])
    
    plt.ion()
    fig=plt.figure(figNum)
    fig.clear()
    
    plt.subplot(231)
    plt.imshow(vxyR)
    plt.title('real($\\varepsilon$), xy-plane')
    plt.colorbar()
    plt.subplot(232)
    plt.imshow(sci.swapaxes(vxzR,0,1))
    plt.title('real($\\varepsilon$), xz-plane')
    plt.colorbar()
    plt.subplot(233)
    plt.imshow(sci.swapaxes(vyzR,0,1))
    plt.title('real($\\varepsilon$), yz-plane')
    plt.colorbar()
    plt.subplot(234)
    plt.imshow(vxyI)
    plt.title('imag($\\varepsilon$), xy-plane')
    plt.colorbar()
    plt.subplot(235)
    plt.imshow(sci.swapaxes(vxzI,0,1))
    plt.title('imag($\\varepsilon$), xz-plane')
    plt.colorbar()
    plt.subplot(236)
    plt.imshow(sci.swapaxes(vyzI,0,1))
    plt.title('imag($\\varepsilon$), yz-plane')
    plt.colorbar()
    
    fig.canvas.draw()
    plt.ioff()
    
    return
    
def sliceFieldArr(Ein,NX,NY,NZ,SliceX,SliceY,SliceZ,figNum,fieldName):
    """ Visualizing field with slices 
    """
    vxy1=Ein[:,:,SliceZ,0]
    vxy2=Ein[:,:,SliceZ,1]
    vxy3=Ein[:,:,SliceZ,2]
    vyz1=Ein[SliceX,:,:,0]
    vyz2=Ein[SliceX,:,:,1]
    vyz3=Ein[SliceX,:,:,2]
    vxz1=Ein[:,SliceY,:,0]
    vxz2=Ein[:,SliceY,:,1]
    vxz3=Ein[:,SliceY,:,2]
    plt.ion()
    fig=plt.figure(figNum)
    fig.clear()
    plt.subplot(331)
    plt.imshow(sci.real(vxy1))
    plt.title('real($'+fieldName+'_{x}$), xy-plane')
    plt.colorbar()
    plt.subplot(332)
    plt.imshow(sci.real(sci.swapaxes(vxz1,0,1)))
    plt.title('real($'+fieldName+'_{x}$), xz-plane')
    plt.colorbar()
    plt.subplot(333)
    plt.imshow(sci.real(sci.swapaxes(vyz1,0,1)))
    plt.title('real($'+fieldName+'_{x}$), yz-plane')
    plt.colorbar()
    plt.subplot(334)
    plt.imshow(sci.real(vxy2))
    plt.title('real($'+fieldName+'_{y}$), xy-plane')
    plt.colorbar()
    plt.subplot(335)
    plt.imshow(sci.real(sci.swapaxes(vxz2,0,1)))
    plt.title('real($'+fieldName+'_{y}$), xz-plane')
    plt.colorbar()
    plt.subplot(336)
    plt.imshow(sci.real(sci.swapaxes(vyz2,0,1)))
    plt.title('real($'+fieldName+'_{y}$), yz-plane')
    plt.colorbar()
    plt.subplot(337)
    plt.imshow(sci.real(vxy3))
    plt.title('real($'+fieldName+'_{z}$), xy-plane')
    plt.colorbar()
    plt.subplot(338)
    plt.imshow(sci.real(sci.swapaxes(vxz3,0,1)))
    plt.title('real($'+fieldName+'_{z}$), xz-plane')
    plt.colorbar()
    plt.subplot(339)
    plt.imshow(sci.real(sci.swapaxes(vyz3,0,1)))
    plt.title('real($'+fieldName+'_{z}$), yz-plane')
    plt.colorbar()
    
    fig.canvas.draw()
    plt.ioff()
    
    return
