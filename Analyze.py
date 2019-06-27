''' Import library '''
import sys# Standard library
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anm
#sys.path.append('/home-4/yshen57@jhu.edu/Library/Python-Library')# My library path on MARCC
sys.path.append('C:\\Python-Library')# My library path on my Lenovo-E540
#sys.path.append('C:\\Users\\56402\\OneDrive\\Research\\Source\\Python')# My library path on my surface pro
from PythonLibrary import *

''' Auxiliary routine '''
def Get_input(source):
    with open(source,'r') as f: data=f.readlines()
    jobtype=str(data[3].strip()); NState=int(data[5].strip()); mass=float(data[7].strip())
    x0=float(data[9].strip()); left=float(data[11].strip()); right=float(data[13].strip())
    maxdx=float(data[15].strip()); maxdt=float(data[17].strip())
    p0=float(data[19].strip()); TotalTime=float(data[21].strip()); dt=float(data[23].strip())
    Scatter=bool(data[25].strip())
    pleft=float(data[27].strip()); pright=float(data[29].strip()); dp=float(data[31].strip())
    skipx=int(data[33].strip())
    return jobtype,NState,mass,x0,left,right,maxdx,maxdt,p0,TotalTime,dt,Scatter,pleft,pright,dp,skipx

def Get_wavefunction(source,lx,NState,lt):
    with open(source,'r') as f: data=f.readlines()
    psy=numpy.empty((lx,NState,lt),dtype=complex)
    index=0
    for i in range(lt):
        for j in range(NState):
            for k in range(lx):
                psy[k,j,i]=float(data[index].strip())+1j*float(data[index+1].strip())
                index=index+2
    return psy

def Get_TR(source,NState):
    with open(source,'r') as f: data=f.readlines()
    lp=int(len(data)/(1+2*NState)); p=numpy.empty(lp)
    tran=numpy.empty((lp,NState)); ref=numpy.empty((lp,NState))
    index=0
    for i in range(lp):
        p[i]=float(data[index].strip())
        index=index+1
        for j in range(NState):
            tran[i,j]=float(data[index].strip()); ref[i,j]=float(data[index+1].strip())
            index=index+2
    return p,tran,ref

def Get_Wigner(source,lwx,lp,NState,lt):
    with open(source,'r') as f: data=f.readlines()
    wigner=numpy.empty((lwx,lp,NState,lt))
    index=0
    for i in range(lt):
        for j in range(NState):
            for k in range(lp):
                for l in range(lwx):
                    wigner[l,k,j,i]=float(data[index].strip())
                    index=index+1
    return wigner

def Get_SMD(source,NState,lt,nSMD):
    with open(source,'r') as f: data=f.readlines()
    SMD=numpy.empty((nSMD,NState,lt))
    index=0
    for ii in range(lt):
        for i in range(NState):
            for j in range(nSMD):
                SMD[j,i,ii]=float(data[index].strip())
                index=index+1
    return SMD

''' Do the job '''
#Read input
jobtype,NState,mass,x0,left,right,maxdx,maxdt,p0,TotalTime,dt,Scatter,pleft,pright,dp,skipx=Get_input('OneDimDVR.in')
if(jobtype=='NewTrajectory'):
    x=General.Get_Array('x.out'); t=General.Get_Array('t.out')
    psy=Get_wavefunction('Psy.out',x.shape[0],NState,t.shape[0]); den=numpy.multiply(numpy.conj(psy),psy)
    if(NState>1): Visualization.Animate2D_MultiLevel(t,x,den,\
        title='Coordinate space density',xlabel='x [a.u.]',ylabel='density',save=True,FileName='xDensity',show=False)
    else: Visualization.Animate2D(t,x,den[:,0,:],\
        title='Coordinate space density',xlabel='x [a.u.]',ylabel='density',save=True,FileName='xDensity',show=False)
elif(jobtype=='TR-p0'):
    p,tran,ref=Get_TR('TR.out',NState)
    with open('TR.txt','w') as f:
        print('p',end='\t',file=f)
        for i in range(NState):
            print('Transmission_'+str(i),end='\t',file=f)
            print('Reflection_'+str(i),end='\t',file=f)
        print(file=f)
        for i in range(p.shape[0]):
            print(p[i],end='\t',file=f)
            for j in range(NState):
                print(tran[i,j],end='\t',file=f)
                print(ref[i,j],end='\t',file=f)
            print(file=f)
    for i in range(NState):
        for j in range(2):
            plt.subplot(NState,2,2*i+1)
            plt.plot(p,tran[:,i])
            plt.title('Transmission on state'+str(i))
            plt.subplot(NState,2,2*i+2)
            plt.plot(p,ref[:,i])
            plt.title('Reflection on state'+str(i))
    plt.show()
elif(jobtype=='SMD'):
    t=General.Get_Array('t.out')
    SMDOrder=2#Currently only support 2
    nSMD=int(SMDOrder*(SMDOrder+3)/2)+1
    SMD=Get_SMD('SMD.out',NState,t.shape[0],nSMD)
    with open('SMD.txt','w') as f:
        print('t/a.u.',end='\t',file=f)
        for i in range(NState):
            print('q'+str(i)+'/a.u.','p'+str(i)+'/a.u.','sigmaq'+str(i)+'/a.u.','covariance'+str(i),'sigmap'+str(i)+'/a.u.','population'+str(i),sep='\t',end='\t',file=f)
        print(file=f)
        for i in range(t.shape[0]):
            print(t[i],end='\t',file=f)
            for k in range(NState):
                for j in range(nSMD):
                    print(SMD[j,k,i],end='\t',file=f)
            print(file=f)
elif(jobtype=='pRepresentation'):
    p=General.Get_Array('p.out'); t=General.Get_Array('t.out')
    phi=Get_wavefunction('Phi.out',p.shape[0],NState,t.shape[0]); den=numpy.multiply(numpy.conj(phi),phi)
    if(NState>1): Visualization.Animate2D_MultiLevel(t,p,den,\
        title='Momentum space density',xlabel='p [a.u.]',ylabel='density',save=True,FileName='pDensity',show=False)
    else: Visualization.Animate2D(t,p,den[:,0,:],\
        title='Momentum space density',xlabel='p [a.u.]',ylabel='density',save=True,FileName='pDensity',show=False)
elif(jobtype=='WignerDistribution'):
    x=General.Get_Array('Wigner_x.out'); p=General.Get_Array('Wigner_p.out'); t=General.Get_Array('t.out')
    wigner=Get_Wigner('Wigner.out',x.shape[0],p.shape[0],NState,t.shape[0])
    [xtemp,ptemp]=numpy.meshgrid(x,p); Xtemp=xtemp.ravel(); Ptemp=ptemp.ravel()
    X=numpy.empty((Xtemp.shape[0],t.shape[0])); P=numpy.empty((Ptemp.shape[0],t.shape[0]))
    for i in range(t.shape[0]):
        X[:,i]=Xtemp; P[:,i]=Ptemp
    WIGNER=numpy.empty((int(x.shape[0]*p.shape[0]),NState,t.shape[0]))
    for k in range(t.shape[0]):
        for j in range(NState):
            index=0
            for i in range(p.shape[0]):
                WIGNER[index:index+x.shape[0],j,k]=wigner[:,i,j,k]
                index=index+x.shape[0]
    for j in range(NState): Visualization.Animate3DSurface(t,X,P,WIGNER[:,j,:],\
        title='Wigner distribution on state '+str(j),xlabel='x [a.u.]',ylabel='p [a.u.]',zlabel='density',\
        colormap='seismic',save=True,FileName='Wigner'+str(j),show=False)
else:
    print('Unkown job type?')