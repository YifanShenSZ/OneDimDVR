''' Import library '''
import sys; import numpy
import matplotlib.pyplot as plt
#sys.path.append('/home-4/yshen57@jhu.edu/Library/Python-Library')# My library path on MARCC
#sys.path.append('C:\\Python-Library')# My library path on my Lenovo-E540
sys.path.append('C:\\Users\\56402\\OneDrive\\Research\\Source\\Python')# My library path on my surface pro
from PythonLibrary import *

''' Auxiliary routine '''
def Get_input(source):
    with open(source,'r') as f: data=f.readlines()
    jobtype=str(data[3].strip()); NState=int(data[5].strip())
    return jobtype,NState

def Get_wavefunction(source,lq,NState,lt):
    with open(source,'r') as f: data=f.readlines()
    psy=numpy.empty((lq,NState,lt),dtype=complex)
    index=0
    for i in range(lt):
        for j in range(NState):
            for k in range(lq):
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

def Get_Wigner(source,lwq,lp,NState,lt):
    with open(source,'r') as f: data=f.readlines()
    wigner=numpy.empty((lwq,lp,NState,lt))
    index=0
    for i in range(lt):
        for j in range(NState):
            for k in range(lp):
                for l in range(lwq):
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
jobtype,NState=Get_input('OneDimDVR.in')
if(jobtype=='NewTrajectory'):# Animate coordinate space density
    q=General.Get_Array('q.out'); t=General.Get_Array('t.out')
    psy=Get_wavefunction('Psy.out',q.shape[0],NState,t.shape[0]); den=numpy.multiply(numpy.conj(psy),psy)
    if(NState>1): Visualization.Animate2D_MultiLevel(t,q,den,\
        title='Coordinate space density',xlabel='q [a.u.]',ylabel='density',save=True,FileName='qDensity',show=False)
    else: Visualization.Animate2D(t,q,den[:,0,:],\
        title='Coordinate space density',xlabel='q [a.u.]',ylabel='density',save=True,FileName='qDensity',show=False)
elif(jobtype=='TR-p0'):# Reformat output and plot
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
elif(jobtype=='pRepresentation'):# Animate momentum space density
    p=General.Get_Array('p.out'); t=General.Get_Array('t.out')
    phi=Get_wavefunction('Phi.out',p.shape[0],NState,t.shape[0]); den=numpy.multiply(numpy.conj(phi),phi)
    if(NState>1): Visualization.Animate2D_MultiLevel(t,p,den,\
        title='Momentum space density',xlabel='p [a.u.]',ylabel='density',save=True,FileName='pDensity',show=False)
    else: Visualization.Animate2D(t,p,den[:,0,:],\
        title='Momentum space density',xlabel='p [a.u.]',ylabel='density',save=True,FileName='pDensity',show=False)
elif(jobtype=='WignerDistribution'):# Animate Wigner distribution
    q=General.Get_Array('Wigner_q.out'); p=General.Get_Array('Wigner_p.out'); t=General.Get_Array('t.out')
    wigner=Get_Wigner('Wigner.out',q.shape[0],p.shape[0],NState,t.shape[0])
    [qtemp,ptemp]=numpy.meshgrid(q,p); Qtemp=qtemp.ravel(); Ptemp=ptemp.ravel()
    Q=numpy.empty((Qtemp.shape[0],t.shape[0])); P=numpy.empty((Ptemp.shape[0],t.shape[0]))
    for i in range(t.shape[0]):
        Q[:,i]=Qtemp; P[:,i]=Ptemp
    WIGNER=numpy.empty((int(q.shape[0]*p.shape[0]),NState,t.shape[0]))
    for k in range(t.shape[0]):
        for j in range(NState):
            index=0
            for i in range(p.shape[0]):
                WIGNER[index:index+q.shape[0],j,k]=wigner[:,i,j,k]
                index=index+q.shape[0]
    for j in range(NState): Visualization.Animate3DSurface(t,Q,P,WIGNER[:,j,:],\
        title='Wigner distribution on state '+str(j),xlabel='q [a.u.]',ylabel='p [a.u.]',zlabel='density',\
        colormap='seismic',save=True,FileName='Wigner'+str(j),show=False)
elif(jobtype=='SMD'):# Reformat output and plot
    t=General.Get_Array('t.out')
    SMDOrder=2#Currently only support 2
    nSMD=int(SMDOrder*(SMDOrder+3)/2)+1
    SMD=Get_SMD('SMD.out',NState,t.shape[0],nSMD)
    if(NState>1):
        with open('SMD.txt','w') as f:
            print('t/a.u.',end='\t',file=f)
            for i in range(NState):
                print('<q>'+str(i)+'/a.u.','<p>'+str(i)+'/a.u.','sigmaq'+str(i)+'/a.u.','correlation'+str(i),'sigmap'+str(i)+'/a.u.','population'+str(i),sep='\t',end='\t',file=f)
            print(file=f)
            for i in range(t.shape[0]):
                print(t[i],end='\t',file=f)
                for k in range(NState):
                    for j in range(nSMD):
                        print(SMD[j,k,i],end='\t',file=f)
                print(file=f)
    else:
        with open('SMD.txt','w') as f:
            print('t/a.u.','<q>/a.u.','<p>/a.u.','sigmaq/a.u.','correlation','sigmap/a.u.',sep='\t',file=f)
            for i in range(t.shape[0]):
                print(t[i],end='\t',file=f)
                for j in range(nSMD-1):
                    print(SMD[j,0,i],end='\t',file=f)
                print(file=f)
        Visualization.Plot2D(t,SMD[0,0,:],title='<q>-t')
else: print('Unkown job type?')