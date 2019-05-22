''' Import libraries '''
import sys# Standard library
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anm
sys.path.append('C:\\Python-Library')# My library path on my Lenovo-E540
#sys.path.append('C:\\Users\\56402\\OneDrive\\Research\\Source\\Python')# My library path on my NState pro
from PythonLibrary import *

''' Auxiliary routines '''
def Get_input(source):
    with open(source,'r') as f:
        data=f.readlines()
    jobtype=str(data[3].strip())
    NState=int(data[5].strip())
    mass=float(data[7].strip())
    x0=float(data[9].strip())
    TotalTime=float(data[11].strip())
    left=float(data[13].strip())
    right=float(data[15].strip())
    maxdx=float(data[17].strip())
    maxdt=float(data[19].strip())
    p0=float(data[21].strip())
    dt=float(data[23].strip())
    Absorbed=bool(data[25].strip())
    p0left=float(data[27].strip())
    p0right=float(data[29].strip())
    dp0=float(data[31].strip())
    return jobtype,NState,mass,x0,TotalTime,left,right,maxdx,maxdt,p0,dt,Absorbed,p0left,p0right,dp0

def Get_ParametersUsed(source):
    with open(source,'r') as f:
        data=f.readlines()
    NGrid=int(data[0].strip())
    dx=float(data[1].strip())
    actualtime=float(data[2].strip())
    lt=int(data[3].strip())
    lp0=int(data[4].strip())
    NState=int(data[5].strip())
    return NGrid,dx,actualtime,lt,lp0,NState

def Get_Grid(source):
    with open(source,'r') as f:
        data=f.readlines()
    ls=len(data)
    s=numpy.empty(ls)
    for i in range(ls):
        s[i]=float(data[i])
    return s

def Get_psy(source,lx,NState,lt):
    with open(source,'r') as f:
        data=f.readlines()
    psy=numpy.empty((lx,NState,lt),dtype=complex)
    index=0
    for i in range(lt):
        for j in range(NState):
            for k in range(lx):
                psy[k,j,i]=float(data[index].strip())+1j*float(data[index+1].strip())
                index=index+2
    return psy

def Get_TR(source,NState):
    with open(source,'r') as f:
        data=f.readlines()
    ld=len(data)
    lk=int(ld/(1+2*NState))
    k=numpy.empty(lk)
    tran=numpy.empty((lk,NState))
    ref=numpy.empty((lk,NState))
    index=-1
    for i in range(lk):
        index=index+1
        k[i]=float(data[index].strip())
        for j in range(NState):
            index=index+1
            tran[i,j]=float(data[index].strip())
            index=index+1
            ref[i,j]=float(data[index+1].strip())
        print(file=f)
    return k,tran,ref

def Get_Wigner(source,lx,lk,NState,lt):
    with open(source,'r') as f:
        data=f.readlines()
    wigner=numpy.empty((lx,lk,NState,lt))
    index=0
    for i in range(lt):
        for j in range(NState):
            for k in range(lk):
                for l in range(lx):
                    wigner[l,k,j,i]=float(data[index].strip())
                    index=index+1
    return wigner

def Get_SMD(source,NState,lt,nSMD):
    with open(source,'r') as f:
        data=f.readlines()
    SMD=numpy.empty((nSMD,NState,lt))
    index=0
    for ii in range(lt):
        for i in range(NState):
            for j in range(nSMD):
                SMD[j,i,ii]=float(data[index].strip())
                index=index+1
    return SMD

def Animate_Density(left,right,x,t,psy,NState,speed=1.0,title='Density',xlabel='x [a.u.]',ylabel='Density [a.u.]',FileName='Density'):
    den=numpy.multiply(numpy.conj(psy),psy)
    fig,ax=plt.subplots(NState,1,squeeze=False)
    line=[]
    for j in range(NState):
       temp, =ax[j,0].plot(x,den[:,NState-1-j,0])
       ax[j,0].set_title(title)
       ax[j,0].set_xlim(left,right)
       ax[j,0].set_xlabel(xlabel)
       ax[j,0].set_ylim(0,numpy.amax(den))
       ax[j,0].set_ylabel(ylabel)
       line.append(temp)
    def animate(i):
        for j in range(NState):
            line[j].set_ydata(den[:,NState-1-j,i])
        return line
    #animate each a.u. time as 0.01s
    ani=anm.FuncAnimation(fig,animate,t.shape[0],interval=40.0/speed,blit=True)
    ani.save(FileName+'.gif')
    plt.show()

''' Do the job '''
#Read input
jobtype,NState,mass,x0,TotalTime,left,right,maxdx,maxdt,p0,dt,Absorbed,p0left,p0right,dp0=Get_input('OneDimDVR.in')
#Read parameters used in DVR evolution
NGrid,dx,actualtime,lt,lp0,NState=Get_ParametersUsed('ParametersUsed.DVR')

if(jobtype=='NewTrajectory'):
    x=Get_Grid('x.DVR')
    t=Get_Grid('t.DVR')
    psy=Get_psy('Psy.DVR',x.shape[0],NState,lt)
    Animate_Density(left,right,x,t,psy,NState,title='Coordinate space density',FileName='xDensity')
elif(jobtype=='TR-p0'):
    k,tran,ref=Get_TR('TR.DVR',NState)
    with open('TR.txt','w') as f:
        print('k',end='\t',file=f)
        for i in range(NState):
            print('Transmission_'+str(i),end='\t',file=f)
            print('Reflection_'+str(i),end='\t',file=f)
        print(file=f)
        for i in range(lk):
            print(k[i],end='\t',file=f)
            for j in range(NState):
                print(tran[i,j],end='\t',file=f)
                print(ref[i,j],end='\t',file=f)
            print(file=f)
    for i in range(NState):
        for j in range(2):
            plt.subplot(NState,2,2*i+1)
            plt.plot(k,tran[:,i])
            plt.title('Transmission on state'+str(i))
            plt.subplot(NState,2,2*i+2)
            plt.plot(k,ref[:,i])
            plt.title('Reflection on state'+str(i))
    plt.show()
elif(jobtype=='SMD'):
    SMDOrder=2#Currently only support 2
    t=Get_Grid('t.DVR')
    nSMD=int(SMDOrder*(SMDOrder+3)/2)+1
    SMD=Get_SMD('SMD.DVR',NState,lt,nSMD)
    with open('SMD.txt','w') as f:
        print('t/a.u.',end='\t',file=f)
        for i in range(NState):
            print('q'+str(i)+'/a.u.','p'+str(i)+'/a.u.','sigmaq'+str(i)+'/a.u.','rho'+str(i),'sigmap'+str(i)+'/a.u.','pop'+str(i),sep='\t',end='\t',file=f)
        print(file=f)
        for i in range(lt):
            print(t[i],end='\t',file=f)
            for k in range(NState):
                for j in range(nSMD):
                    print(SMD[j,k,i],end='\t',file=f)
            print(file=f)
elif(jobtype=='pRepresentation'):
    k=Get_Grid('k.DVR')
    t=Get_Grid('t.DVR')
    psy=Get_psy('Phi.DVR',k.shape[0],NState,lt)
    Animate_Density(k[0],k[k.shape[0]-1],k,t,psy,NState,title='Momentum space density',xlabel='p [a.u.]',FileName='pDensity')
elif(jobtype=='WignerRepresentation'):
    x=Get_Grid('x.DVR')
    k=Get_Grid('k.DVR')
    t=Get_Grid('t.DVR')
    wigner=Get_Wigner('Wigner.DVR',NGrid,NGrid,NState,lt)
    for i in range(t.shape[0]):
        ax=plt.subplot(111,projection='3d')
        ax.scatter(x[j],y[i],z[j,i])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.savefig('WignerAtTime'+str(i)+'.jpg')
        ax.close()
else:
    print('Unkown job type?')