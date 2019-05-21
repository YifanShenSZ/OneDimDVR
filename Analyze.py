''' Import libraries '''
import sys# Standard library
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anm
sys.path.append('C:\\Python-Library')# My library path on my Lenovo-E540
#sys.path.append('C:\\Users\\56402\\OneDrive\\Research\\Source\\Python')# My library path on my surface pro
from PythonLibrary import *

''' Auxiliary routines '''
def Get_input(source):
    with open(source,'r') as f:
        data=f.readlines()
    jobtype=str(data[1].strip())
    steptype=str(data[3].strip())
    surface=int(data[5].strip())
    mass=float(data[7].strip())
    x0=float(data[9].strip())
    totaltime=float(data[11].strip())
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
    QHDOrder=int(data[33].strip())
    return jobtype,steptype,surface,mass,x0,totaltime,left,right,maxdx,maxdt,p0,dt,Absorbed,p0left,p0right,dp0,QHDOrder

def Get_ParametersUsed(source):
    with open(source,'r') as f:
        data=f.readlines()
    NGrid=int(data[0].strip())
    dx=float(data[1].strip())
    actualtime=float(data[2].strip())
    lt=int(data[3].strip())
    lp0=int(data[4].strip())
    surface=int(data[5].strip())
    return NGrid,dx,actualtime,lt,lp0,surface

def Get_Grid(source):
    with open(source,'r') as f:
        data=f.readlines()
    ls=len(data)
    s=numpy.empty(ls)
    for i in range(ls):
        s[i]=float(data[i])
    return s

def Get_psy(source,lx,surface,lt):
    with open(source,'r') as f:
        data=f.readlines()
    psy=numpy.empty((lx,surface,lt),dtype=complex)
    index=0
    for i in range(lt):
        for j in range(surface):
            for k in range(lx):
                psy[k,j,i]=float(data[index].strip())+1j*float(data[index+1].strip())
                index=index+2
    return psy

def Get_TR(source,surface):
    with open(source,'r') as f:
        data=f.readlines()
    ld=len(data)
    lk=int(ld/(1+2*surface))
    k=numpy.empty(lk)
    tran=numpy.empty((lk,surface))
    ref=numpy.empty((lk,surface))
    index=-1
    for i in range(lk):
        index=index+1
        k[i]=float(data[index].strip())
        for j in range(surface):
            index=index+1
            tran[i,j]=float(data[index].strip())
            index=index+1
            ref[i,j]=float(data[index+1].strip())
        print(file=f)
    return k,tran,ref

def Get_Wigner(source,lx,lk,surface,lt):
    with open(source,'r') as f:
        data=f.readlines()
    wigner=numpy.empty((lx,lk,surface,lt))
    index=0
    for i in range(lt):
        for j in range(surface):
            for k in range(lk):
                for l in range(lx):
                    wigner[l,k,j,i]=float(data[index].strip())
                    index=index+1
    return wigner

def Get_QHD(source,surface,lt,nqhd):
    with open(source,'r') as f:
        data=f.readlines()
    qhd=numpy.empty((nqhd,surface,lt))
    index=0
    for ii in range(lt):
        for i in range(surface):
            for j in range(nqhd):
                qhd[j,i,ii]=float(data[index].strip())
                index=index+1
    return qhd

def Animate_Density(left,right,x,t,psy,surface):
    den=numpy.multiply(numpy.conj(psy),psy)
    fig,ax=plt.subplots(surface,1,squeeze=False)
    line=[]
    for j in range(surface):
       temp, =ax[j,0].plot(x,den[:,surface-1-j,0])
       ax[j,0].set_xlim(left,right)
       ax[j,0].set_title(str(0))
       ax[j,0].set_ylim(0,numpy.amax(den))
       line.append(temp)
    def animate(i):
        for j in range(surface):
            print(i)
            ax[j,0].set_title(str(i))
            line[j].set_ydata(den[:,surface-1-j,i])
        return line
    #animate each a.u. time as 0.01s
    ani=anm.FuncAnimation(fig,animate,t.shape[0],interval=10*(t[1]-t[0]),blit=True)
    ani.save('Density.gif')
    plt.show()

''' Do the job '''
#Read input
jobtype,steptype,surface,mass,x0,totaltime,left,right,maxdx,maxdt,p0,dt,Absorbed,p0left,p0right,dp0,SMDOrder=Get_input('input')
#Read parameters used in DVR evolution
NGrid,dx,actualtime,lt,lp0,surface=Get_ParametersUsed('ParametersUsed.DVR')

if(jobtype=='NewTrajectory'):
    x=Get_Grid('x.DVR')
    t=Get_Grid('t.DVR')
    psy=Get_psy('Psy.DVR',NGrid,surface,lt)
    Animate_Density(left,right,x,t,psy,surface)
elif(jobtype=='TR-p0'):
    k,tran,ref=Get_TR('TR.DVR',surface)
    with open('TR.txt','w') as f:
        print('k',end='\t',file=f)
        for i in range(surface):
            print('Transmission_'+str(i),end='\t',file=f)
            print('Reflection_'+str(i),end='\t',file=f)
        print(file=f)
        for i in range(lk):
            print(k[i],end='\t',file=f)
            for j in range(surface):
                print(tran[i,j],end='\t',file=f)
                print(ref[i,j],end='\t',file=f)
            print(file=f)
    for i in range(surface):
        for j in range(2):
            plt.subplot(surface,2,2*i+1)
            plt.plot(k,tran[:,i])
            plt.title('Transmission on state'+str(i))
            plt.subplot(surface,2,2*i+2)
            plt.plot(k,ref[:,i])
            plt.title('Reflection on state'+str(i))
    plt.show()
elif(jobtype=='SMD'):
    t=Get_Grid('t.DVR')
    nSMD=int(SMDOrder*(SMDOrder+3)/2)+1
    SMD=Get_SMD('SMD.DVR',surface,lt,nSMD)
    with open('SMD.txt','w') as f:
        print('t',end='\t',file=f)
        for i in range(surface):
            print('x'+str(i),'p'+str(i),'sigmax'+str(i),'rho'+str(i),'sigmap'+str(i),'pop'+str(i),sep='\t',end='\t',file=f)
        print(file=f)
        for i in range(lt):
            print(t[i],end='\t',file=f)
            for k in range(surface):
                for j in range(nSMD):
                    print(SMD[j,k,i],end='\t',file=f)
            print(file=f)
elif(jobtype=='pRepresentation'):
    k=Get_Grid('k.DVR')
    t=Get_Grid('t.DVR')
    psy=Get_psy('Phi.DVR',NGrid,surface,lt)
    Animate_Density(k[0],k[NGrid-1],k,t,psy,surface)
elif(jobtype=='WignerRepresentation'):
    x=Get_Grid('x.DVR')
    k=Get_Grid('k.DVR')
    t=Get_Grid('t.DVR')
    wigner=Get_Wigner('Wigner.DVR',NGrid,NGrid,surface,lt)
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