import matplotlib.pyplot as plt
import matplotlib.animation as anm
import numpy

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
    QHDOrder=int(data[23].strip())
    p0left=float(data[25].strip())
    p0right=float(data[27].strip())
    dp0=float(data[29].strip())
    return jobtype,steptype,surface,mass,x0,totaltime,left,right,maxdx,maxdt,p0,QHDOrder,p0left,p0right,dp0

def Get_ParametersUsed(source):
    with open(source,'r') as f:
        data=f.readlines()
    NGrid=int(data[0].strip())
    NAbsorbGrid=int(data[1].strip())
    lx=int(data[2].strip())
    dx=float(data[3].strip())
    actualtime=float(data[4].strip())
    lt=int(data[5].strip())
    dt=float(data[6].strip())
    lp0=int(data[7].strip())
    surface=int(data[8].strip())
    return NGrid,NAbsorbGrid,lx,dx,actualtime,lt,dt,lp0,surface

def Get_Grid(source):
    with open(source,'r') as f:
        data=f.readlines()
    ls=len(data)
    s=numpy.empty(ls)
    for i in range(ls):
        s[i]=float(data[i])
    return s

def Get_TR(source,surface):
    with open(source,'r') as f:
        data=f.readlines()
    ld=len(data)
    lk=int(ld/(1+4*surface))
    k=numpy.empty(lk)
    tran=numpy.empty((lk,surface))
    ref=numpy.empty((lk,surface))
    atran=numpy.empty((lk,surface))
    aref=numpy.empty((lk,surface))
    index=-1
    for i in range(lk):
        index=index+1
        k[i]=float(data[index].strip())
        for j in range(surface):
            index=index+1
            tran[i,j]=float(data[index].strip())
            ref[i,j]=float(data[index+1].strip())
            atran[i,j]=float(data[index+2].strip())
            aref[i,j]=float(data[index+3].strip())
            index=index+3
        print(file=f)
    return k,tran,ref

def Get_psy(source,lx,lt,surface):
    with open(source,'r') as f:
        data=f.readlines()
    psy=numpy.empty((lx,lt,surface),dtype=complex)
    index=0
    for ii in range(surface):
        for i in range(lt):
            for j in range(lx):
                psy[j,i,ii]=float(data[index].strip())+1j*float(data[index+1].strip())
                index=index+2
    return psy

def Animate_Density(left,right,x,t,psy,surface):
    den=numpy.multiply(numpy.conj(psy),psy)
    fig,ax=plt.subplots(surface,1,squeeze=False)
    line=[]
    for j in range(surface):
       temp, =ax[j,0].plot(x,den[:,0,surface-1-j])
       ax[j,0].set_xlim(left,right)
       ax[j,0].set_title(str(0))
       ax[j,0].set_ylim(0,numpy.amax(den))
       line.append(temp)
    def animate(i):
        for j in range(surface):
            print(i)
            ax[j,0].set_title(str(i))
            line[j].set_ydata(den[:,i,surface-1-j])
        return line
    #animate each a.u. time as 0.01s
    ani=anm.FuncAnimation(fig,animate,t.shape[0],interval=10*(t[1]-t[0]),blit=True)
    plt.show()

def Get_QHD(source,lt,surface,nqhd):
    with open(source,'r') as f:
        data=f.readlines()
    qhd=numpy.empty((nqhd,lt,surface))
    index=0
    for ii in range(surface):
        for i in range(lt):
            for j in range(nqhd):
                qhd[j,i,ii]=float(data[index].strip())
                index=index+1
    return qhd