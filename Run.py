from Analyze import *

#Read input
jobtype,steptype,surface,mass,x0,totaltime,left,right,maxdx,maxdt,p0,QHDOrder,p0left,p0right,dp0=Get_input('input')
#Read parameters used in DVR evolution
NGrid,NAbsorbGrid,lx,dx,actualtime,lt,dt,lp0,surface=Get_ParametersUsed('ParametersUsed.DVR')

if(jobtype=='TR-p0'):
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
elif(jobtype=='NewTrajectory'):
    x=Get_Grid('x.DVR')
    t=Get_Grid('t.DVR')
    psy=Get_psy('Psy.DVR',lx,lt,surface)
    Animate_Density(left,right,x[NAbsorbGrid:NAbsorbGrid+NGrid],t,psy[NAbsorbGrid:NAbsorbGrid+NGrid,:,:],surface)
elif(jobtype=='QHD'):
    step=int(1/dt)
    t=Get_Grid('t.DVR')
    nqhd=int(QHDOrder*(QHDOrder+3)/2+1)
    qhdformer=Get_QHD('QHDformer.DVR',lt,surface,nqhd)
    qhdlatter=Get_QHD('QHDlatter.DVR',lt,surface,nqhd)
    qhdtotal=Get_QHD('QHD.DVR',lt,surface,nqhd)
    for j in range(surface):
        for i in range(2):
            plt.plot(t,qhdtotal[i,:,j])
            plt.title(str(i)+' on surface '+str(j))
            plt.show()
    with open('QHDformer.txt','w') as f:
        print('t',end='\t',file=f)
        for i in range(surface):
            print('x'+str(i),'p'+str(i),'xx'+str(i),'pp'+str(i),'xp'+str(i),'pr'+str(i),sep='\t',end='\t',file=f)
        print(file=f)
        for i in range(0,lt,step):
            print(t[i],end='\t',file=f)
            for k in range(surface):
                for j in range(nqhd):
                    print(qhdformer[j,i,k],end='\t',file=f)
            print(file=f)
    with open('QHDlatter.txt','w') as f:
        print('t',end='\t',file=f)
        for i in range(surface):
            print('x'+str(i),'p'+str(i),'xx'+str(i),'pp'+str(i),'xp'+str(i),'pr'+str(i),sep='\t',end='\t',file=f)
        print(file=f)
        for i in range(0,lt,step):
            print(t[i],end='\t',file=f)
            for k in range(surface):
                for j in range(nqhd):
                    print(qhdlatter[j,i,k],end='\t',file=f)
            print(file=f)
    with open('QHD.txt','w') as f:
        print('t',end='\t',file=f)
        for i in range(surface):
            print('x'+str(i),'p'+str(i),'xx'+str(i),'pp'+str(i),'xp'+str(i),'pr'+str(i),sep='\t',end='\t',file=f)
        print(file=f)
        for i in range(0,lt,step):
            print(t[i],end='\t',file=f)
            for k in range(surface):
                for j in range(nqhd):
                    print(qhdtotal[j,i,k],end='\t',file=f)
            print(file=f)
elif(jobtype=='pRepresentation'):
    x=Get_Grid('k.DVR')
    t=Get_Grid('t.DVR')
    psy=Get_psy('Phi.DVR',NGrid,lt,surface)
    Animate_Density(x[0],x[NGrid-1],x,t,psy,surface)
elif(jobtype=='WignerRepresentation'):
    print('in coming ...')
else:
    print('Unkown job type?')