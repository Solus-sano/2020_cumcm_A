import random
from mpl_toolkits import mplot3d
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import shuffle,randint
from math import *
import time
import copy
import operator
"""
温度:  [170.93892305,185.12134781,228.98774734,264.77149086,25]
速度:  88.42446571224875
用时:  945s
"""

def warm_area(T,x):
    """回焊炉中环境温度分布"""
    def liner(x1,y1,x2,y2,x):
        return y1+(y2-y1)*(x-x1)/(x2-x1)
    if x<25:
        return liner(0,25,25,T[0],x)
    if 25<=x and x<197.5:
        return T[0]
    if 197.5<=x and x<202.5:
        return liner(197.5,T[0],202.5,T[1],x)
    if 202.5<=x and x<=233:
        return T[1]
    if 233<=x and x<238:
        return liner(233,T[1],238,T[2],x)
    if 238<=x and x<268.5:
        return T[2]
    if 268.5<=x and x<273.5:
        return liner(268.5,T[2],273.5,T[3],x)
    if 273.5<=x and x<344.5:
        return T[3]
    if 344.5<=x and x< 349.5:
        return liner(344.5,T[3],349.5,T[4],x)
    if x>=349.5:
        return T[4]

def heat_trans(k,beta_lst,T,v,train=True):
    """
    热传导数值解,k为方程系数,beta_lst为电子版与环境间热交换系数
    T为温区设置,v为速度
    """
    m=20;n=10000;t0=400;x0=0.015

    dt=t0/(n-1);dx=x0/(m-1)
    a=np.zeros((n,m))
    for xx in range(0,m):
        a[0][xx]=25

    r=k**2*dt/dx**2
    s=1-2*r

    for t in range(1,n):

        real_t=t*dt
        real_x=real_t*v/60
        if real_x<197.5:
            beta=beta_lst[0]
        elif 197.5<=real_x and real_x<233:
            beta=beta_lst[1]
        elif 233<=real_x and real_x<268:
            beta=beta_lst[2]
        elif 268<=real_x and real_x<344.5:
            beta=beta_lst[3]
        else:
            beta=beta_lst[4]


        for xx in range(1,m-1):
            # a[t][xx]=a[(t-1)][xx]+(dt/(dx**2))*(a[(t-1)][xx+1]+a[(t-1)][xx-1]-2*a[(t-1)][xx])*(k**2)
            a[t][xx]=s*a[t-1][xx]+r*(a[t-1][xx-1]+a[t-1][xx+1])
        
        a[t][0]=(1/(1+beta*dx))*(beta*dx*warm_area(T,dt*t*v/60)+a[t][1])
        a[t][m-1]=a[t][0]

    if train==False:
        # for x in [111.25,217.75,253.25,304]:
        #     idx=int((x*60/v)/dt)
        #     print(a[idx][m//2])

        plt.figure()
        ax=plt.axes(projection='3d')
        # for idx,i in enumerate(np.linspace(0,x0,m)):
        #     ax.plot([i for j in range(n)],np.linspace(0,t0*7/6,n),a[:,idx])
        X,Y=np.meshgrid(np.linspace(0,x0,m),np.linspace(0,t,n))
        ax.plot_surface(X,Y,a,cmap='viridis')
        # plt.plot(list(np.linspace(0,t0*7/6,n)),a[:,m//2])
        # print(np.linspace(0,t0,n))
        plt.grid()

        plt.figure()
        plt.plot()
        t_init=data[:,0]
        y_init=data[:,1]
        x_init=t_init*7/6
        # plt.plot(t_init,y_init)
        plt.plot(np.linspace(0,t0,n),a[:,m//2])
        plt.plot(np.linspace(0,t0,n),a[:,0])
        plt.xlabel('t')
        plt.grid()
    return a[:,m//2],dt

def check(tem_lst,dt):
    """检查温度时间序列是否满足制程限制"""
    t=data[:,0].reshape((-1,))
    flag=1

    d_tem=np.diff(tem_lst)/dt
    d_tem=np.hstack((d_tem,d_tem[-1:]))
    # print(np.amax(d_tem),'   ',np.amin(d_tem))
    if np.amax(d_tem)>3 or np.amin(d_tem)<-3:
        flag=0
    if np.amax(tem_lst)>250 or np.amax(tem_lst)<240:
        flag=0

    t_down,t_up=[0,400]
    for idx,tem in enumerate(tem_lst):
        if tem<150 and d_tem[idx]>0:
            t_down=idx*dt
        if tem<190 and d_tem[idx]>0:
            t_up=idx*dt
    if t_up-t_down>120 or t_up-t_down<60:
        flag=0

    cnt=0
    for idx,tem in enumerate(tem_lst):
        if tem>217:
            cnt+=1
    if cnt*dt>90 or cnt*dt<40:
        flag=0

    return flag


def Loss(T,v):
    tem_lst,dt=heat_trans(0.00072,[1509,25046,139716,80178,232],T,v)

    d_tem=np.diff(tem_lst)/dt
    d_tem=np.hstack((d_tem,d_tem[-1:]))
    idx_down,idx_up=[0,0]
    for idx,tem in enumerate(tem_lst):
        if tem<217 and d_tem[idx]>0:
            idx_down=idx
        if d_tem[idx]>0:
            idx_up=idx
    idx_mid=np.argmax(tem_lst)

    S=0
    symmetry_loss=0
    for idx in range(idx_down,idx_up):
        S+=(tem_lst[idx]-217)*dt
        symmetry_loss+=(tem_lst[idx]-tem_lst[2*idx_mid-idx])**2
    
    # print(S,' ',symmetry_loss)
    return S+(1-check(tem_lst,dt))*inf+symmetry_loss*0.001

def heredity_atc(w,g,p):#(种群个数，进化代数,变异率)
    print('-------------------------')
    print('种群个数：%d'%(w))
    print('进化代数：%d'%(g))
    print('变异率：%f'%(p))
    print('-------------------------')
    
    
    """生成初始解"""
    J=[]
    S_lst=[]
    def rand_synthesis(idx):
        gen_tmp=np.empty((6,))
        gen_tmp[0]=165+random.random()*20
        gen_tmp[1]=185+random.random()*20
        gen_tmp[2]=225+random.random()*20
        gen_tmp[3]=245+random.random()*20
        gen_tmp[4]=65+random.random()*35
        return gen_tmp[idx]
    for i in range(w):
        gen_tmp=np.empty((5,))
        for j in range(5):
            gen_tmp[j]=rand_synthesis(j)
        J.append(gen_tmp)
    # print(J)
    """迭代"""
    for k in range(g):
        A=copy.deepcopy(J)
        c1=np.arange(w);shuffle(c1)#交叉操作的染色体配对组

        """交叉"""
        for i in range(0,w,2):
            idx_lst=list(set(randint(5,size=(5,))))
           
            l1=A[c1[i]]
            l2=A[c1[i+1]]
            for idx in idx_lst:
                tmp=l1[idx]
                l1[idx]=l2[idx]
                l2[idx]=tmp

            A[c1[i]]=l2;A[c1[i+1]]=l1
        
        """变异"""
        B=copy.deepcopy(J)
        for i in range(w):
            if random.random()>p: continue;

            l=B[i]
            idx_lst=list(set(randint(5,size=(5,))))
            for idx in idx_lst:
                l[idx]=rand_synthesis(idx)

            B[i]=l

        M=J+A+B
        """选择"""
        lst=[]
        for idx,item in enumerate(M):
            lst.append((idx,Loss(np.hstack((M[idx][:-1],np.array([25]))),M[idx][-1])))
        lst=sorted(lst,key=operator.itemgetter(1))

        J_tmp=[]

        for i in range(w):
            J_tmp.append(M[lst[i][0]])

        J=J_tmp
        # print(J)
        print("第%d次进化。"%(k+1))
        S_lst.append(lst[0][1])
    plt.figure()
    plt.plot(S_lst)
    return Loss(np.hstack((J[0][:-1],np.array([25]))),J[0][-1]),J[0]

if __name__=='__main__':
    data = pd.read_excel(r"附件.xlsx").values
    inf=9999999999999
    T=[182,203,237,254,25]

    t_begin=time.time()
    S,ans=heredity_atc(50,20,0.1)
    print('面积: %.3f'%(S))
    print("温度: ",ans[:-1])
    print("速度: ",ans[-1])
    print("用时: ",time.time()-t_begin)
    tem_lst,dt=heat_trans(0.00072,[1509,25046,139716,80178,232],np.hstack((ans[:-1],np.array([25]))),ans[-1],train=False)
    plt.show()
    
