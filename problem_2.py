from mpl_toolkits import mplot3d
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from math import *
import time

"""
最大速度: 79.175
"""

def warm_area(x):
    """回焊炉中温度分布"""
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

def heat_trans(k,beta_lst,v,train=True):
    """热传导数值解,k为方程系数"""
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
        
        a[t][0]=(1/(1+beta*dx))*(beta*dx*warm_area(dt*t*v/60)+a[t][1])
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




if __name__=='__main__':
    data = pd.read_excel(r"附件.xlsx").values
    T=[182,203,237,254,25]

    v_up,v_down=[100,65]

    while(abs(v_down-v_up)>1e-5):
        tem_lst,dt=heat_trans(0.00072,[1509,25046,139716,80178,232],(v_down+v_up)/2)
        if check(tem_lst,dt):
            v_down=(v_down+v_up)/2
        else:
            v_up=(v_down+v_up)/2
        print(v_down,' ',v_up)
    print('最大速度: %.3f'%(v_up))
    tem_lst,dt=heat_trans(0.00072,[1509,25046,139716,80178,232],(v_down+v_up)/2,train=False)
    plt.show()
    
