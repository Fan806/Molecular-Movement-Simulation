#-------------------------------------------------------------------------------
#-*-coding:cp936-*-
# Name:        casual
# Purpose:
#
# Author:      mac
#
# Created:     08/12/2016
# Copyright:   (c) mac 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from random import randint, uniform
from math import sqrt,acos,asin,pi,sin,cos
from tkinter import*
import matplotlib
matplotlib.use('TkAgg')
from PIL import ImageTk,Image
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.spatial.distance import pdist,squareform
import scipy.integrate as integrate
import matplotlib.animation as animation
font={'family':'fantasy','color':'Snow','weight':'normal','size':15}
font1={'family':'fantasy','color':'Snow','weight':'normal','size':10}
font3={'family':'fantasy','color':'Ivory','weight':'normal','size':15}
font4={'family':'fantasy','color':'Ivory','weight':'normal','size':10}
font2=('Comic Sans MS',11,'normal')
root=Tk()
root.title('Simulation of A Physical System')
image1=Image.open('4.jpg')
background_image=ImageTk.PhotoImage(image1)
Label(root,image=background_image).place(x=0,y=0,relwidth=1,relheight=1)
root.geometry('1700x1000')

intro_text="""
In the percolation project, grids were generated according to a specified probability distribution to experimentally determine the percolation probability. This approach is an example of a Monte Carlo method, an algorithm that relies on repeated random sampling to compute results. Monte Carlo methods are often used for simulating physical and mathematical systems when other methods are impossible or computationally infeasible."""
notice_text="""
Please be cautious when you want to change the initial state, because we will clear the graph generated before.
You can use this menu just by a click of the right mouse button.
If the N you give to us is greater than upper bound or less than lower bound, then N will be set as the lower bound."""
help_text="""
Execute the following steps, then you will get four pictures which will show you the analysis of the result:
1. Choose the initial state. If you give all the energy to system, the particles will move at first. If you give all the energy to demon, the particles will not move at first.
2. Input the total energy and steps.
3. Change the number of particles.
Click the button 'Virtual Laboratry', which will show you an animation of the simution of molecular motion.
The last button 'System Energy' is prepared for you to check all the data used in the simulation."""
info_text="""
Group Name:    Keep Going
Author:             Xie Fan
Time:               2016.12.30"""
def introduction():
    win1=Tk()
    win1.title('Introduction')
    Message(win1,text=intro_text,aspect=400).pack()
def notice():
    win2=Tk()
    win2.title('Notice')
    Message(win2,text=notice_text,aspect=400).pack()
def for_help():
    win3=Tk()
    win3.title('Help')
    Message(win3,text=help_text,aspect=400).pack()
def information():
    win4=Tk()
    win4.title('Information')
    win4.geometry('250x80')
    Message(win4,text=info_text,aspect=400).pack()
menubar=Menu(root)
filemenu=Menu(menubar,tearoff=0)
filemenu.add_command(label='introduction',command=introduction)
filemenu.add_separator()
filemenu.add_command(label='help',command=for_help)
filemenu.add_command(label='notice',command=notice)
filemenu.add_separator()
filemenu.add_command(label='information',command=information)
menubar.add_cascade(label='Menu',menu=filemenu)
root['menu']=menubar
def popup(event):
    menubar.post(event.x_root,event.y_root)
root.bind('<Button-3>',popup)

f1=Figure(figsize=(6.2,4.5),dpi=95,facecolor=(0,0,0),edgecolor=(0,0,0))
f1.patch.set_alpha(0.7)
pic1=f1.add_subplot(111)

f2=Figure(figsize=(6.2,4.5),dpi=95,facecolor=(0,0,0))
f2.patch.set_alpha(0.7)
pic2=f2.add_subplot(111)

f3=Figure(figsize=(6.2,4.5),dpi=95,facecolor=(0,0,0))
f3.patch.set_alpha(0.7)
pic3=f3.add_subplot(111)

f4=Figure(figsize=(6.2,4.5),dpi=95,facecolor=(0,0,0))
f4.patch.set_alpha(0.7)
pic4=f4.add_subplot(111)

global a,b,c,d,f
a=[]
b=[]
c=[]
d=[]
f=[]
def ideal_gas(
        N,                              # Number of particles
        totalenergy,                    # total of demon and system energy
        steps,                          # number of simulation steps
        state,                          # Initial state
        visuals=True                    # plot histogram?
    ):
    global a,b,c,f,g
    v0=sqrt(2.0*totalenergy/N)
    deltaV=v0/10
    k=0
    demono=[]
    x=[]
    g=[]
    if state==1:
        v=[v0]*N
        systemenergy=float(totalenergy)
        demonenergy=0.0
    elif state==2:
        v=[0.0]*N
        systemenergy=0.0
        demonenergy=float(totalenergy)
    else:
        Label(root,text='Please choose an initial state!',bg='#0e0906',fg='Silver',height=3,width=30,font=font2).place(x=10,y=730)
        Label(root,bg='#0e0906',height=2,width=30).place(x=10,y=770)
        return
    pic4.clear()
    for k in range(steps):
        for i in range(N):
            dv=uniform(-deltaV,deltaV)                            #a random velocity change(temporary)
            rani=randint(0,N-1)                                   #choose a random number
            Vnew=v[rani]+dv                                       #assume the v[rani] is changed
            deltaenergy=0.5*(Vnew**2-v[rani]**2)                  #calculate the delta energy
            if deltaenergy<demonenergy:                           #under the condition (delta energy<demon energy):
                v[rani]=Vnew                                      #the assuption of the change of volecity make effact
                systemenergy=systemenergy+deltaenergy             #calculation of the new system energy
                demonenergy=demonenergy-deltaenergy               #calculation of the demon energy
                visuals=True
        g.append(systemenergy)
        x.append(k+1)
        demono.append(demonenergy)                               #after one step, store the domon energy
    systemenergy=sum(g)/len(g)
    pic4.plot(x,demono,'y-o',color='#F8CC53',alpha=1)
    pic4.set_title("The Step Demon Energies Over Time",fontdict=font3)
    pic4.set_xlabel("steps",fontdict=font4)
    pic4.set_ylabel("demon energies",fontdict=font4)
    canvas4.show()
    a.append(N)
    b.append(systemenergy)
    Final_Molecule_Velocity(v)
    Demon_Energy(demono)
    e=[]
    e.append(N)
    e.append(state)
    e.append(totalenergy)
    e.append(steps)
    e.append(systemenergy)
    c.append(e)
    f.append(totalenergy)

def Final_Molecule_Velocity(v):
    pic1.clear()
    pic1.set_xlabel('v',fontdict=font1)
    pic1.set_ylabel('frequency',fontdict=font1)
    pic1.set_title("Final Particle Velocity Distribution",fontdict=font)
    v=np.asarray(v)
    n,bins,patches=pic1.hist(v,bins=10,facecolor='#F8CC53',edgecolor='#F7A303')
    canvas1.draw()

def Demon_Energy(demono):
    pic2.clear()
    pic2.set_xlabel('Demon Energy',fontdict=font4)
    pic2.set_ylabel('frequency',fontdict=font4)
    demono=np.asarray(demono)
    n,bins,patches=pic2.hist(demono,bins=10,facecolor='#F8CC53',edgecolor='#F7A303')
    pic2.set_title("Demon Energy Histogram",fontdict=font3)
    canvas2.draw()

def on_click():
    Label(root,bg='#0e0906',height=3,width=40).place(x=10,y=730)
    Label(root,bg='#0e0906',height=3,width=40).place(x=10,y=770)
    state=state_text.get()
    N=N_text.get()
    totalenergy=energy_text.get()
    steps=steps_text.get()
    global a,b,d,f
    if d==[]:
        d.append(state)
    elif d==[0]:
        d=[]
        d.append(state)
        if d[0]==0:
            N_text.set(0)
        else:
            N_text.set(low_b.get())
            N=N_text.get()
            a=[]
        Label(root,text="You've give us the initial state.",bg='#0e0906',fg='Silver',height=3,width=30,font=font2).place(x=10,y=730)
        Label(root,text="Let's start a simulated collision!",bg='#0e0906',fg='Silver',height=2,width=30,font=font2).place(x=10,y=770)
    else:
        if state not in d:
            Label(root,text="You've changed the state.",bg='#0e0906',fg='Silver',height=3,width=30,font=font2).place(x=10,y=730)
            Label(root,text="Let's start a new simulation!",bg='#0e0906',fg='Silver',height=2,width=30,font=font2).place(x=10,y=770)
            N_text.set(low_b.get())
            N=N_text.get()
            d=[]
            b=[]
            a=[]
            f=[]
    if len(a)>=1:
        if N<=a[len(a)-1]:
            win5=Tk()
            Label(win5,text='You give us a number which is less or equal to the previous one.',font=font2).place(x=5,y=10)
            Label(win5,text='Do you want to plot a new graph?',font=font2).place(x=5,y=45)
            win5.geometry('520x130')
            def yes():
                global a,b,f
                a=[]
                b=[]
                f=[]
                ideal_gas(N,totalenergy,steps,state,visuals=True)
                pic3.clear()
                pic3.plot(a,b,'c-o',color='#668B8B',label='System Energy',linewidth=0.8)
                pic3.plot(a,f,'m-o',color='#7FFFD4',label='Total Energy',linewidth=3.5)
                pic3.legend(loc=1,fontsize=10)
                pic3.set_title("N vs. Energy",fontdict=font)
                pic3.set_xlabel('N',fontdict=font1)
                pic3.set_ylabel('Energy',fontdict=font1)
                pic3.grid(True)
                canvas3.draw()
                win5.destroy()
            def no():
                N_text.set(a[len(a)-1])
                win5.destroy()
            Button(win5,text='Yes',command=yes,font=font2).place(x=180,y=85)
            Button(win5,text='No',command=no,font=font2).place(x=245,y=85)
        else:
            ideal_gas(N,totalenergy,steps,state,visuals=True)
    else:
        ideal_gas(N,totalenergy,steps,state,visuals=True)
    pic3.clear()
    pic3.plot(a,b,'c-o',color='#F8CC53',label='System Energy',linewidth=0.8)
    pic3.plot(a,f,'m-o',color='#7FFFD4',label='Total Energy',linewidth=3.5)
    pic3.legend(loc=1,fontsize=10)
    pic3.set_title("N vs. Energy",fontdict=font)
    pic3.set_xlabel('N',fontdict=font1)
    pic3.set_ylabel('Energy',fontdict=font1)
    pic3.grid(True)
    canvas3.draw()

def show_systemenergy():
    win=Tk()
    win.title('System Energy')
    for i in range(len(c)):
        Label(win,text='N='+str(c[i][0])+';   '+'state='+str(c[i][1])+';   '+'total energy='+str(c[i][2])+';   '+'steps='+str(c[i][3])+'systemenergy='+str(c[i][4]),bg='#0e0906',fg='Silver',compound='center',width=90,font=font2).pack(anchor=NW)

class particlebox:
    def __init__(self,init_state,init_d,bounds=[-2.7,2.7,-2.2,2.2],size=0.04):
        self.init_state=init_state
        self.init_d=init_d
        self.size=size
        self.state=self.init_state.copy()
        self.demon=self.init_d.copy()
        self.time_elapsed=0
        self.bounds=bounds
    def step(self,dt,deltaV):
        #step once by dt seconds
        self.time_elapsed+=dt

        #update positions
        self.state[:,:2]+=dt*self.state[:,2:]
        self.demon[:2]+=dt*self.demon[2:]

        #find pairs of particles undergoing a collision
        d=squareform(pdist(self.state[:,:2]))
        ind1,ind2=np.where(d<=2*self.size)
        unique=(ind1<ind2)
        ind1=ind1[unique]
        ind2=ind2[unique]
        v_d2=sqrt((self.demon[2])**2+(self.demon[3])**2)
        def ran():
            global p
            p=randint(-1,1)
            if p==0:
                ran()
            return p
        global dv_x,dv_y,dv
        #update velocities of colliding pairs
        if 0.5*v_d2**2!=energy_text.get():
          for i1,i2 in zip(ind1,ind2):
            r1=self.state[i1,:2]
            r2=self.state[i2,:2]
            v1=self.state[i1,2:]
            v2=self.state[i2,2:]
            r_r=r1-r2
            r_x=r_r[0]
            r_y=r_r[1]
            r_=sqrt(r_x**2+r_y**2)
            sinseta=r_y/r_
            cosseta=r_x/r_
            v1=self.state[i1,2:]
            v2=self.state[i2,2:]
            dv1=uniform(-deltaV,deltaV)
            dv2=uniform(-deltaV,deltaV)
            dv=dv1+dv2
            if v_d2>=dv:
             if v_d2>0:
                dv_x=(dv/v_d2)*self.demon[2]
                dv_y=(dv/v_d2)*self.demon[3]
             else:
                dv_x=uniform(-dv,dv)
                dv_y=sqrt(dv**2-dv_x**2)*ran()
             if v2[0]*cosseta+v2[1]*sinseta>v1[0]*cosseta+v1[1]*sinseta:
                v1_n=v2[0]*cosseta+v2[1]*sinseta-dv1
                v1_t=v1[1]*cosseta-v1[0]*sinseta
                v2_n=v1[0]*cosseta+v1[1]*sinseta-dv2
                v2_t=v2[1]*cosseta-v2[0]*sinseta
                v1[0]=v1_n*cosseta-v1_t*sinseta
                v1[1]=v1_n*sinseta+v1_t*cosseta
                v2[0]=v2_n*cosseta-v2_t*sinseta
                v2[1]=v2_n*sinseta+v2_t*cosseta
                self.state[i1,2:]=v1
                self.state[i2,2:]=v2
                self.demon[2]-=dv_x
                self.demon[3]-=dv_y
        else:
            count=[]
            for i in range(N_text.get()):
                dv=uniform(0,deltaV)
                count.append(dv)
                self.state[i,2:][0]=uniform(-dv,dv)
                self.state[i,2:][1]=sqrt(dv**2-self.state[i,2:][0]**2)
            dv=sum(count)
            dv_x=uniform(-dv,dv)
            dv_y=sqrt(dv**2-dv_x**2)*ran()
            self.demon[2]-=dv_x
            self.demon[3]-=dv_y

        #check for crossing boundary
        crossed_x1=(self.state[:,0]<self.bounds[0]+self.size)
        crossed_x2=(self.state[:,0]>self.bounds[1]-self.size)
        crossed_y1=(self.state[:,1]<self.bounds[2]+self.size)
        crossed_y2=(self.state[:,1]>self.bounds[3]-self.size)

        self.state[crossed_x1,0]=self.bounds[0]+self.size
        self.state[crossed_x2,0]=self.bounds[1]-self.size
        self.state[crossed_y1,1]=self.bounds[2]+self.size
        self.state[crossed_y2,1]=self.bounds[3]-self.size

        self.state[crossed_x1|crossed_x2,2]*=-1
        self.state[crossed_y1|crossed_y2,3]*=-1

        c1=(self.demon[0]<self.bounds[0]+self.size)
        c2=(self.demon[0]>self.bounds[1]-self.size)
        c3=(self.demon[1]<self.bounds[2]+self.size)
        c4=(self.demon[1]>self.bounds[3]-self.size)

        if c1:
            self.demon[0]=self.bounds[0]+self.size
        if c2:
            self.demon[0]=self.bounds[1]-self.size
        if c3:
            self.demon[1]=self.bounds[2]+self.size
        if c4:
            self.demon[1]=self.bounds[3]-self.size

        if c1 or c2:
            self.demon[2]*=-1
        if c3 or c4:
            self.demon[3]*=-1

def Lab():
    N=N_text.get()
    if N!=0:
        state=state_text.get()
        totalenergy=energy_text.get()
        steps=steps_text.get()
        v0=sqrt(2.0*totalenergy/N)
        global rect,box,dt,deltaV,fig,ax
        deltaV=v0/10
        np.random.seed(0)
        axi=np.random.uniform(-2,2,(N,2))
        d_axis=np.random.uniform(-2,2,size=2)
        if state==1:
            v1=[v0,0]*N
            demonenergy=0.0
        elif state==2:
            v1=[0,0]*N
            demonenergy=float(totalenergy)
        v1=np.asarray(v1,dtype=float).reshape(N,2)
        init_state=np.hstack((axi,v1))
        v_d=sqrt(2*demonenergy)
        v_d=[v_d,0]
        v_d=np.asarray(v_d,dtype=float)
        init_d=np.hstack((d_axis,v_d))
        box=particlebox(init_state,init_d)
        dt=1./fps.get()
        fig=plt.figure()
        fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
        ax=fig.add_subplot(111,aspect='equal',autoscale_on=False,xlim=(-3.2,3.2),ylim=(-2.4,2.4))
        particles1,=ax.plot([],[],'co',color='#F8CC53',ms=6)
        particles2,=ax.plot([],[],'yo',color='#6A5ACD',ms=6)
        rect=plt.Rectangle(box.bounds[::2],box.bounds[1]-box.bounds[0],box.bounds[3]-box.bounds[2],ec='none',lw=2,fc='none')
        ax.add_patch(rect)
        def init():
            """initialize animation"""
            global box,rect
            particles1.set_data([],[])
            particles2.set_data([],[])
            rect.set_edgecolor('#0e0906')
            return particles1,particles2,rect
        def animate(i):
            """perform animation step"""
            global box,rect,dt,ax,fig
            box.step(dt,deltaV)
            ms=int(fig.dpi*2*box.size*fig.get_figwidth()/np.diff(ax.get_xbound())[0])
            #update pieces of the animation
            rect.set_edgecolor('#0e0906')
            particles1.set_data(box.state[:,0],box.state[:,1])
            particles1.set_markersize(ms)
            particles2.set_data(box.demon[0],box.demon[1])
            particles2.set_markersize(ms)
            return particles1,particles2,rect
        ani=animation.FuncAnimation(fig,animate,frames=600,interval=10,blit=False,init_func=init)
        plt.show()
    else:
        Label(root,text="Please don't let N be zero!",bg='#0e0906',fg='Silver',font=font2,height=3,width=25).place(x=10,y=730)
        Label(root,bg='#0e0906',height=2,width=30).place(x=10,y=770)

def ok():
    global N_text,fps
    if low_b.get() < up_b.get():
        Label(root, text='Please choose the number of particles:', bg='#0e0906', fg='Silver', font=font2).place(x=10,y=500,anchor=NW)
        N_text = IntVar()
        Spinbox(root, relief=SUNKEN, from_=low_b.get(), to_=up_b.get(), increment=incre.get(), width=5,textvariable=N_text, command=on_click, bg='#0e0906', fg='Silver', font=font2).place(x=300, y=500)
        N_text.set(low_b.get()-incre.get())
        Label(root, text="Please give us the frames per second:", bg='#0e0906', fg='Silver', font=font2).place(x=10,
                                                                                                               y=560)
        fps = IntVar()
        Entry(root, relief=SUNKEN, width=8, textvariable=fps, bg='#0e0906', fg='Silver', font=font2).place(x=300, y=560)
        fps.set(100)

        Button(root, text="Simulated Collision", command=Lab, bg='#0e0906', fg='Silver', font=font2).place(x=10, y=610)
        Button(root, text="System Energy", command=show_systemenergy, bg='#0e0906', fg='Silver', font=font2).place(x=10,
                                                                                                                   y=670)
    else:
        win7 = Tk()
        Label(win7,
              text="The lower bound is greater or equal to the upper bound. Please choose the lower bound and the upper bound again.").pack()
        Button(win7, text="OK", command=win7.destroy).pack()
        win7.mainloop()


canvas1=FigureCanvasTkAgg(f1,root)
canvas1.get_tk_widget().place(x=410,y=30)

canvas2=FigureCanvasTkAgg(f2,root)
canvas2.get_tk_widget().place(x=1070,y=30)

canvas3=FigureCanvasTkAgg(f3,root)
canvas3.get_tk_widget().place(x=1070,y=520)

canvas4=FigureCanvasTkAgg(f4,root)
canvas4.get_tk_widget().place(x=410,y=520)

Label(root,text="Please choose the initial state:",bg='#0e0906',fg='Silver',font=font2).place(x=10,y=50,anchor=NW)
state_text=IntVar()
Radiobutton(root,variable=state_text,text='Give all the energy to system',value=1,bg='#0e0906',fg='Silver',font=font2).place(x=10,y=100)
Radiobutton(root,variable=state_text,text='Give all the energy to demon',value=2,bg='#0e0906',fg='Silver',font=font2).place(x=10,y=150)

Label(root,text='Please input the total energy:',bg='#0e0906',fg='Silver',font=font2).place(x=10,y=200,anchor=NW)
energy_text=IntVar()
Entry(root,relief=SUNKEN,width=8,textvariable=energy_text,bg='#0e0906',fg='Silver',font=font2).place(x=300,y=200)
energy_text.set(500)

Label(root,text='Please input the steps:',bg='#0e0906',fg='Silver',font=font2).place(x=10,y=250,anchor=NW)
steps_text=IntVar()
Entry(root,relief=SUNKEN,width=8,textvariable=steps_text,bg='#0e0906',fg='Silver',font=font2).place(x=300,y=250)
steps_text.set(3000)

Label(root,text="Please choose a lower bound:",bg='#0e0906',fg='Silver',font=font2).place(x=10,y=300)
low_b=IntVar()
Entry(root,relief=SUNKEN,width=8,textvariable=low_b,bg='#0e0906',fg='Silver',font=font2).place(x=300,y=300)
low_b.set(50)
Label(root,text="Please choose an upper bound:",bg='#0e0906',fg='Silver',font=font2).place(x=10,y=350)
up_b=IntVar()
Entry(root,relief=SUNKEN,width=8,textvariable=up_b,bg='#0e0906',fg='Silver',font=font2).place(x=300,y=350)
up_b.set(500)
Label(root,text="Please choose the increment:",bg='#0e0906',fg='Silver',font=font2).place(x=10,y=400)
incre=IntVar()
Entry(root,relief=SUNKEN,width=8,textvariable=incre,bg='#0e0906',fg='Silver',font=font2).place(x=300,y=400)
incre.set(50)
Button(root,text="OK",command=ok,bg='#0e0906',fg='Silver',font=font2).place(x=10,y=450)

Label(root,text="Hello!",bg='#0e0906',fg='Silver',height=3,width=30,font=font2).place(x=10,y=730)
Label(root,text="Welcome to our world!",bg='#0e0906',fg='Silver',height=3,width=30,font=font2).place(x=10,y=770)
root.mainloop()
