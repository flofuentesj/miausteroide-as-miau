import numpy as np 
import matplotlib.pyplot as plt

#sin asteroides

def f(u, t , ms=1000.0): #mj = 1
    xs, ys, Vxs, Vys, xj, yj, Vxj, Vyj = u
    
    a=((xs-xj)**2 + (ys - yj)**2)**(1.5)
    dxsj = xs - xj
    dysj = ys - yj
    
    return  np.array([Vxs, Vys, -dxsj/a, -dysj/a, #ec sol
                      Vxj, Vyj, ms*dxsj/a, ms*dysj/a]) # ec jupi
        
tmax=1.0
h=0.00005
t=np.arange(0,tmax,h)

u= np.empty([t.size,8]) # 8 pq ecuaciones

#condiciones iniciales 8
u[0]=[0.0, 0.0, 0.0, 0.0, # condiciones iniciales para el sol (en el origen y con v_0 = 0)
      1.0, 0.0, 0.0, 31.785  # condiciones iniciales para el jupiter (al costado derecho del sol y con velocidad arriba valor incierto try and errror)
     ]


for n in range(t.size-1):
    K1= f(u[n],t[n])
    K2= f(u[n]+0.5*h*K1,t[n]+0.5*h)
    K3= f(u[n]+0.5*h*K2,t[n]+0.5*h)
    K4= f(u[n]+h*K3,t[n]+h)
    u[n+1] = u[n] + (h/6)*(K1 + 2*K2 + 2*K3 + K4)
    
#plt.plot(u)

plt.figure(figsize =(6,6))
xs = u[:, 0]
ys = u[:, 1]
xj = u[:, 4]
yj = u[:, 5]



#plt.scatter(xs,ys,color="orange",marker=".", label="Sol")
#plt.plot(xj,yj, color="tan", label="Jupiter" )
#plt.legend()
#plt.grid(True)
#plt.title("Gráfico de trayectoria entre el Sol y Jupiter ")
#plt.show()
#plt.savefig("Trayectoria_sj.pdf")


for n in range(100):
    plt.scatter(u[n,0],u[n,1] , label ="Sol")
    plt.scatter(u[n,4],u[n,5], label ="Jupiter")
    plt.title("Gráfico de trayectoria entre el Sol y Jupiter a distintos tiempos")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()
    plt.savefig("Trayectoria_sj_"+str(n)+".png")


    
