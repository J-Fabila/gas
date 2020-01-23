import numpy as np
import random
from vpython import *

N_particulas=100
k=1.38E-23   #Konstante de Boltzmann (K)
time_step=0.01
criterio_tiempo=1
tiempo=0
tiempo_total=1000
radios_atomicos = {'H' : 0.46 , 'O' : 0.74, 'N' : 0.74, 'S' : 1.04, 'C' : 0.77, 'P' : 1.10 } # A
masas_atomicas = {'H' : 1.1 , 'O' : 16 , 'N' : 14, 'S' : 32, 'C' : 12 , 'P' : 31 }  # uma
lennard_jones={'H' : [43,3.3] , 'O' : [162,3.9], 'N' : [131,4.2] } #lenard_jones={Simbolo: [E_Disociación [J], Radio_equilibrio [A]]}

class celda:
    def __init__(self,xmin=-10.0,xmax=10.0,ymin=-10.0,ymax=10.0,zmin=-10.0,zmax=10.0):
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        self.zmin=zmin
        self.zmax=zmax

class atomo:
    def __init__(self, Simbolo, caja,T=300):
        self.simbolo = Simbolo                          # Simbolo atómico
        self.R = radios_atomicos[Simbolo]               # Radio atómico [A]
        self.m = masas_atomicas[Simbolo]*1.66054e-27    # Masa atómica [uma]
        self.e = lennard_jones[Simbolo][0]              # Energía disociación [J]
        self.r = lennard_jones[Simbolo][1]              # Radio de equilibrio [A]
        self.x = random.uniform(caja.xmin , caja.xmax)  # Posición X [A]
        self.y = random.uniform(caja.ymin , caja.ymax)  # Posición Y [A]
        self.z = random.uniform(caja.zmin , caja.zmax)  # Posición Z [A]

        v=np.sqrt(np.sqrt(8.0*k*T/(np.pi*(self.m))))/3.0 # Velocidad de acuerdo con Boltzmann

        self.vx = v                                     # Velocidad X [A/s]
        self.vy = v                                     # Velocidad X [A/s]
        self.vz = v                                     # Velocidad X [A/s]
        self.ax = 0                                     # Aceleración X [A/s^2]
        self.ay = 0                                     # Aceleración Y [A/s^2]
        self.az = 0                                     # Aceleración Z [A/s^2]

def distancia_atomica(atom_1,atom_2):
    return np.sqrt( np.power(atom_1.vx-atom_2.vx,2)+np.power(atom_1.vy-atom_2.vy,2)+np.power(atom_1.vz-atom_2.vz,2))

def l_j(atom_1,atom_2):
    r_sep=distancia_atomica(atom_1,atom_2)
    return 4*atom_1.e*( (12.0*np.power((atom_1.r),12)/np.power(r_sep,13)) - (6.0*np.power((atom_1.r),6)/np.power(r_sep,7)) )

caja=celda()

gas=[] #Lista vacia

for i in range(N_particulas):
    gas.append(atomo('H',caja))

while criterio_tiempo==0:
    for i in range(N_particulas):  ##Calcula las fuerzas para todas las particulas
        fx=0 ; fy=0 ; fz=0 ;
        for j in range(N_particulas):
            if (i!=j):
                F=l_j(gas[i],gas[j])
                fx=fx+F*(gas[i].x-gas[j].x)
                fy=fy+F*(gas[i].y-gas[j].y)
                fz=fz+F*(gas[i].z-gas[j].z)
            gas[i].ax=fx/gas[i].m
            gas[i].ay=fy/gas[i].m
            gas[i].az=fz/gas[i].m

    for i in range(N_particulas):  ## EULER
        gas[i].vx=gas[i].vx+(time_step*gas[i].ax)
        gas[i].vy=gas[i].vy+(time_step*gas[i].ay)
        gas[i].vz=gas[i].vz+(time_step*gas[i].az)

        gas[i].x=gas[i].x+(time_step*gas[i].vx)
        gas[i].y=gas[i].y+(time_step*gas[i].vy)
        gas[i].z=gas[i].z+(time_step*gas[i].vz)

    tiempo=tiempo+1;
    if(tiempo>tiempo_total):
        criterio_tiempo=1;
