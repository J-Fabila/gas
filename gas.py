import numpy as np
import random
from vpython import *
## # # # # # #  #  #  #   #   #     #    #    #    #   #  #  ##
#                                     Declaración de constantes
## # # # # # #  #  #  #   #   #     #    #    #    #   #  #  ##

N_particulas=100
k=1.38E-23   #Konstante de Boltzmann (K)
time_step=0.01
criterio_tiempo=1
tiempo=0
tiempo_total=1000
radios_atomicos = {'H' : 0.46 , 'O' : 0.74, 'N' : 0.74, 'S' : 1.04, 'C' : 0.77, 'P' : 1.10 } # A
masas_atomicas = {'H' : 1.1 , 'O' : 16 , 'N' : 14, 'S' : 32, 'C' : 12 , 'P' : 31 }  # uma
lennard_jones={'H' : [43,3.3] , 'O' : [162,3.9], 'N' : [131,4.2] } #lenard_jones={Simbolo: [E_Disociación [J], Radio_equilibrio [A]]}

## # # # # # #  #  #  #   #   #     #    #    #    #   #  #  ##
#                                         Declaración de clases
## # # # # # #  #  #  #   #   #     #    #    #    #   #  #  ##
"""
self.e = lennard_jones[Simbolo][0]              # Energía disociación [J]
self.r = lennard_jones[Simbolo][1]              # Radio de equilibrio [A]
"""

class Atomo:
    def __init__(self, _Simbolo, _x, _y, _z):
        self.Simbolo = _Simbolo                          # Simbolo atómico
        self.R = Radios[_Simbolo]               # Radio atómico [A]
        self.M =  Masas[_Simbolo] #*1.66054e-27    # Masa atómica [uma]
        self.x = _x  # Posición X [A]
        self.y = _y  # Posición Y [A]
        self.z = _z  # Posición Z [A]
#### ALTERNATIVA: self.x = np.zeros(3) self.x[0] == X ; self.x[1] == Y ;self.x[1] == Z
        self.vx = _v_x                                     # Velocidad X [A/s]
        self.vy = _v_y                                     # Velocidad X [A/s]
        self.vz = _v_z                                     # Velocidad X [A/s]
#### ALTERNATIVA: self.v = np.zeros(3) self.v[0] == X ; self.v[1] == Y ;self.v[1] == Z
        self.ax = 0                                     # Aceleración X [A/s^2]
        self.ay = 0                                     # Aceleración Y [A/s^2]
        self.az = 0                                     # Aceleración Z [A/s^2]
#### ALTERNATIVA: self.a = np.zeros(3) self.a[0] == X ; self.a[1] == Y ;self.a[1] == Z
    def __repr__(self):
        return self.Symbol+" "+str(self.x)+" "+str(self.y)+" "+str(self.z)+"\n"


class Estructura_Atomica:
    def __init__(self, N):
        self.Nat=N   # Número de átomos
        self.atomos=[]
    def Initializacion_Boltzmann(self, Temperatura, Celda):
        #Implementar inicializacion de velocidades
        #v=np.sqrt(np.sqrt(8.0*k*T/(np.pi*(self.m))))/3.0 # Velocidad de acuerdo con Boltzmann
        pass
    def Leer_XYZ(self, file):
        with open(file,"r") as f:
            data=f.readlines()
            contador=1
            for line in data:
                words = line.split()
                if(contador==1):
                    self.Nat=words[0]
                if(contador>2):
                    atom=Atomo(words[0],words[1],words[2],words[3])
                    self.atomos.append(atom)
                contador=contador+1
        print(self.atoms)

class celda:
    def __init__(self,xmin=-10.0,xmax=10.0,ymin=-10.0,ymax=10.0,zmin=-10.0,zmax=10.0):
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        self.zmin=zmin
        self.zmax=zmax

## # # # # # #  #  #  #   #   #     #    #    #    #   #  #  ##
#                                      Declaración de funciones
## # # # # # #  #  #  #   #   #     #    #    #    #   #  #  ##

def distancia_atomica(atom_1,atom_2):
    return np.sqrt( np.power(atom_1.vx-atom_2.vx,2)+np.power(atom_1.vy-atom_2.vy,2)+np.power(atom_1.vz-atom_2.vz,2))
## A lo mejor habria k vectorizarla

def Fuerza_l_j(atom_1,atom_2):
    r_sep=distancia_atomica(atom_1,atom_2)
    return 4*atom_1.e*( (12.0*np.power((atom_1.r),12)/np.power(r_sep,13)) - (6.0*np.power((atom_1.r),6)/np.power(r_sep,7)) )
## A lo mejor tambien habria k vectorizarla
"""''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Otra opción es implementar las fuerzas como vectores (devuelvan vectores)
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"""

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

caja=celda()
gas=Estructura_Atomica()
gas.Initializacion_Boltzmann()
Dinamica_molecular(gas,celda)
