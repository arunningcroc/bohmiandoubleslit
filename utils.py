import numpy as np
import fastpotential as ts
#Find the nearest index in an array corresponding to a value.
#Used to convert particle position to grid point (i,j)
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
#Self explanatory
def x_deriv(psi,i,j,dx):
    return (psi[i+1,j]-psi[i,j])/dx
def y_deriv(psi,i,j,dx):
    return (psi[i,j+1]-psi[i,j])/dx

def create_wall(Npoints,x,y):
    V = np.zeros((Npoints,Npoints))
    Vdrop = np.zeros((Npoints,Npoints))
    sigmapot = 1/250
    factorpot = 50000
    yloc = np.linspace(-1,-0.13,40)
    yloc2 = np.linspace(-0.01,0.01,10)
    yloc3 = np.linspace(0.13,1,40)
    xloc = 0.3
    xdroploc = 0.9
    droploc = np.linspace(-1,1,100)
    for i in range(len(yloc)):
        V = V+ ts.create_gaussian(Npoints,x,y,xloc,yloc[i],factorpot,sigmapot)
    
    for i in range(len(yloc2)):
        V = V+ ts.create_gaussian(Npoints,x,y,xloc,yloc2[i],factorpot,sigmapot)
    for i in range(len(yloc3)):
        V = V+ ts.create_gaussian(Npoints,x,y,xloc,yloc3[i],factorpot,sigmapot)
    
    for i in range(len(droploc)):
        Vdrop = Vdrop + ts.create_gaussian(Npoints,x,y,xdroploc,droploc[i],-factorpot/10,sigmapot)
    return V+Vdrop