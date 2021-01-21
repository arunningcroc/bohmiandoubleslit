import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import fastpotential as ts
from utils import find_nearest,x_deriv, y_deriv, create_wall
#If you have ffmpeg, you can use this and the two lines at the bottom to save an animation
#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
 
# Update psi using a single RK4 step with timestep dt.

def update_psi_rk4(psi, dt):
    k1 = -(1. / (2.j)) * ts.laplace(psi,Npoints,dx) + (1. / 1.j)*V*(psi)
    k2 = -(1. / (2.j)) * ts.laplace(psi + 0.5 * dt * k1,Npoints,dx)+ (1. / 1.j)*V*(psi + 0.5 * dt * k1)
    k3 = -(1. / (2.j)) * ts.laplace(psi + 0.5 * dt * k2,Npoints,dx)+ (1. / 1.j)*V*(psi + 0.5 * dt * k2)
    k4 = -(1. / (2.j)) * ts.laplace(psi + dt * k3,Npoints,dx)+ (1. / 1.j)*V*(psi + dt * k3)
     
    psinew = psi + dt * 1. / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return psinew

def update_pos_rk4(pos,psi,x,y,dt,dx):
    #Runge-Kutta 4th order to update the position.
    #We find the indices closest to our double value (x,y).
    xcoord = find_nearest(x,pos[0])
    ycoord = find_nearest(y,pos[1])
    #Nan to num is used here to avoid potential division by zero.
    #This part updates the x coordinate.
    k1   = np.nan_to_num((x_deriv(psi,xcoord,ycoord,dx)/psi[xcoord,ycoord]).imag )
    k2   = np.nan_to_num((x_deriv(psi+0.5*dt*k1,xcoord,ycoord,dx)/(psi[xcoord,ycoord]+0.5*dt*k1)).imag)
    k3   = np.nan_to_num((x_deriv(psi+0.5*dt*k2,xcoord,ycoord,dx)/(psi[xcoord,ycoord]+0.5*dt*k2)).imag)
    k4   = np.nan_to_num((x_deriv(psi+dt*k3,xcoord,ycoord,dx)/(psi[xcoord,ycoord]+dt*k3)).imag)
    xnew = pos[0] + dt * 1. / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    #Updating the y part with RK4..
    k1   = np.nan_to_num((y_deriv(psi,xcoord,ycoord,dx)/psi[xcoord,ycoord]).imag )
    k2   = np.nan_to_num((y_deriv(psi+0.5*dt*k1,xcoord,ycoord,dx)/(psi[xcoord,ycoord]+0.5*dt*k1)).imag)
    k3   = np.nan_to_num((y_deriv(psi+0.5*dt*k2,xcoord,ycoord,dx)/(psi[xcoord,ycoord]+0.5*dt*k2)).imag)
    k4   = np.nan_to_num((y_deriv(psi+dt*k3,xcoord,ycoord,dx)/(psi[xcoord,ycoord]+dt*k3)).imag)
    ynew = pos[1] + dt * 1. / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return [xnew, ynew]

#Parameters here. We are using atomic units. 
L         = 1.
Npoints   = 201   
sigma     = 1./4.
x         = np.linspace(-L, L, Npoints)
y         = np.linspace(-L, L, Npoints)
dx        = x[1]-x[0]
time_unit = 2.4188843265857e-17
timestep  = 0.00001
psi       = np.zeros((Npoints,Npoints), dtype=complex)
kx        = -1000.0
V = create_wall(Npoints,x,y)
# Initialize and normalize psi.
for i in range(Npoints):
    for j in range(Npoints):
        psi[i,j] = np.exp(-((x[i]+0.2)**2+(y[j])**2) / 2 / sigma / sigma)*(np.cos(kx*x[i])+1j*np.sin(kx*x[i]))  
norm = ts.get_norm(psi,Npoints,dx)
#print(norm)
psi = psi/np.sqrt(norm)
#Initial particle position
pos = [-0.15,0.10]

# Set up figure.
fig, ax = plt.subplots()
line = ax.imshow(np.abs(psi)**2,cmap="gist_rainbow")

plt.xlabel(r'$y$')
plt.ylabel(r'$x$')
textheight = abs(np.max(psi))**2
plt.title(r'A Bohmian particle along with its wave function')
rk4_steps_per_frame = 4
 

#Animate everything

def animate(i):
    global psi
    global pos
    for q in range(rk4_steps_per_frame):
        psinew = update_psi_rk4(psi, timestep)
        posnew = update_pos_rk4(pos, psi, x,y,timestep, dx )
        psi = psinew
        pos = posnew

    ax.patches = []
    #Our particle!
    electron = plt.Circle((find_nearest(x,pos[1]),find_nearest(x,pos[0])),2,color="black")
    ax.add_patch(electron)
    currentnorm = ts.get_norm(psi,Npoints,dx)
    #If the norm changes from 1 significantly, the simulation is probably in trouble.
    print(i,currentnorm)
    plotted = abs(psi)**2 + V.real
    line.set_data(plotted)  # update the data
    line.set_clim(vmax=np.amax(np.abs(psi)**2))
    line.set_clim(vmin=0)
    return line

def init():
    return line

#show or animate

ani = animation.FuncAnimation(fig, animate, np.arange(1, 530), init_func=init,
                              interval=25, save_count=530)
#plt.show()
#FFwriter=animation.FFMpegWriter(fps=60, extra_args=['-vcodec', 'libx264'])
#ani.save('this.mp4', writer = FFwriter)

plt.show()
