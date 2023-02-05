import numpy
import math 
import matplotlib.pyplot as plt 
from scipy.special import sph_harm 
from matplotlib.widgets import Slider, Button
import scipy.special
from scipy.special import assoc_laguerre
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap
from skimage import measure
import skimage

#Definir la función de onda en terminos de n,l,m y x,y,z 

def atmhidrogeno(n, l, m, X, Y, Z):
	R = numpy.sqrt(X**2 + Y**2 + Z**2)
	theta = numpy.arccos(Z/R)
	phi = numpy.arctan2(Y, X)

	rho = 2.*R/n 
	s_harm = sph_harm(m, l, phi, theta)
	l_poly = scipy.special.genlaguerre(n-l-1, 2*l+1,)(rho)

	prefactor = numpy.sqrt((2./n)**3*math.factorial(n-l-1)/2.*n*math.factorial(n+1))
	wf = prefactor*numpy.exp(-rho/2.)*rho**1*s_harm*l_poly
	wf = numpy.nan_to_num(wf)
	return wf 

#Definir grilla de graficación 3D

dz = 0.5 
zmin = -10
zmax = 10
x = numpy.arange(zmin, zmax, dz)
y = numpy.arange(zmin, zmax, dz)
z = numpy.arange(zmin, zmax, dz)
X, Y, Z = numpy.meshgrid(x, y, z)

#Elegir los valores de los numeros cuanticos 

n = 3
l = 2
m = 0

data = atmhidrogeno(n, l, m, X, Y, Z)
data = abs(data)**2 
R = numpy.sqrt(X**2 + Y**2 + Z**2)

#Grafica de la densidad de probabilades de la funcion de onda en en palano XZ y Y variable 

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.15, bottom=0.15)
im = plt.imshow(data[int((0-zmin)/dz),:,:], vmin=0, vmax=numpy.max(data), extent=[zmin,zmax,zmin,zmax])
plt.colorbar()
sli = Slider(plt.axes([0.25, 0.025, 0.65, 0.03]), "Y", z[0], z[len(z)-1], valinit=0)
ax.set_title("Orbital del hidrogeno xz"+str("%.2f"%sli.val)+"): n="+str(n)+", l="+str(l)+", m="+str(m))

def update(val):
	index = int((sli.val-zmin)/ dz)
	im.set_data(data[index,:,:])
	ax.set_title("Orbital del hidrogeno xz"+str("%.2f"%sli.val)+"): n="+str(n)+", l="+str(l)+", m="+str(m))

sli.on_changed(update)
plt.show()

#Grafica de los orbitales en 3D

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.set_xlim([0,len(x)])
ax.set_ylim([0,len(y)])
ax.set_zlim([0,len(z)])
max_val = numpy.max(data)

verts, faces, normals, values = measure.marching_cubes(data, max_val/2, spacing = (1,1,1))
result = ax.plot_trisurf(verts[:,0], verts[:,1], faces, verts[:,2], cmap ='Spectral', lw=0)

sli = Slider(plt.axes([0.25, 0.025, 0.65, 0.03]), "iso", 0, max_val, valinit = max_val/2)
ax.set_title("Orbitales atomo de hidrogeno("+str("%.5f"%sli.val)+"): n="+str(n)+", l="+str(l)+", m="+str(m))

def update (val):
	ax.clear()
	verts, faces, normals, values = measure.marching_cubes(data, sli.val, spacing = (1,1,1))
	result = ax.plot_trisurf(verts[:,0], verts[:,1], faces, verts[:,2], cmap ='Spectral', lw=0)
	ax.set_xlim([0,len(x)])
	ax.set_ylim([0,len(y)])
	ax.set_zlim([0,len(z)])
	ax.set_title("Orbitales atomo de hidrogeno("+str("%.5f"%sli.val)+"): n="+str(n)+", l="+str(l)+", m="+str(m))

sli.on_changed(update)
plt.show()