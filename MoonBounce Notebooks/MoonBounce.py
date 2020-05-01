#!/usr/bin/env python
# coding: utf-8

# In[2]:


from urllib.parse import urlencode
from urllib.request import urlretrieve
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, Angle
from IPython.display import Image
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, get_body,get_moon
import astropy
import pandas as pd
from astropy.coordinates import Angle, Latitude, Longitude
import geopandas
import numpy as np


# In[3]:


coordinates = EarthLocation.of_site('GBT')
print(coordinates)


# In[4]:


#Takes in length of obervation in hours, minute. 
#Returns an isot time format.
def timerange(day, hours, minutes):
    times = []
    nextday = str(day+1)
    currday = str(day)
    counter = 0
    for h in range(0, hours):
        if (h == hours):
            mins = minutes
        else:
            mins = 60
            for m in range(0, mins):
                  for s in range(0,60):
                    if ((22+h)%24 == 0) and (counter == 0):
                        currday = nextday
                        counter += 1
                    times = np.append(times,('2019-09-0'+currday+'T'+str((22+h)%24)+":"+str((20+m)%60)+":"+str(s)))
    return times
#Takes in an array of isot times.
#Returns an array of moon sky coordinates for times.
def getmooncords(times):
    timeslist = Time(times, format='isot', scale='utc')
    return get_moon(timeslist)

times = timerange(1, 1, 0)
mooncord = get_moon(Time(times, format='isot', scale='utc'), coordinates)
mooncord


# In[5]:


mooncord[0]


# In[6]:


ra = []
declination = []
distance = []
for i in range(len(mooncord)):
    declination = np.append(declination, mooncord[i].dec)
    ra = np.append(ra, mooncord[i].ra)
    distance = np.append(distance, mooncord[i].distance)
np.set_printoptions(suppress=True)
declination = declination.value
ra = ra.value
distance = distance.value


# In[52]:


import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
xs = np.arange(7)
ys = xs**2
fig = plt.figure(figsize=(5, 10))
ax = plt.subplot(2, 1, 1)

# If we want the same offset for each text instance,
# we only need to make one transform.  To get the
# transform argument to offset_copy, we need to make the axes
# first; the subplot command above is one way to do this.
trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, x=0.05, y=0.10, units='dots')
num =len(ra)

for x, y in zip(ra[:num], distance[:num]):
    plt.plot(x, y, 'ro')
    plt.text(x, y, '', transform=trans_offset)


# offset_copy works for polar plots also.
ax = plt.subplot(2, 1, 2, projection='polar')
trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, y=6, units='dots')

for x, y in zip(ra[:num], distance[:num]):
    plt.polar((x * np.pi)/180, y, 'ro')
    plt.text((x * np.pi)/180, y, '',
             transform=trans_offset,
             horizontalalignment='center',
             verticalalignment='bottom')

plt.show()


# In[50]:


xs = np.arange(7)
ys = xs**2

fig = plt.figure(figsize=(5, 10))
ax = plt.subplot(2, 1, 1)

# If we want the same offset for each text instance,
# we only need to make one transform.  To get the
# transform argument to offset_copy, we need to make the axes
# first; the subplot command above is one way to do this.
trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, x=0.05, y=0.10, units='dots')
num =1000

for x, y in zip(declination[:num], distance[:num]):

    plt.plot(x, y, 'ro')
    plt.text(x, y, '', transform=trans_offset)


# offset_copy works for polar plots also.
ax = plt.subplot(2, 1, 2, projection='polar')
trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, y=6, units='dots')

for x, y in zip(declination[:num], distance[:num]):
    plt.polar((x * np.pi)/180, y, 'ro')
    plt.text((x * np.pi)/180, y, '',
             transform=trans_offset,
             horizontalalignment='center',
             verticalalignment='bottom')
plt.show()


# In[55]:


import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import itertools

# The matrix for completness
A = np.matrix([[4,3,4], [3,1,3], [4,3,4]])
A = scipy.linalg.hessenberg(A)

num_points = 100
theta = np.linspace(0, 2 * np.pi, num_points)
phi = np.linspace(0, np.pi, num_points)

THETA, PHI = np.meshgrid(theta, phi)


X = np.sin(PHI) * np.cos(THETA)
Y = np.sin(PHI) * np.sin(THETA)
Z = np.cos(PHI)

# Calculate RQS for points on unit sphere:
RQS = np.empty(PHI.shape)
for theta_pos, phi_pos in itertools.product(range(num_points), range(num_points)):
    x = np.array([X[theta_pos, phi_pos],
                  Y[theta_pos, phi_pos],
                  Z[theta_pos, phi_pos]])
    RQS[theta_pos, phi_pos] = np.dot(x, np.dot(A, np.conj(x).T))

# normalize in range 0..1
maxRQS = abs(RQS).max()
N = (RQS+maxRQS)/(2*maxRQS)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,facecolors=cm.jet(N),linewidth=0, antialiased=False, shade=False)
plt.show()


# In[61]:


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 10 * np.outer(np.cos(u), np.sin(v))
y = 10 * np.outer(np.sin(u), np.sin(v))
z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='grey')


# In[47]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
num = 300
x = 10 * np.outer(np.cos(np.radians(declination[:num])), np.radians(np.sin(ra[:num])))
y = 10 * np.outer(np.sin(np.radians(declination[:num])), np.sin(np.radians(ra[:num])))
z = 10 * np.outer(np.ones(np.size(declination[:num])), np.cos(np.radians(ra[:num])))
                 
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='grey')


# In[48]:


function addTube(radius, segments, P1, P2) {
  // Q = P1→P2 moved to origin
  var Qx = P2[0] - P1[0];
  var Qy = P2[1] - P1[1];
  var Qz = P2[2] - P1[2];
  (/, Create, vectors, U, and, V, that, are, (1), mutually, perpendicular, and, (2), perpendicular, to, Q)
  if (Qx != 0) {  // create a perpendicular vector on the XY plane
    // there are an infinite number of potential vectors; arbitrarily select y = 1
    var Ux = -Qy/Qx;
    var Uy = 1;
    var Uz = 0;
    // to prove U is perpendicular:
    // (Qx, Qy, Qz)·(Ux, Uy, Uz) = Qx·Ux + Qy·Uy + Qz·Uz = Qx·-Qy/Qx + Qy·1 + Qz·0 = -Qy + Qy + 0 = 0
    }
  else if (Qy != 0) {  // create a perpendicular vector on the YZ plane
    var Ux = 0;
    var Uy = -Qz/Qy;
    var Uz = 1;
    }
  else {  // assume Qz != 0; create a perpendicular vector on the XZ plane
    var Ux = 1;
    var Uy = 0;
    var Uz = -Qx/Qz;
    }
  (/, The, cross, product, of, two, vectors, is, perpendicular, to, both,, so, to, find, V:)
  (/, (Vx,, Vy,, Vz), =, (Qx,, Qy,, Qz)×(Ux,, Uy,, Uz), =, (Qy×Uz, -, Qz×Uy,, Qz×Ux, -, Qx×Uz,, Qx×Uy, -, Qy×Ux))
  var Vx = Qy*Uz - Qz*Uy;
  var Vy = Qz*Ux - Qx*Uz;
  var Vz = Qx*Uy - Qy*Ux;
  (/, normalize, U, and, V:)
  var Ulength = Math.sqrt(Math.pow(Ux, 2) + Math.pow(Uy, 2) + Math.pow(Uz, 2));
  var Vlength = Math.sqrt(Math.pow(Vx, 2) + Math.pow(Vy, 2) + Math.pow(Vz, 2));
  Ux /= Ulength;
  Uy /= Ulength;
  Uz /= Ulength;
  Vx /= Vlength;
  Vy /= Vlength;
  Vz /= Vlength;
  for (var i = 0; i < segments; i++) {
    var θ = 2*Math.PI*i/segments;  // theta
    var dx = radius*(Math.cos(θ)*Ux + Math.sin(θ)*Vx);
    var dy = radius*(Math.cos(θ)*Uy + Math.sin(θ)*Vy);
    var dz = radius*(Math.cos(θ)*Uz + Math.sin(θ)*Vz);
    drawLine(P1[0] + dx, P1[1] + dy, P1[2] + dz,  // point on circle around P1
             P2[0] + dx, P2[1] + dy, P2[2] + dz)  // point on circle around P2 


# In[53]:


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# prepare the sphere surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
phi = np.linspace(0,2*np.pi, 50)
theta = np.linspace(0, np.pi, 25)
x=np.outer(np.cos(phi), np.sin(theta))
y=np.outer(np.sin(phi), np.sin(theta))
z=np.outer(np.ones(np.size(phi)), np.cos(theta))

# prepare function to plot
PHI=np.outer(phi,np.ones(np.size(theta)))
THETA=np.outer(np.ones(np.size(phi)),theta)
data = PHI/np.pi

# plot
surface=ax.plot_surface(x, y, z, cstride=1, rstride=1, facecolors=cm.jet(data),cmap=plt.get_cmap('jet'))

# add colorbar
m = cm.ScalarMappable(cmap=surface.cmap,norm=surface.norm)
m.set_array(data)
plt.colorbar(m)
plt.show()


# In[6]:


from ipyleaflet import *

m = Map(center=(52, 10), zoom=8, basemap=basemaps.Hydda.Full)
m


# In[10]:


import geopandas
geoplot.polyplot(world, figsize=(8, 4))


# In[ ]:




