import matplotlib.pyplot as plt
import numpy as np

d = 6771000
p = 10**5

v = 7672
p = 100

t, x, y, d, vx, vy, mx, my = [], [], [], [], [], [], [], []
with open("orbit.txt") as f:
    for line in f:
        j, ti, xi, yi, di, vxi, vyi, mxi, myi = map(float, line.split())
        t.append(ti)
        x.append(xi)
        y.append(yi)
        d.append(di)
        vx.append(vxi)
        vy.append(vyi)
        mx.append(mxi)
        my.append(myi)

t = np.array(t)
x = np.array(x)
y = np.array(y)
d = np.array(d)
vx = np.array(vx)
vy = np.array(vy)
mx = np.array(mx)
my = np.array(my)

# ORBIT
plt.subplot(121)
plt.scatter(x, y, s=1, color="r")
plt.axis('equal')
plt.gca().set_xlim(0, max(max(x), max(mx)))
plt.title("Orbit")

# EARTH
earth = plt.Circle((0, 0), 10*6.371e+6, color='blue', fill=False)
plt.gca().add_artist(earth)

# MOON
plt.scatter(mx, my, s=1, color="g")

# DISTANCE
plt.subplot(222)
plt.plot(t, d)
plt.gca().set_xlim(0, max(t))

plt.title("Distance")

# SPEED
plt.subplot(224)
plt.plot(t, np.sqrt(vx**2 + vy**2))
#  plt.axis('equal')
plt.title("Speed")

# plt.axis('off')
plt.tight_layout()

#plt.savefig("graph.png")
plt.show()

"""
start = 120000
for i in range(start, start + 20000, 500):
    print(t[i], np.sqrt(vx[i]**2 + vy[i]**2))

print(np.where(t == 228802.0))

start = 132500
for i in range(start - 500, start + 500, 50):
    print(t[i], np.sqrt(vx[i]**2 + vy[i]**2))
"""

# print((228752 - 96302.0)/(60*60*24))
