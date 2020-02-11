import matplotlib.pyplot as plt
import numpy as np

d = 6771000
p = 10**5

v = 7672
p = 100

t, x, y, d, vx, vy = [], [], [], [], [] ,[]
with open("orbit.txt") as f:
    for line in f:
        j, ti, xi, yi, di, vxi, vyi = map(float, line.split())
        t.append(ti)
        x.append(xi)
        y.append(yi)
        d.append(di)
        vx.append(vxi)
        vy.append(vyi)

t = np.array(t)
x = np.array(x)
y = np.array(y)
d = np.array(d)
vx = np.array(vx)
vy = np.array(vy)

plt.figure(figsize=(10, 5))

# ORBIT
plt.subplot(131)
plt.scatter(x, y, s=1, color="r")
plt.axis('equal')
plt.title("Orbit")

# EARTH
earth = plt.Circle((0, 0), 6.371e+6, color='blue', fill=False)
plt.gca().add_artist(earth)

# DISTANCE
plt.subplot(235)
plt.plot(t, d)
plt.gca().set_xlim(0, max(t))

plt.title("Distance")

# SPEED
plt.subplot(232)
plt.plot(t, np.sqrt(vx**2 + vy**2))
plt.axis('equal')
plt.title("Speed")

# VELOCITY
plt.subplot(233)
plt.plot(t, vx)
plt.axis('equal')
plt.title("Vx")

plt.subplot(236)
plt.plot(t, vy)
plt.axis('equal')
plt.title("Vy")

# plt.axis('off')
plt.tight_layout()

plt.savefig("graph.png")

