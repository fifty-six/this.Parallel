import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import pandas as pd

EARTH = 10
SPEED = 400

data = pd.read_csv("orbit.txt", delimiter=" ", header=None, index_col=0)
# print(data)

t, x, y, d, vx, vy, mx, my = data.T.values

def update_plot(i, _x, _y, scatter):
    scatter.set_offsets(np.c_[_x[:i * SPEED], _y[:i * SPEED]])
    return scatter

fig, ax = plt.subplots(figsize = (10, 6))

# Remove axis for what seems to be literally nothing
plt.axis("off")

# ROCKET ORBIT
orbit_plot = fig.add_subplot(121)
orbit_plot.grid()
scatter_rocket = orbit_plot.scatter(x, y, s=1, color="r")
orbit_plot.axis('equal')
orbit_plot.title.set_text("Orbit")

# EARTH
earth = plt.Circle((0, 0), EARTH*6.371e+6, color='blue', fill=False)
fig.gca().add_artist(earth)

# MOON
scatter_moon = orbit_plot.scatter(mx, my, s=1, color="g")

anim_rocket = animation.FuncAnimation(fig, update_plot, frames=range(len(t / SPEED)), fargs=(x, y, scatter_rocket), interval=50)
anim_moon = animation.FuncAnimation(fig, update_plot, frames=range(len(t / SPEED)), fargs=(mx, my, scatter_moon), interval=50)

# DISTANCE
dist_plot = fig.add_subplot(222)
dist_plot.grid()
dist_plot.plot(t, d)
fig.gca().set_xlim(0, max(t))
dist_plot.title.set_text("Distance")

# SPEED
speed_plot = fig.add_subplot(224)
speed_plot.grid()
speed_plot.plot(t, np.sqrt(vx**2 + vy**2))
speed_plot.title.set_text( "Speed")

fig.tight_layout(pad=2)

# plt.savefig("graph.png")
plt.show()
