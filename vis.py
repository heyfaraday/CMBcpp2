from numpy import genfromtxt, zeros

from math import pi
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

map = genfromtxt('bin/out.dat')

N = 512
projection = 'cyl'  # 'cyl', 'moll', 'ortho'
save_as_png = False
save_as_svg = False

inside_map = zeros((int(N + 1), int(N / 2 + 1)))
x = zeros((int(N + 1), int(N / 2 + 1)))
y = zeros((int(N + 1), int(N / 2 + 1)))

for i in range(0, int(N + 1)):
    for j in range(0, int(N / 2 + 1)):
        x[i][j] = 2.0 * i / N * pi
        y[i][j] = 2.0 * j / N * pi

for i in range(0, int(N + 1) * int(N / 2 + 1)):
    inside_map[int(map[i][1] - 1)][int(map[i][2] - 1)] = map[i][0]

rad = 180.0 / pi

fig = plt.figure(figsize=(8, 4))
fig.subplots_adjust(
    left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)

ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.axis('off')

cmb_map = Basemap(projection=projection, lon_0=0.0, lat_0=0.0, resolution='l')
cmb_map.contourf(
    x * rad, y * rad, inside_map, 512, cmap=plt.cm.jet, latlon=True)

if save_as_png:
    plt.savefig('out.png', dpi=300)
if save_as_svg:
    plt.savefig('out.svg')

plt.show()