from numpy import genfromtxt, zeros, size

from math import pi
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

map_from_file = genfromtxt('../data/CMB_smooth.dat')
points_from_file = genfromtxt('../data/CMB_extrema_smooth.dat').T
# sing_from_file = genfromtxt('../bin/hfi_sing.dat').T

N = 4096
projection = 'cyl'
save_as_png = False
save_as_svg = False

inside_map = zeros((int(N + 1), int(N / 2 + 1)))
x = zeros((int(N + 1), int(N / 2 + 1)))
y = zeros((int(N + 1), int(N / 2 + 1)))

for i in range(0, int(N + 1)):
    for j in range(0, int(N / 2 + 1)):
        x[i][j] = 2.0 * i / N * pi - pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

for i in range(0, int(N + 1) * int(N / 2 + 1)):
    inside_map[int(map_from_file[i][0])][int(map_from_file[i][1])] = map_from_file[i][2]

rad = 180.0 / pi

fig = plt.figure(figsize=(8, 4))
fig.subplots_adjust(
    left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)

ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.axis('off')

cmbmap = Basemap(projection=projection, lon_0=0, resolution='l')
cmbmap.contourf(x * rad, y * rad, inside_map, 300, cmap=plt.cm.jet, latlon=True)

points_marker = ''

for i in range(0, 2000):
    if points_from_file[3][i] == 0:
        points_marker = '.'
    elif points_from_file[3][i] == 1:
        points_marker = 'x'
    elif points_from_file[3][i] == 2:
        points_marker = '+'
    cmbmap.scatter((points_from_file[0][i] - pi) * rad, (points_from_file[1][i] - pi / 2.0) * rad,
                   marker=points_marker, color='red')

cmbmap.contour(x * rad, y * rad, inside_map)
#
# sing_marker = ''
#
# for i in range(0, size(sing_from_file[0])):
#     if sing_from_file[2][i] == 1:
#         sing_marker = '.'
#     elif sing_from_file[2][i] == 2:
#         sing_marker = 'x'
#     elif sing_from_file[2][i] == 3:
#         sing_marker = '+'
#     cmbmap.scatter((sing_from_file[0][i] - pi) * rad, (sing_from_file[1][i] - pi / 2.0) * rad,
#                    marker=sing_marker, color='black')

# cmb_map = Basemap(projection=projection, lon_0=0.0, lat_0=0.0, resolution='l')
# cmb_map.contourf(x * rad, y * rad, inside_map, 1024, cmap=plt.cm.jet, latlon=True)

if save_as_png:
    plt.savefig('out.png', dpi=300)
if save_as_svg:
    plt.savefig('out.svg')

plt.show()
