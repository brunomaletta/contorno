import osmnx as ox
import matplotlib.pyplot as plt
from shapely import Polygon
import geopy.distance
import sys
import matplotlib.animation as animation
import matplotlib.patches as patches
import numpy as np

enclosing_coordinates = [
    (-43.9442396, -19.9389905), (-43.9374161, -19.9405336), (-43.9268053, -19.9394141),
    (-43.9258611, -19.9389300), (-43.9239407, -19.9318799), (-43.9200354, -19.9276336),
    (-43.9186621, -19.9227113), (-43.9189732, -19.9223583), (-43.9237797, -19.9211781),
    (-43.9288008, -19.9134111), (-43.9308929, -19.9129269), (-43.9338541, -19.9146014),
    (-43.9368796, -19.9139255), (-43.9387894, -19.9121905), (-43.9393902, -19.9119686),
    (-43.9402807, -19.9119888), (-43.9411283, -19.9122914), (-43.9439607, -19.9138247),
    (-43.9471257, -19.9172543), (-43.9477909, -19.9175267), (-43.9537132, -19.9171333),
    (-43.9566743, -19.9182731), (-43.9586914, -19.9260097), (-43.9570498, -19.9304679),
    (-43.9538205, -19.9339173), (-43.9460313, -19.9358135), (-43.9442396, -19.9389905),
]
pl = Polygon(enclosing_coordinates)
G = ox.graph_from_polygon(pl, network_type = 'walk')

def dist(coord1, coord2):
	return geopy.distance.geodesic(coord1, coord2).m

n = int(input())
coords = []
neigh = []
for i in range(n):
	a, b, c = (float(x) for x in input().split())
	coords.append((a, b));
	neigh.append(int(round(c)))

tour = [int(x) for x in input().split()]

fig, ax = plt.subplots()
ox.plot_graph(G, ax = ax, figsize = (40, 40), node_size=1, show=False)

final_ans = 0

print ("""<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.0"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        creator="Bruno Monteiro"
        xmlns="http://www.topografix.com/GPX/1/0"
        xsi:schemaLocation="http://www.topografix.com/GPX/1/0 http://www.topografix.com/GPX/1/0/gpx.xsd">
<trk><trkseg>""")

for i in range(len(tour)):
	j = (i+1)%len(tour)
	x_values = [coords[tour[i]][0], coords[tour[j]][0]]
	y_values = [coords[tour[i]][1], coords[tour[j]][1]]
	final_ans += dist(coords[tour[i]], coords[tour[j]])
	print("<trkpt lon=\"" + str(coords[tour[i]][0]) + "\" lat=\"" + str(coords[tour[i]][1]) + "\"></trkpt>")
print("""</trkseg></trk>
</gpx>
""")

def hex_to_RGB(hex_str):
    """ #FFFFFF -> [255,255,255]"""
    #Pass 16 to the integer function for change of base
    return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]
def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1))/255
    c2_rgb = np.array(hex_to_RGB(c2))/255
    mix_pcts = [x/(n-1) for x in range(n)]
    rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
    return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]


color_dark_blue = "#090979"
color_red = "#FF0000"
color_green = "#00FF00"
color_blue = "#0000FF"
colors = get_color_gradient(color_dark_blue, color_red, len(tour))
colors_neigh = get_color_gradient(color_red, color_blue, max(neigh)+1)

#for i in range(n):
#	c = colors_neigh[neigh[i]]
#	ax.plot(coords[i][0], coords[i][1], color=c, marker='o')

(x_text, y_text) = enclosing_coordinates[0]
x_text -= 0.015
text = ax.text(x_text, y_text, str(int(final_ans)/float(1000)) + " km")
x_values = [coords[tour[i]][0]]
y_values = [coords[tour[i]][1]]
ball = ax.scatter(x_values, y_values, color = colors[0])

final_ans = 0
def animate(i):
	global final_ans
	if i >= len(tour):
		return
	j = (i+1)%len(tour)
	x_values = [coords[tour[i]][0], coords[tour[j]][0]]
	y_values = [coords[tour[i]][1], coords[tour[j]][1]]
	final_ans += dist(coords[tour[i]], coords[tour[j]])
	ax.plot(x_values, y_values, color=colors[i], lw=1.5)

	#text.set_text(f'{i*100//len(tour)}%')
	#text.set_text(f'{int(final_ans)/float(1000)} km')
	text.set_text('{:.1f} km'.format(int(final_ans)/float(1000)))
	ball.set_offsets(np.c_[x_values[-1:], y_values[-1:]])
	ball.set_facecolor([colors[i]])
	ball.set_edgecolor([colors[i]])

anim = animation.FuncAnimation(fig, animate, frames=len(tour)+30*3, interval=10)

anim.save("tour.gif",writer="imagemagick")
plt.savefig('final_output.pdf')
