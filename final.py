import osmnx as ox
import matplotlib.pyplot as plt
from shapely import Polygon
import geopy.distance
import sys

### Download map

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
for i in range(n):
	a, b = (float(x) for x in input().split())
	coords.append((a, b));

tour = [int(x) for x in input().split()]

axis = plt.axes()
ox.plot_graph(G, ax = axis, figsize = (40, 40), node_size=1, show=False)

final_ans = 0

#fig, ax = plt.subplots()

for i in range(len(tour)):
	j = (i+1)%len(tour)
	x_values = [coords[tour[i]][0], coords[tour[j]][0]]
	y_values = [coords[tour[i]][1], coords[tour[j]][1]]
	plt.plot(x_values, y_values, color="red", lw=0.8)
	final_ans += dist(coords[tour[i]], coords[tour[j]])

#def update(frame):
#    # for each frame, update the data stored on each artist.
#    x = t[:frame]
#    y = z[:frame]
#    # update the scatter plot:
#    data = np.stack([x, y]).T
#    scat.set_offsets(data)
#    # update the line plot:
#    line2.set_xdata(t[:frame])
#    line2.set_ydata(z2[:frame])
#    return (scat, line2)
#
#ani = FuncAnimation(fig, update, frames=40, interval=30)
#
#ani.save("test.gif",writer="imagemagick")
plt.savefig('final_output.pdf')
print(final_ans)
