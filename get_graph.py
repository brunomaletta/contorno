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

# Limit distance, in meters
DIST_LIM = 5.0

# Epsilon
EPS = 1e-9

def dist(coord1, coord2):
	return geopy.distance.geodesic(coord1, coord2).m

coords = []

started = False
read_nums = []
points = []
for line in sys.stdin:
	if line.find('coordinates') > 0:
		started = True
		continue
	if not started:
		continue
	
	line = line.replace(',', '')
	#print(line, end='')
	try:
		num = float(line)
		read_nums.append(num)
	except:
		pass

	if len(read_nums) == 2:
		points.append((read_nums[0], read_nums[1]))
		read_nums.clear()

f = [i for i in range(len(points))]
def find(i):
	if f[i] == i:
		return i
	f[i] = find(f[i])
	return f[i]
def union(i, j):
	i = find(i)
	j = find(j)
	if i != j:
		f[i] = j

for i in range(len(points)):
	for j in range(len(points)):
		if j <= i:
			continue;
		if dist(points[i], points[j]) < DIST_LIM:
			union(i, j)

edges = []
for i in range(len(points)-1):
	id1 = find(i)
	id2 = find(i+1)
	edges.append((id1, id2))

for (lat, long) in points:
	print(lat, long)
print()
for (a, b) in edges:
	print(a, b, dist(points[a], points[b]))

axis = plt.axes()
ox.plot_graph(G, ax = axis, figsize = (40, 40), node_size=1, show=False)

for (a, b) in edges:
	x_values = [points[a][0], points[b][0]]
	y_values = [points[a][1], points[b][1]]
	plt.plot(x_values, y_values, color="red", lw=0.8)

plt.savefig('output.pdf')
