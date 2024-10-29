import osmnx as ox
import matplotlib
import sys
import geopy.distance

neigh = []
neigh.append('Lourdes')
neigh.append('Santa Efigênia')
neigh.append('Centro')
neigh.append('Carlos Prates (HC)')
neigh.append('Floresta (HC)')
neigh.append('Savassi')
neigh.append('Colégio Batista (HC)')
neigh.append('Santo Agostinho')
neigh.append('Funcionários')
neigh.append('Barro Preto')
neigh.append('Boa Viagem')
neigh.append('Bonfim (HC)')
neigh.append('Lagoinha (HC)')
for i in range(len(neigh)):
	neigh[i] += ', Belo Horizonte, Brazil'

pols = []
for i in neigh:
	g = ox.graph_from_place(i)
	x = g.nodes().data('x')
	y = g.nodes().data('y')
	coord_list = []

	coords = []
	for j in g.nodes():
		coords.append((x[j], y[j]))
	pols.append(coords)

#G = ox.graph_from_place(neigh, network_type='walk')
#ox.plot_graph(G, filepath='neigh.svg', save=True, show=False, close=True)

def dist(coord1, coord2):
	return geopy.distance.geodesic(coord1, coord2).m

nodes = True
for line in sys.stdin:
	line = line[:-1]
	if len(line) == 0:
		nodes = False
	if not nodes:
		print(line)
		continue

	nums = [float(i) for i in line.split()]
	neigh = -1
	min_dist = 100000000
	for i in range(len(pols)):
		for j in pols[i]:
			d = dist(j, (nums[0], nums[1]))
			if d < min_dist:
				min_dist = d
				neigh = i
	line += " " + str(neigh)
	print(line)
