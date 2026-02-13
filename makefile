all:
	python3 get_graph.py < polyline.txt > graph.txt
	python3 get_neighborhoods.py < graph.txt > graph_neigh.txt
	g++ -std=c++17 -O2 -o tour tour.cpp
	./tour < graph_neigh.txt > tour.txt
	rm tour
	python3 final.py < tour.txt > tour.gpx
	ffmpeg -i tour.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" video.mp4
