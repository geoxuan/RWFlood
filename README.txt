RWFlood:
A new and faster internal memory method to compute the drainage network, that is, the flow direction and 
accumulation on terrains represented by raster elevation matrix.

Compile:
g++ -O3 hydrog.cpp -o hydrog

Run:
./hydrog < terrain_name.hgt  nrows

Output:
As a result it is given a file named flow.hgt containing the accumulated flow with 2 bytes.
