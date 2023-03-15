# Sparse-maximal-2-planar-graphs
Code repository for the paper titled "The number of edges in maximal 2-planar graphs" by Michael Hoffmann and Meghana M Reddy

You need a C++ compiler, such g++ (https://gcc.gnu.org/), as well as the Boost Graph Library
(BGL, cf. https://www.boost.org/doc/libs/1_80_0/libs/graph/doc/) to build these programs.

There are two main files: iterative_nested_decagons.cpp and gadget-X.cpp.
The build commands are:

g++ -std=c++11 -O2 gadget-X.cpp -ogadget-X 
g++ -std=c++11 -O2 iterative_nested_decagons.cpp -oiterative_nested_decagons 

If the BGL headers are not in your system search path, but in directory DIR add the option -IDIR.
For instance, if the headers are in /opt/local/include, then use the commands:

g++ -std=c++11 -I/opt/local/include -O2 gadget-X.cpp -ogadget-X 
g++ -std=c++11 -I/opt/local/include -O2 iterative_nested_decagons.cpp -oiterative_nested_decagons 

The computation for the gadget X takes about 2-3 min., the one for G_k about 1.5 days.

The .graphml files (see http://graphml.graphdrawing.org/) can, for
instance, be viewed with the freely available graph editor/viewer yEd
(https://www.yworks.com/products/yed). To get a plane layout in yEd,
load the graphml file and then select from the menu Layout ->
Orthogonal -> Classic, and press ok (or play with the parameters).
