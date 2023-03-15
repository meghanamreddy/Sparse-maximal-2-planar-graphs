#include "hds.h"
#include <sstream>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

// K_9 minus (P_3 + P_2 + P_2) => unique 2-plane drawing
const std::size_t n = 9;
const Edges edges = {
  {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},{0,8},
  {1,2},{1,3},{1,4},{1,5},{1,6},{1,7},{1,8},
  //{2,3},
  {2,4},{2,5},{2,6},{2,7},{2,8},
  {3,4},{3,5},{3,6},{3,7},{3,8},
  //{4,5},
  {4,6},{4,7},{4,8},
  {5,6},{5,7},{5,8},
  {6,7},//{6,8}
  //,{7,8}
};

int main()
{
  std::vector< Drawing<2> > solutions;
  Drawing<2> d(n);
  d.add_first_edge(edges[0][0], edges[0][1]);
  for (auto e = edges.begin() + 1;;) {
    std::size_t u = (*e)[0];
    std::size_t v = (*e)[1];
    HdsPath p = d.first_path(u, v);
    if (p.empty()) {
    BACKUP:
      // no way to add uv -> do previous edges differently
      do {
	if (--e == edges.begin()) goto END;
	u = (*e)[0];
	assert(u == d.edges.back().u);
	v = (*e)[1];
	assert(v == d.edges.back().v);
	p = d.edges.back().built;
	d.remove_edge();
      } while (!d.next_path(p, v));
    }
    d.add_edge(p, v);
    if (++e == edges.end()) {
      bool is_new = true;
      if (is_new) {
	solutions.push_back(d);
	std::cout << "Drawing #" << solutions.size() << ":\n" << d;
	std::ofstream of;
	std::ostringstream filename;
	filename << "gadget-X-" << solutions.size() << ".graphml";
	of.open(filename.str());
	d.graphml_output(of);
	of.close();
      }
      goto BACKUP;
    }
  }
 END:
  std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
  return 0;
}
