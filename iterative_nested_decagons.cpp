// This program computes all 2-plane drawings of a family G_k of
// graphs, which can roughly be described as a sequence of ten-cycles
// any two consecutive of which are connected by a braided perfect
// matching. A more precise description follows.
//
// For an integer k >= 2 the graph G_k has 10k+140 vertices and
// 20k+630 edges. The vertices are distributed over a sequence of k
// pairwise vertex-disjoint ten-cycles D_1,...,D_k, as well as 20
// copies of a special gadget graph X that has 9 vertices and 32
// edges. These copies of X are (perfectly) matched with the edges of
// D_1 and D_k such that one edge of the copy of X is glued with the
// corresponding edge of D_1 or D_k. In this way, every copy of X adds
// 7 vertices and 31 edges, which leads to the overall count of
// 10k+20*7=10k+140 vertices and 10k+(k-1)10+20*31=20k+610 edges. The
// remaining 20 edges are added to cyclically connect the vertices of
// D_1 and D_k, respectively, by connecting all pairs of vertices that
// have distance exactly two along each of these cycles.
//
// The gadget graph X is analyzed separately (in another program
// gadget-X.cpp), to show that for a specific edge e of X (which is
// exactly the edge that we glue to the corresponding edge of D_1 or
// D_k) in every 2-plane drawing of X both incident faces are
// triangles that are bounded the uncrossed edge e and two doubly
// crossed edges. As a consequence, we can regard the edge e as
// uncrossable in every 2-plane drawing of G_k: As G_k is 3-connected
// and e can be crossed at most twice, all vertices outside of this
// copy of X have to be in the same face of the corresponding X
// subdrawing. Therefore, we do not consider the copies of X here and
// instead mark the edges of D_1 and D_k as uncrossable.
//
// Actually, we do not explicitly consider the length two chords that
// connect the vertices of D_1 and D_k, respectively, either. As the
// cycle edges are uncrossable, all remaining vertices have to be on
// one side of the cycle. We pick this side w.l.o.g. due to symmetry,
// and just block the other side with an uncrossable star, centered at
// an additional dummy vertex.
//
// So instead of the graphs G_k as defined above we consider the
// graphs G_k^-, which are derived from G_k by
//   (1) removing all X gadgets,
//   (2) removing the length two edges that connect the vertices of
//       the first and the last cycle D_1 and D_k, and
//   (3) marking the edges of D_1 and D_k as uncrossable.
//
// In order to handle the infinite family G_k^- in a finite
// computation, we iteratively compute all 2-plane drawings of the
// graphs (G_i^-)', for i=2,3,..., where (G_i^-)' is identical to
// G_i^-, except that the edges of D_i in (G_i^-)' are {\em
// unconstrained}, that is, we allow them to be crossed. (The edges of
// the first cycle D_1 are marked uncrossable in (G_i^-)' as well.)
// Considering (G_i^-)' makes sense because every 2-plane drawing of
// G_k^- contains a 2-plane subdrawing of (G_i^-)', for all 2 <= i <=
// k-1, but not necessarily of all these G_i^-.
//
// In other words, in order to enumerate all 2-plane drawings of
// G_k^-, it suffices to consider all 2-plane drawings of (G_{k-1}^-)'
// and try all possible ways to extend them to 2-plane drawings of
// (G_k^-)'. Drawings that admit such an extension we call {\em
// extensible}. Whenever we discover a (new) extensible drawing, we
// also check whether it is {\em completeable} to a 2-plane drawing of
// G_k^-, that is, where the added ten-cycle is prescribed to be
// crossing-free. (Clearly, every drawing that is completeable is also
// extensible.)
//
// The restriction to drawings of G_k^- does not suffice to make the
// computation finite. After all, we know that there exists at least
// one 2-plane drawing of G_k and, thus, of G_k^-, for every k.  To
// further restrict the search space, let us consider a subdrawing of
// (G_k^-)' in a drawing of (G_l^-)', for some l > k. The cycle D_k of
// (G_k^-)', let us call it the {\em last} cycle of G_k, forms a cut
// set in G_l that separates the vertices of D_1,...,D_{k-1} from the
// vertices of D_{k+1},...,D_l. Therefore, all remaining vertices of
// (G_k^-)', which are not on its last cycle, and all the edges and
// faces of the drawing of (G_k^-)' are only relevant for the
// extensibility of this drawing if they either contain a vertex of
// the last cycle, or they are "close enough" to such a vertex. This
// notion of "close enough" is made more precise by a classifications
// of the faces of the drawing into active, passive, transit, and
// irrelevant faces.
//
// The reasoning behind this classification and the proofs that show
// that it really captures what we need can be found in the
// paper. Here we just provide the definitions along with some
// intuition. For the definition we consider flow networks in the dual
// multigraph of the drawing (whose vertices are the faces, and they
// are connected via each shared boundary edge; the capacity of such a
// network edge is 2-#(crossings on that edge in the drawing)),
// augmented by the vertices of the last cycle (connected to all their
// incident faces, as well as to a single artificial source, with unit
// capacity edges).
//
// Then a face f of the drawing is a {\em potential final face} (or
// PFF, for short) if the network with sink f has a flow of value at
// least 10. The existence of a PFF is necessary for a drawing to be
// completeable. (Intuition: As the edges of D_l are uncrossable, all
// vertices of D_l lie in the same face of the drawing of
// (G_k^-)'. Moreover, there is a collection of pairwise
// vertex-disjoint paths that connect the vertices of D_l to those of
// D_k.) Our classification of faces, as detailed next, is based on
// this notion of a pff. A face in a drawing of (G_k^-)' is ...
//
//  * {\em active} if it has a flow of value >= 3 to some pff;
//
//  * {\em passive} if it is not active, it has a vertex v of the last
//    cycle D_k on its boundary, and in the dual network it has a path
//    P of length <= 2 to some active face such that the first edge of
//    P is not incident to v;
//
//  * {\em transit} if it is neither active nor passive, but in the
//    dual network it has at least one edge to an active face and at
//    least one other edge that leads to an active or passive face;
//
//  * {\em irrelevant} if it is neither active, nor passive, nor transit.
//
// In the paper we show that the notion "irrelevant" is justified in
// the sense that in every 2-plane drawing of (G_l^-)' no vertex or
// edge of (G_l^-)' - (G_k^-)' interacts with any irrelevant
// face. Therefore, we can remove all irrelevant faces from each
// drawing we discover and proceed to compute with the remaining, {\em
// relevant} part of the drawing only. The hope is that this culling
// suffices to cut down the search space to a finite number of
// drawings to consider. (Intuition: Given that at most two crossings
// per edge are allowed, the vertices and edges of the last cycle
// cannot be distributed too wildly, but only in some small,
// constant-size neighborhood of a PFF.) 
//
// Output format: The output is provided in a plain text format, which
// lists the rotation system of each discovered drawing. That is, for
// each vertex, the circular sequence of neighbors is listed. The
// rotation system by itself does not suffice to uniquely determine
// the drawing. Thus, with each neighbor in the circular sequence, if
// the edge to this neighbor has crossings, we also list the crossed
// edges in parentheses, in the order and orientation (left to right)
// they are crossed along the way. As an example, here is such a
// description of a (the) 1-plane drawing of K_6:
//
// { #vertices = 6, #edges = 15, #halfedges = 42, #crossings = 3
// 	0 : 1(x2-5) 2 3(x4-2) 4 5 
// 	1 : 0(x5-2) 5 4(x3-5) 3 2 
// 	2 : 0 5(x1-0) 1 3 4(x0-3) 
// 	3 : 0(x2-4) 2 1 5(x4-1) 4 
// 	4 : 0 2(x3-0) 3 1(x5-3) 5 
// 	5 : 0 4 3(x1-4) 1 2(x0-1) 
// }
//
// In addition, every drawing is also saved into a graphml file (see
// http://graphml.graphdrawing.org/). These files can, for instance,
// be viewed with the freely available graph editor/viewer yEd
// (https://www.yworks.com/products/yed). To get a plane layout in
// yEd, load the graphml file and then select from the menu Layout ->
// Orthogonal -> Classic, and press ok (or play with the parameters).
// To preserve the drawing - and not just the graph - in the graphml
// representation, we also add the dual, that is, a new vertex for
// each face, connected to all original vertices on the face
// boundary. The avoid too much clutter, these edges are shown with
// zero width. (Color and edge/vertex size parameters use a yEd
// specific graphml extension. So they only work with yEd probably.)

#include "hds.h"
#include <sstream>
#include <fstream>
#include <deque>
#include <vector>
#include <iostream>

// Return true iff the two drawings are isomorphic, in the following
// sense. We assume that in both drawings, the last ten-cycle (that
// is, D_k) is formed by the vertices [0,...,9]. We only allow
// isomorphisms that fix this ten-cycle (as a whole) and also keep the
// parity (odd/even) of its indices.
// PRE: Both d1 and d2 are plane (no crossings) and connected.
bool are_isomorphic(const Drawing<2>& d1, const Drawing<2>& d2) {
  if (d1.crossings.size() > 0 || d2.crossings.size() > 0)
    throw std::runtime_error("no crossings here, please");
  if (d1.vertices.size() != d2.vertices.size()) return false;
  if (d1.edges.size() != d2.edges.size()) return false;
  if (d1.halfedges.size() != d2.halfedges.size()) return false;

  const std::size_t n = d1.vertices.size();
  for (int mir = 0; mir < 2; ++mir) // mirror?
    for (int sh = 0; sh < 10; sh += 2) { // shift for cycle 0,...,9
      // map phi: V[d1] -> V[d2]
      std::vector<std::size_t> phi(n, n); // n == unitialized
      // inverse map
      std::vector<std::size_t> psi(n, n); // n == unitialized
      // fix vertices [0,9]
      if (mir == 0)
	for (std::size_t i = 0; i < 10; ++i) {
	  phi[i] = (i+sh) % 10;
	  psi[(i+sh) % 10] = i;
	}
      else
	for (std::size_t i = 0; i < 10; ++i) {
	  phi[i] = (10-i+sh) % 10;
	  psi[(10-i+sh) % 10] = i;
	}
      // dfs starting from the halfedge 1,0
      auto start = d1.vertices[0].halfedge;
      while (start->twin->vertex->label != 1) start = start->next->twin;
      std::vector<const HdsHalfedge*> stack(1, start);
      std::vector<bool> done(n, false); // vertex visited by dfs?
      done[0] = true;
      struct Contradiction {}; // exception to end this try
      try {
	while (!stack.empty()) {
	  auto h1 = stack.back();
	  stack.pop_back();
	  std::size_t u = h1->vertex->label;
	  std::size_t v = h1->twin->vertex->label;
	  if (phi[u] == n || phi[v] == n)
	    throw std::runtime_error("internal error in isomorphism test");
	  auto h2 = d2.vertices[phi[u]].halfedge;
	  // try to locate phi[v] in the rotation at phi[u]
	  for (;;) {
	    if (h2->twin->vertex->label == phi[v]) break;
	    h2 = h2->next->twin;
	    if (h2 == d2.vertices[phi[u]].halfedge) throw Contradiction();
	  }
	  // check if the rotations agree
	  auto end1 = h1;
	  auto end2 = h2;
	  do {
	    if (phi[h1->twin->vertex->label] == n) {
	      if (psi[h2->twin->vertex->label] != n) throw Contradiction();
	      phi[h1->twin->vertex->label] = h2->twin->vertex->label;
	      psi[h2->twin->vertex->label] = h1->twin->vertex->label;
	    } else if (phi[h1->twin->vertex->label] != h2->twin->vertex->label)
	      throw Contradiction();
	    if (!done[h1->twin->vertex->label]) {
	      stack.push_back(h1->twin);
	      done[h1->twin->vertex->label] = true;
	    }
	    h1 = h1->next->twin;
	    if (mir == 0) h2 = h2->next->twin; else h2 = h2->twin->prev;
	  } while (h1 != end1 && h2 != end2);
	  // does the degree match?
	  if (h1 != end1 || h2 != end2) throw Contradiction();
	}
	// this looks correct, let's verify the map
	for (std::size_t i = 0; i < n; ++i)
	  if (phi[i] == n || psi[phi[i]] != i)
	    std::runtime_error("disconnected drawing in isomorphism test");
	return true;
      } catch (Contradiction) {}
    }
  return false;
}

//  maxflow using BGL
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

struct Network {
  typedef boost::adjacency_list_traits<boost::vecS,
				       boost::vecS,
				       boost::directedS>
  boost_graph_traits;
  
  typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::directedS,
    boost::no_property,
    boost::property<boost::edge_capacity_t, int,
      boost::property<boost::edge_residual_capacity_t, int,
        boost::property<boost::edge_reverse_t, boost_graph_traits::edge_descriptor>
      >
    >
  >
  Graph;
  
  Network(std::size_t n) : G(n) {}
  
  void add_edge(int from, int to, int capacity) {
    auto c_map = boost::get(boost::edge_capacity, G);
    auto r_map = boost::get(boost::edge_reverse, G);
    const auto e = boost::add_edge(from, to, G).first;
    const auto rev_e = boost::add_edge(to, from, G).first;
    c_map[e] = capacity;
    c_map[rev_e] = 0; // reverse edge has no capacity!
    r_map[e] = rev_e;
    r_map[rev_e] = e;
  }

  int flow(int from, int to) { return boost::push_relabel_max_flow(G, from, to); }

  Graph G;
};

// Store the marked drawing d, where the faces are colored with their
// activity level, into a graphml file. The parameters have the
// following semantics:
//
//  * concatenation of fileprefix and index is used as a filename;
//  * fc denotes the number of faces in d;
//  * vact maps vertex index to color (only two options: vact[i] ==
//    vact.size() means gray = irrelevant, otw it is blue for vertices
//    and orange for crossings);
//  * active maps face index to activity level (3 means PFF, shown
//    purple, 2 means active, shown cyan, 1 means passive, shown
//    white, 0 means transit, shown gray, and anything else means
//    irrelevant, show black);
//  * fhedge maps face index to some halfedge on the boundary of the face.
template < typename V, typename A, typename F >
void store_marked(const Drawing<2>& d,
		  const std::string& fileprefix,
		  std::size_t index,
		  std::size_t fc,
		  const V& vact,
		  const A& active,
		  const F& fhedge)
{
  std::ofstream of;
  std::ostringstream filename;
  filename << fileprefix << "-" << index << ".graphml";
  of.open(filename.str());
  d.graphml_output_header(of, vact);
  // add dual so that we (1) fix the embedding and (2) can see the
  // status of faces
  std::size_t fec = 0;
  for (std::size_t i = 0; i < fc; ++i) {
    std::string color = "000000";
    if (active[i] == 0) color = "C0C0C0";
    else if (active[i] == 1) color = "FFFFFF";
    else if (active[i] == 2) color = "CCFFFF";
    else if (active[i] == 3) color = "FF00FF";
    // add vertex
    of << "    <node id=\"f" << i << "\">\n"
       << "      <data key=\"d5\">\n"
       << "        <y:ShapeNode>\n"
       << "          <y:NodeLabel>" << i << "</y:NodeLabel>\n"
       << "          <y:Geometry height=\"20.0\" width=\"20.0\"/>\n"
       << "          <y:Fill color=\"#" << color << "\"/>\n"
       << "          <y:Shape type=\"diamond\"/>\n"
       << "        </y:ShapeNode>\n"
       << "      </data>\n"
       << "    </node>\n";
    
    // add edges
    auto e = fhedge[i];
    do {
      std::size_t u = e->vertex->label;
      of << "    <edge id=\"g" << fec++ << "\" source=\"";
      if (u < d.vertices.size())
	of << "v" << u;
      else
	of << "c" << u - d.vertices.size();
      of << "\" target=\"f" << i << "\">\n"
	 << "      <data key=\"d9\">\n"
	 << "        <y:PolyLineEdge>\n"
	 << "          <y:LineStyle color=\"#FFFFFF\" "
	 <<                "type=\"line\" width=\"1\"/>\n"
	 << "          <y:BendStyle smoothed=\"true\"/>\n"
	 << "        </y:PolyLineEdge>\n"
	 << "      </data>\n"
	 << "    </edge>\n";
      e = e->next;
    } while (e != fhedge[i]);
  }
  d.graphml_output_footer(of);
  of.close();
}

// Store the drawing d into a graphml file. The concatenation of
// fileprefix and index is used as a filename. No specific colors for
// vertices and faces are used.
void store_reduced(const Drawing<2>& d,
		   const std::string& fileprefix,
		   std::size_t index)
{
  // create map halfedge -> face
  std::vector<int> face(d.halfedges.size(), -1);
  // ... and the reverse: face -> some halfedge on boundary
  std::vector<const HdsHalfedge*> fhedge;
  std::size_t fc = 0;
  for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
    if (face[i->label] != -1) continue;
    const HdsHalfedge* j = &*i;
    fhedge.push_back(j);
    do {
      face[j->label] = fc;
      j = j->next;
    } while (j->label != i->label);
    ++fc;
  }

  // all vertices in default color
  std::vector<std::size_t> vact(d.vertices.size() + d.crossings.size(), 0);
  // all faces plain white
  std::vector<std::size_t> active(fc, 1);

  store_marked(d, fileprefix, index, fc, vact, active, fhedge);
}

// This function is used when extracting the relevant parts of a
// drawing to create the reduced drawing that is then used for further
// processing. Vertices are classified as relevant or irrelevant by
// the parameter vact: vact[v] == vact.size() for a vertex v <=> v is
// irrelevant. Irrelevant vertices do not appear in the reduced
// drawing, so for a halfedge e that points to a relevant vertex we
// want to know what is the other endpoint (i.e., source vertex) of
// the corresponding halfedge in the reduced drawing (which should be
// a relevant vertex, and it need not be the source of e in the
// original, non-reduced drawing). If only one side of e has a
// relevant face, then we look for the source along this relevant
// face. To find the face incident to a halfedge and to determine if
// it is relevant we use the maps face and active, respectively.
//
// The function returns a pair, the 1st component is a halfedge that
// points to the desired source vertex, the second component specifies
// the number of crossings that the edge should have in the reduced
// drawing (which is set to two <-> uncrossable, unless the faces on
// both sides are relevant).
template < typename F, typename A, typename V > inline
std::pair<const HdsHalfedge*,int>
prev_active(const HdsHalfedge* e, const F& face, const A& active, const V& vact) {
  if (active[face[e->label]] >= 0) {
    int c = (active[face[e->twin->label]] < 0 ? 2 : e->edge->ncr);
    do e = e->prev; while (vact[e->vertex->label] == vact.size());
    return std::make_pair(e, c);
  }
  e = e->twin;
  if (active[face[e->label]] < 0) throw std::runtime_error("no active face");
  while (vact[e->vertex->label] == vact.size()) e = e->next;
  return std::make_pair(e, (active[face[e->twin->label]] < 0 ? 2 : e->edge->ncr));
}

int main()
{
  std::vector< Drawing<2> > solutions; // store all solutions found

  // We build D_1 as a universal starting configuration: a decagon,
  // where all edges are marked as uncrossable. As the cycle cannot be
  // crossed, all remaining vertices are on one side of it. By
  // symmetry, we may pick this side without loss of generality. So we
  // block the other side with an uncrossable plane star, centered at
  // a dummy vertex (which preserves symmetry of the drawing).
  {
    Drawing<2> d(11);
    // first the cycle; store the cycle edges to build the star later
    std::vector<HdsHalfedge*> cycle(10, 0); // [i] -> edge to i
    cycle[1] = d.add_first_edge(0, 1, 2);
    for (std::size_t i = 2; i <= 9; ++i)
      cycle[i] = d.add_edge(HdsPath({cycle[i-1], 0}), i, 2);
    cycle[0] = d.add_edge(HdsPath({cycle[9], cycle[1]->twin}), 0, 2);
    // then the star
    auto e = d.add_edge(HdsPath({cycle[0], 0}), 10, 2);
    for (std::size_t i = 9; i > 0; --i)
      e = d.add_edge(HdsPath({cycle[i], e}), 10, 2);
    solutions.push_back(d);
    store_reduced(d, "reduced", 0);
  }
 
  // #discarded solutions (isomorphic to a known solution)
  std::size_t discarded = 0;
  // we process the solutions in order, this is the current one we process
  std::size_t solcount = 0;
  // #full solutions (i.e., 2-plane drawings of G_k); we claim that
  // this never exceeds one...
  std::size_t fullsolcount = 0;
  
  // main loop
  while (solcount < solutions.size()) {
    std::cout << "=== Considering solution #" << solcount << " (>= " 
	      << solutions.size() - solcount - 1 << " more to go)" << std::endl;
    // remember current #solutions to see if we found any new ones this iteration
    std::size_t solsofar = solutions.size();
    std::size_t fullsolsofar = fullsolcount;
    
    // We run the computation twice: First we try to add another
    // unconstrained cycle and record all new descendants. Second,
    // assuming that some completion was found in the first round, we
    // try to add another uncrossable cycle.

    // record if we found a valid solution (even if it might not be
    // new because it is isomorphic to one that we found earlier)
    bool is_extensible = false;
    
    for (int constrained = 0; constrained < 2; ++constrained) {
      // see what drawings we get when adding a new decagon to solutions[solcount]
      Drawing<2> d = solutions[solcount];
      std::size_t nm10 = d.vertices.size();
      d.add_vertices(10);
      typedef std::vector<std::size_t> Edge;
      std::vector<Edge> edges = {
	{0,nm10+8},
	{nm10+8,nm10+9,2},{nm10+9,nm10,2},
	{nm10,nm10+1,2},{nm10+1,nm10+2,2},
	{nm10+2,nm10+3,2},{nm10+3,nm10+4,2},
	{nm10+4,nm10+5,2},{nm10+5,nm10+6,2},
	{nm10+6,nm10+7,2},{nm10+7,nm10+8,2},
	{2,nm10},{4,nm10+2},{6,nm10+4},{8,nm10+6},
	{1,nm10+3},{3,nm10+5},{5,nm10+7},{7,nm10+9},{9,nm10+1}
      };

      // exception to end the current try
      struct done_w_current_drawing {};
      
      // if we can't add a cycle without constraints, we won't be able
      // to add one with
      if (constrained == 0 || is_extensible) 
	try {
	  // edge adding loop
	  for (auto e = edges.begin();;) {
	    std::size_t u = (*e)[0];
	    std::size_t v = (*e)[1];
	    int pcr = 0; // prescribed #crossings (== 2 ==> uncrossable)
	    if (constrained == 1 && e->size() == 3) pcr = (*e)[2];
	    HdsPath p = d.first_path(u, v, pcr);
	    if (p.empty()) {
	    BACKUP:
	      // no way to add uv -> do previous edges differently
	      do {
		if (e == edges.begin()) throw done_w_current_drawing();
		--e;
		u = (*e)[0];
		assert(u == d.edges.back().u);
		v = (*e)[1];
		assert(v == d.edges.back().v);
		pcr = 0;
		if (constrained == 1 && e->size() == 3) pcr = (*e)[2];
		p = d.edges.back().built;
		//std::cerr << "remove edge " << u << "-" << v << std::endl;
		d.remove_edge();
	      } while (!d.next_path(p, v, pcr));
	    }
	    d.add_edge(p, v, pcr);
	    if (++e == edges.end()) {
	      // for the constrained cycle this is good enough
	      if (constrained == 1) {
		std::cout << "*** found solution #" << fullsolcount << std::endl;
		std::ofstream of;
		std::ostringstream filename;
		filename << "solution-" << fullsolcount++ << ".graphml";
		of.open(filename.str());
		d.graphml_output(of);
		of.close();
		goto BACKUP;
	      }
	      
	      // to be valid, there needs to be a face that all vertices
	      // of the unconstrained cycle C_u = cub,...,cue-1 can access
	      // (via an edge that has no crossing restrictions)
	      std::size_t cue = d.vertices.size();
	      std::size_t cub = cue - 10;
	      
	      // first create map halfedge -> face
	      std::vector<int> face(d.halfedges.size(), -1);
	      std::size_t fc = 0;
	      for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
		if (face[i->label] != -1) continue;
		const HdsHalfedge* j = &*i;
		do {
		  face[j->label] = fc;
		  j = j->next;
		} while (j->label != i->label);
		++fc;
	      }
	      
	      // create flow network for the dual graph
	      //
	      // vertices 0,...,fc-1 correspond to the faces of d
	      // vertices fc,...,fc+9 correspond to the vertices of C_u
	      // vertex fc+10 is the artificial single source
	      int source = fc + 10;
	      Network g(source + 1);
	      // also build a separate adjacency list representation of
	      // the dual graph (the flow network has artifacts, such as
	      // reverse edges and the artiticial source, etc.)
	      std::vector<std::vector<std::size_t> > dualg(fc);
	      // edges from source to vertices
	      for (std::size_t i = fc; i < fc + 10; ++i) g.add_edge(source, i, 1);
	      // edges between faces and from vertices to incident faces
	      for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
		if (i->vertex->label >= cub && i->vertex->label < cue)
		  g.add_edge(fc + i->vertex->label - cub, face[i->label], 1);
		int cap = 2 - i->edge->ncr; // edge capacity
		if (cap > 0) {
		  g.add_edge(face[i->label], face[i->twin->label], cap);
		  dualg[face[i->label]].push_back(face[i->twin->label]);
		}
	      }
	      
	      // The final cycle will be uncrossable, so it will go into a
	      // single face of the current drawing. Is there such a
	      // possible final face? It has to have a path from each
	      // vertex of C_u.
	      std::vector<int> pff; // possible final faces
	      for (std::size_t i = 0; i < fc; ++i)
		if (g.flow(source, i) >= 10) pff.emplace_back(i);
	      if (pff.empty()) goto BACKUP;
	      
	      // now do the face pruning: remove everything that is
	      // irrelevant for the next iteration
	      
	      // faces are still possibly relevant if they have a flow >=3 to
	      // a possible final face
	      // active = 3 -> possible final face
	      // active = 2 -> active face
	      // active = 1 -> passive face
	      // active = 0 -> transit face
	      // active = -1 -> irrelevant face	      
	      std::vector<int> active(fc, -1);
	      for (auto i = pff.begin(); i != pff.end(); ++i) active[*i] = 3;
	      for (std::size_t i = 0; i < fc; ++i) {
		if (active[i] != -1) continue;
		for (auto j = pff.begin(); j != pff.end(); ++j)
		  if (g.flow(i, *j) >= 3) { active[i] = 2; break; }
	      }
	      
	      // compute map face -> (some) halfedge on its boundary
	      std::vector<const HdsHalfedge*> fhedge(fc, 0);
	      for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i)
		if (fhedge[face[i->label]] == 0) fhedge[face[i->label]] = &*i;
	      
	      // mark passive faces, i.e., faces that have a vertex v
	      // from [cub,cue) on their boundary and a path P of
	      // length <= 2 in the dual to some active face so that P
	      // starts with an edge whose primal is not incident to v
	      for (std::size_t i = 0; i < fc; ++i) {
		if (fhedge[i] == 0) throw std::runtime_error("no edge for face");
		if (active[i] != -1) continue;
		auto e = fhedge[i];
		do {
		  std::size_t v = e->vertex->label;
		  if (v >= cub && v < cue)
		    for (auto f = e->next->next; f != e; f = f->next)
		      if (f->edge->ncr < 2 &&
			  f->vertex->label != v &&
			  f->twin->vertex->label != v)
			{
			  // do two step bfs from i via f
			  int fn = face[f->twin->label];
			  if (active[fn] >= 2) { active[i] = 1; break; }
			  for (auto j = dualg[fn].begin();
			       j != dualg[fn].end();
			       ++j)
			    if (active[*j] >= 2) { active[i] = 1; break; }
			}
		  e = e->next;
		} while (active[i] == -1 && e != fhedge[i]);
	      }
	      
	      // mark transit faces, i.e., faces adjacent to two
	      // (active||passive) faces
	      for (std::size_t i = 0; i < fc; ++i) {
		if (active[i] > 0) continue;
		int an = 0;
		int pn = 0;
		for (auto j = dualg[i].begin(); j != dualg[i].end(); ++j) 
		  if (active[*j] >= 2) ++an; else if (active[*j] == 1) ++pn;
		if (an >= 2 || (an == 1 && pn >= 1)) active[i] = 0;
	      }
	      
	      // determine relevant vertices
	      // ensure that the vertices cub...cue are mapped to 0...9 in relv
	      std::vector<const HdsHalfedge*> relv(10, 0); // relevant
	      std::vector<const HdsHalfedge*> prelv; // possibly relevant
	      for (auto i = d.vertices.begin(); i != d.vertices.end(); ++i) {
		auto e = i->halfedge;
		std::size_t nirf = 0; // #incident relevant faces
		auto a = e;
		do {
		  if (active[face[e->label]] >= 0) { ++nirf; a = e; }
		  e = e->next->twin;
		} while (e != i->halfedge);
		
		if (i->label >= cub && i->label < cue) {
		  // if no incident face is relevant, then this is
		  // actually not a valid solution
		  if (nirf == 0) goto BACKUP;
		  relv[(i->label) - cub] = a;
		} else if (nirf >= 2)
		  relv.push_back(a);
		else if (nirf == 1)
		  prelv.push_back(a);
	      }
	      // ... and crossings
	      for (auto i = d.crossings.begin(); i != d.crossings.end(); ++i) {
		auto e = i->halfedge;
		std::size_t nirf = 0;
		auto a = e;
		do {
		  if (active[face[e->label]] >= 0) { ++nirf; a = e; }
		  e = e->next->twin;
		} while (e != i->halfedge);
		if (nirf >= 2)
		  relv.push_back(a);
		else if (nirf == 1)
		  prelv.push_back(a);
	      }
	      
	      // at this point we know that this is a valid solution
	      is_extensible = true;
	      
	      // let's ensure that we don't contract triangles to lenses
	      std::vector<std::size_t> fsize(fc, 0); // size of faces
	      // map: vertex label in d -> vertex label in the reduced
	      // drawing nd
	      std::size_t total_nv = d.vertices.size() + d.crossings.size();
	      std::vector<std::size_t> vact(total_nv, total_nv);
	      // vact[] : total_nv == uninitialized 
	      std::size_t vic = 0;
	      for (auto i = relv.begin(); i != relv.end(); ++i) {
		auto j = *i;
		do {
		  ++fsize[face[j->label]];
		  j = j->next->twin;
		} while (j != *i);
		vact[(*i)->vertex->label] = vic++;
	      }
	      for (auto i = prelv.begin(); i != prelv.end(); ++i)
		if (fsize[face[(*i)->label]] < 3) {
		  ++fsize[face[(*i)->label]];
		  relv.push_back(*i);
		  vact[(*i)->vertex->label] = vic++;
		}
	      
	      // extract relevant parts of the drawing into a new drawing
	      Drawing<2> nd(relv.size() + 1);
	      std::deque<const HdsHalfedge*> todo(1, relv[0]);
	      // status of vertices: -1 -> not considered yet, 0 -> in queue
	      // and present in nd, 1 -> handled and all incident edges
	      // present in nd
	      std::vector<int> done(d.vertices.size() + d.crossings.size(), -1);
	      // Connected regions of irrelevant faces are merged into
	      // "black" regions. Record them here. 
	      std::vector<HdsHalfedge*> black;
	      while (!todo.empty()) {
		auto i = todo.front();
		todo.pop_front();
		std::size_t u = i->vertex->label;
		if (active[face[i->label]] < 0)
		  throw std::runtime_error("irrelevant face");
		// find a place to attach u to
		HdsHalfedge* last = nd.vertices[vact[u]].halfedge;
		if (nd.edges.empty()) {
		  // create first edge
		  auto p = prev_active(i, face, active, vact);
		  std::size_t v = p.first->vertex->label;
		  last = nd.add_first_edge(vact[u], vact[v], p.second)->twin;
		  if (last->vertex->label != 0)
		    throw std::runtime_error("wrong first edge in new drawing");
		  todo.push_back(p.first);
		  done[v] = 0;
		} else if (last == 0) {
		  todo.push_back(i);
		  continue;
		}
		
		// now act[u] is connected in nd and last points to
		// act[u]; build neighborhood of act[u] in nd
		done[u] = 1;
		std::size_t wnd = last->twin->vertex->label;
		// find wnd in d -> w
		std::size_t w;
		auto wi = i;
		do {
		  auto e = prev_active(wi, face, active, vact).first;
		  w = e->vertex->label;
		  if (vact[w] == wnd) break;
		  wi = wi->next->twin;
		} while (wi != i);
		if (vact[w] != wnd) throw std::runtime_error("no w");
		// go over edges incident to u in d
		for (auto j = wi->next->twin; j != wi; j = j->next->twin) {
		  if (active[face[j->label]] < 0 &&
		      active[face[j->twin->label]] < 0)
		    continue;
		  auto p = prev_active(j, face, active, vact);
		  std::size_t v = p.first->vertex->label;
		  if (last->next->vertex->label == vact[v]) {
		    // uv is already present in nd
		    last = last->next->twin;
		    continue;
		  }
		  
		  HdsPath path(2, last);
		  if (!nd.find_target(path, vact[v]))
		    throw std::runtime_error("cannot find v");
		  auto ne = nd.add_edge(path, vact[v], p.second);
		  if (active[face[j->label]] < 0)
		    black.push_back(ne->twin);
		  else if (active[face[j->twin->label]] < 0)
		    black.push_back(ne);
		  last = last->next->twin;
		  if (done[v] < 0) {
		    done[v] = 0;
		    todo.push_back(p.first);
		  }
		}
	      } // while (!todo.empty())
	      
	      // add an uncrossable star in each black face
	      std::vector<bool> blackdone(nd.halfedges.size(), false);
	      std::size_t nbf = 0; // #black faces on >=4 vertices
	      for (auto i = black.begin(); i != black.end(); ++i) {
		if (blackdone[(*i)->label]) continue;
		auto j = *i;
		std::size_t bc = 0;
		do {
		  blackdone[j->label] = true;
		  ++bc;
		  j = j->next;
		} while (j != *i);
		if (bc >= 4) {
		  // For all large black regions (that is, not
		  // triangles) we add an uncrossable star to ensure
		  // that they are not used anymore. We need an
		  // additional vertex (as a star center) for each
		  // these and have accounted for one of those
		  // only. As it turns out, we never need more than
		  // one for these graphs here, so let's not bother
		  // with handling more.
		  if (++nbf >= 2)
		    throw std::runtime_error("too many large black regions");
		  j = *i;
		  auto x = nd.add_edge(HdsPath({j,0}), nd.vertices.size()-1, 2);
		  for (;;) {
		    j = j->next->twin->next;
		    if (j == *i) break;
		    nd.add_edge(HdsPath({j,x}), nd.vertices.size()-1, 2);
		  } while (j != *i);
		}
	      }
	      
	      // is this a new solution?
	      for (auto x = solutions.begin(); x != solutions.end(); ++x) 
		if (are_isomorphic(*x, nd)) { 
		  std::cout << "--- Discard drawing #" << discarded
			    << ", isomorphic to solution #"
			    << (x - solutions.begin()) << std::endl;
		  std::ofstream of;
		  std::ostringstream filename;
		  filename << "discard-" << discarded << ".graphml";
		  of.open(filename.str());
		  d.graphml_output(of);
		  of.close();
		  store_marked(d, "discard-marked",
			       discarded, fc, vact, active, fhedge);
		  store_reduced(nd, "discard-reduced", discarded);
		  {
		    std::ofstream of;
		    std::ostringstream filename;
		    filename << "discard-iso-" << discarded << ".graphml";
		    of.open(filename.str());
		    x->graphml_output(of);
		    of.close();
		  }
		  ++discarded;
		  goto BACKUP;
		}
	      
	      // we really have a valid solution -> record it
	      if (!nd.is_valid()) throw std::runtime_error("nd is invalid");
	      solutions.push_back(nd);
	      std::cout << "Drawing #" << solutions.size()-1 << ":\n"
			<< d << std::endl;
	      if (pff.size() > 1)
		std::cout << "!!! Final face is not unique for this drawing ---"
			  << std::endl;
	      
	      {
		std::ofstream of;
		std::ostringstream filename;
		filename << "drawing-" << solutions.size()-1 << ".graphml";
		of.open(filename.str());
		d.graphml_output(of);
		of.close();
	      }
	      store_marked(d, "marked", solutions.size()-1,
			   fc, vact, active, fhedge);
	      store_reduced(nd, "reduced", solutions.size()-1);
	      
	      goto BACKUP;
	    } // if (++e == edges.end())
	  } // for (auto e = edges.begin();;)
	} // try block
	catch (done_w_current_drawing) {}
      
      if (constrained == 0)
	std::cout << "+++ Solution #" << solcount << " has "
		  << solutions.size() - solsofar
		  << " proper children" << std::endl;
      else if (fullsolcount > fullsolsofar)
	std::cout << "*** Solution #" << solcount
		  << " can be completed with a constrained cycle." << std::endl;
      else 
	std::cout << "--- Solution #" << solcount
		  << " cannot be completed with a constrained cycle."
		  << std::endl;
      
    } // for (int constrained = 0; constrained < 2; ++constrained)
    ++solcount;
  } // while (solcount < solutions.size())
  std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
  return 0;
}
