#ifndef HDS_HPP
#define HDS_HPP

#include <vector>
#include <list>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <string>

// classic Halfedge Data Structure

struct HdsVertex;

struct HdsEdge;

struct HdsHalfedge {
  HdsVertex* vertex; // target vertex
  HdsHalfedge* prev; // previous halfedge along the face
  HdsHalfedge* next; // next halfedge along the face
  HdsHalfedge* twin; // twin halfedge
  HdsEdge* edge; // the underlying graph edge
  std::size_t label; // unique label for each halfedge
};

struct HdsVertex {
  HdsVertex(HdsHalfedge* h, std::size_t l) : halfedge(h), label(l) {}
  HdsHalfedge* halfedge; // some halfedge pointing to *this
  std::size_t label; // unique label for each vertex
};

// Creation paths:
//
// To be able to remove edges, we store how they were inserted. This
// is represented by a so-called creation path, which consists of
// between two and four halfedges: The first halfedge points to the
// source vertex and the insertion position in the rotation at that
// vertex; symmetrically, the last halfedge determines the target
// vertex and the position of the inserted edge in the rotation there;
// the zero, one, or two halfedges in between are crossed along the
// way from source to target.
//
// To insert an edge to an isolated vertex we use the convention that
// the last halfedge of the creation path is 0 (null pointer).
//
// Specifically, for a creation path p of length l, we have:
//
// (1) l >= 2;
// (2) either l == 2 and p[1] == 0 or
//   (a) p[0] and p[1] are on the same face of the drawing and p[1] is
//     not incident to the vertex that p[0] points to;
// [only for kplane >= 1:]
// (3) If l >= 3, then either l == 3 and p[2] == 0 or
//   (a) the underlying edge of p[1] has strictly fewer than kplane
//       crossings in the drawing;
//   (b) p[1]->twin and p[2] are on the same face of the drawing and
//       p[2] is not incident to the vertex that p[0] points to;
// [only for kplane >= 2:]
// (4) If l >= 4, then either l == 4 and p[3] == 0 or
//   (a) the underlying edge of p[2] has strictly fewer than kplane
//       crossings in the drawing;
//   (b) p[2]->twin and p[3] are on the same face of the drawing and
//       p[3] is not incident to the vertex that p[0] points to;
//   (c) the underlying edges of p[1] and p[2] are distinct.
//
// When building a possible creation path, we deal with objects that
// are not creation paths yet, but we want to express that, at least,
// they are valid up to this point. We call such an object a creation
// path stub; the conditions for it are the same as for a creation
// path, with one exception: the last halfedge may be incident to the
// vertex that p[0] points to.

typedef std::vector<HdsHalfedge*> HdsPath;

struct HdsEdge {
  HdsEdge(std::size_t u_, std::size_t v_,
	  const HdsPath& p, std::size_t c, std::size_t l_)
    : u(u_), v(v_), ncr(c), built(p), label(l_)
  {}
  std::size_t u, v; // indices of the two endpoints
  std::size_t ncr; // #crossings
  HdsPath built; // how the edge is built into the drawing
  std::size_t label; // unique label for each edge
};

// Most of the actual functionality is built into the Drawing
// class. It stores the containers of all objects and defines the
// interface to interact with them.
//
// Vertices are static and created upon initialization, although not
// yet connected by edges. Crossings may appear and disappear when
// edges are inserted and removed, respectively; hence, they are
// stored separately, in a list. Halfedges and edges are also stored
// in lists. The advantage of lists over, say, vectors/arrays is that
// pointers remain valid under insertions and removals of other
// elements. Insertion into a vector may lead to reallocation, which
// invalidates *all* pointers into it.
//
// The parameter kplane should be 0, 1 or 2.

template < int kplane >
struct Drawing {
  // n == #vertices
  Drawing(std::size_t n) {
    vertices.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
      vertices.emplace_back((HdsHalfedge*)(0), i);
  }

  // copy constructor; many pointers, we have to recompute all of them
  Drawing(const Drawing& d)
    : vertices(d.vertices),
      crossings(d.crossings),
      halfedges(d.halfedges),
      edges(d.edges)
  {
    // build indices
    std::vector<HdsHalfedge*> ind(halfedges.size());
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
      ind[i->label] = &*i;
    std::vector<HdsEdge*> ind_e(edges.size());
    for (auto i = edges.begin(); i != edges.end(); ++i) 
      ind_e[i->label] = &*i;
    std::vector<HdsVertex*> ind_c(crossings.size());
    for (auto i = crossings.begin(); i != crossings.end(); ++i) 
      ind_c[i->label - vertices.size()] = &*i;
    
    // fix pointers
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      if (vertices[i].label != i)
	throw std::runtime_error("wrong vertex label");
      if (vertices[i].halfedge != 0)
	vertices[i].halfedge = ind[vertices[i].halfedge->label];
    }
    std::size_t count = 0;
    for (auto i = crossings.begin(); i != crossings.end(); ++i,++count) {
      if (i->label < vertices.size())
	throw std::runtime_error("crossing label too small");
      if (i->label - vertices.size() != count)
	throw std::runtime_error("wrong crossing label");
      if (ind_c[i->label - vertices.size()] != &*i)
	throw std::runtime_error("pointer error");
      if (i->halfedge == 0)
	throw std::runtime_error("crossing without halfedge");
      i->halfedge = ind[i->halfedge->label];
    }
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
      if (i->vertex->label < vertices.size())
	i->vertex = &(vertices[i->vertex->label]);
      else 
	i->vertex = ind_c[i->vertex->label - vertices.size()];
      i->prev = ind[i->prev->label];
      i->next = ind[i->next->label];
      i->twin = ind[i->twin->label];
      i->edge = ind_e[i->edge->label];
    }
    for (auto i = edges.begin(); i != edges.end(); ++i) 
      for (auto x = i->built.begin(); x != i->built.end(); ++x)
	if (*x != 0) *x = ind[(*x)->label];
  }

  // no copy-assignment allowed
  Drawing& operator=(const Drawing&) = delete;
  
  // add x new vertices to the drawing
  // PRE: The drawing must not contain any crossings.
  void add_vertices(std::size_t x) {
    if (!crossings.empty())
      throw std::runtime_error("no crossings allowed when adding vertices");

    // As vertices are stored in a vector, inserting into it
    // potentially invalidates all vertex pointers. Thus, we have to
    // recompute them.
    
    // store vertex labels for halfedges
    std::vector<std::size_t> hind(halfedges.size());
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i)
      hind[i->label] = i->vertex->label;
    
    // expand the drawing
    std::size_t n = vertices.size();
    vertices.reserve(n+x);
    for (std::size_t i = 0; i < x; ++i)
      vertices.emplace_back((HdsHalfedge*)(0), n+i);

    // build index
    std::vector<HdsVertex*> vind(vertices.size(), 0);
    for (auto i = vertices.begin(); i != vertices.end(); ++i)
      vind[i->label] = &*i;
    
    // consolidate pointers
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
      i->vertex = vind[hind[i->label]];
  }

  // create and return a new pair of twin halfedges that correspond to
  // the edge e
  HdsHalfedge* new_twin(HdsEdge* e) {
    halfedges.emplace_back();
    halfedges.back().label = halfedges.size()-1;
    HdsHalfedge* nhe = &(halfedges.back());
    halfedges.emplace_back();
    halfedges.back().label = halfedges.size()-1;
    HdsHalfedge* nhet = &(halfedges.back());
    nhe->edge = nhet->edge = e;
    nhe->twin = nhet;
    return nhet->twin = nhe;
  }

  // split h into a path of two halfedges by inserting a crossing on
  // it, so that h points to this new crossing. The underlying edge of
  // the graph remains the same (for both halfedges along the path).
  void split_edge(HdsHalfedge* h) {
    crossings.emplace_back(h, vertices.size() + crossings.size());
    HdsVertex* nv = &(crossings.back());
    HdsHalfedge* nhe = new_twin(h->edge);
    HdsHalfedge* nhet = nhe->twin;
    nhe->vertex = h->vertex;
    nhe->prev = h;
    nhet->vertex = nv;
    nhet->next = h->twin;
    if (h->vertex->halfedge == h) h->vertex->halfedge = nhe;
    h->vertex = nv;
    nv->halfedge = h;
    if (h->next == h->twin) {
      nhe->next = nhet;
      nhet->prev = nhe;
    } else {
      nhe->next = h->next;
      nhe->next->prev = nhe;
      nhet->prev = h->twin->prev;
      nhet->prev->next = nhet;
    }
    nhe->prev->next = nhe;
    nhet->next->prev = nhet;
    if (++(h->edge->ncr) > kplane)
      throw std::runtime_error("too many crossings per edge");
  }

  // opposite of split_edge above: contract the path of two halfedges
  // formed by h and h->next.
  // PRE: h points to the last crossing created.
  void join_vertex(HdsHalfedge* h) {
    if (h->edge->ncr == 0) throw std::runtime_error("internal error");
    if (h->vertex != &(crossings.back()))
      throw std::runtime_error("only the last crossing can be removed");
    HdsHalfedge* nhe = h->next;
    HdsHalfedge* nhet = nhe->twin;
    if (nhe != &(halfedges.back()) && nhet != &(halfedges.back()))
      throw std::runtime_error("only the last halfedge can be removed");
    h->vertex = nhe->vertex;
    if (nhe->vertex->halfedge == nhe) nhe->vertex->halfedge = h;
    if (nhe->next == nhet) {
      h->next = h->twin;
      h->twin->prev = h;
    } else {
      h->next = nhe->next;
      h->next->prev = h;
      h->twin->prev = nhet->prev;
      h->twin->prev->next = h->twin;
    }
    --(h->edge->ncr);
    crossings.pop_back();
    halfedges.pop_back();
    halfedges.pop_back();
  }

private:

  // helper function to insert edges: add and return a new halfedge
  // stub that corresponds to the edge e, and add it in the rotation
  // of s->vertex right behind s.
  //
  // It is a "stub" because it is only properly connected to the rest
  // of the HDS at its source s, the other endpoint must be addressed
  // by the calling function.
  HdsHalfedge* create_starting_stub_at(HdsHalfedge* s, HdsEdge* e) {
    // create new edge twin
    HdsHalfedge* nhe = new_twin(e);
    HdsHalfedge* nhet = nhe->twin;
    // insert at source
    nhet->vertex = s->vertex;
    nhet->next = s->next;
    nhet->next->prev = nhet;
    nhe->prev = s;
    return nhe->prev->next = nhe; 
  }

  // helper function to remove edges and the inverse of
  // create_starting_stub_at above: detach the edge (stub) right
  // behind s in the rotation of s->vertex
  void detach_stub(HdsHalfedge* s) {
    s->next = s->next->twin->next;
    s->next->prev = s;
  }

public:
  // add a first edge to the drawing, with pcr prescribed (artifical)
  // crossings; return the halfedge that points to v
  HdsHalfedge* add_first_edge(std::size_t u,
			      std::size_t v,
			      std::size_t pcr = 0)
  {
    if (!edges.empty())
      throw std::runtime_error("inserting another first edge");
    // create new edge twin (the creation path for this edge will be
    // built manually once the two halfedges have been created)
    edges.emplace_back(u, v, HdsPath(), pcr, 0);
    HdsHalfedge* nhe = new_twin(&(edges.back()));
    HdsHalfedge* nhet = nhe->twin;
    nhe->prev = nhe->next = nhet;
    nhet->prev = nhet->next = nhe;
    nhe->vertex = &(vertices[u]);
    nhet->vertex = &(vertices[v]);
    nhe->vertex->halfedge = nhe;
    nhet->vertex->halfedge = nhet;
    edges.back().built.emplace_back(nhe);
    edges.back().built.emplace_back(nhet);
    return nhet;
  }

  // add an edge described by p (which determines source and potential
  // crossings) with target v to the drawing, with pcr prescribed
  // (artifical) crossings; return the halfedge that points to v
  HdsHalfedge* add_edge(const HdsPath& p,
			std::size_t v,
			std::size_t pcr = 0)
  {
    if (p.size() > 2 + pcr + kplane)
      throw std::runtime_error("inserting an edge with too many crossings");
    if (p.size() < 2) throw std::runtime_error("edge not properly described");
    //std::cerr << "add: " << p[0]->vertex->label << ">" << v << " (cr = "
    // << pcr << ")" << std::endl;
    edges.emplace_back(p[0]->vertex->label, v, p, pcr, edges.size());
    HdsHalfedge* nhe = create_starting_stub_at(p[0], &(edges.back()));
    for (std::size_t i = 1; i < p.size(); ++i) {
      HdsHalfedge* nhet = nhe->twin;
      if (i < p.size()-1) {
	split_edge(p[i]);
	++(edges.back().ncr);
      } else if (p[i] == 0) {
	// this is a new vertex
	if (vertices[v].halfedge != 0)
	  throw std::runtime_error("vertex is already present");
	nhe->next = nhet;
	nhet->prev = nhe;
	nhe->vertex = &(vertices[v]);
	nhe->vertex->halfedge = nhe;
	edges.back().built.back() = nhe;
	break;
      }
      nhe->next = p[i]->next;
      nhe->next->prev = nhe; 
      nhet->prev = p[i];
      nhet->prev->next = nhet; 
      nhe->vertex = p[i]->vertex;
      if (i < p.size()-1) 
	nhe = create_starting_stub_at(nhe->next->twin, &(edges.back()));
    }
    return nhe;
  }

  // remove the last edge (i.e., the edge that was inserted last)
  void remove_edge() {
    HdsPath p = edges.back().built;
    if (p.size() < 2) throw std::runtime_error("no creation path");
    if (p[0]->next->edge != &(edges.back()))
      throw std::runtime_error("creation path mismatch");
    HdsVertex* v = p.back()->vertex;
    if (v->halfedge->next == v->halfedge->twin) 
      // this is the only edge incident to v
      v->halfedge = 0;
    else 
      detach_stub(p.back());
    for (std::size_t i = p.size()-2;; --i) {
      if (i == 0) {
	detach_stub(p[i]);
	halfedges.pop_back();
	halfedges.pop_back();
	break;
      }
      detach_stub(p[i]->next->twin->next->twin);
      halfedges.pop_back();
      halfedges.pop_back();
      detach_stub(p[i]);
      join_vertex(p[i]);
    }
    edges.pop_back();
  }

  // Try to extend the creation path stub p by walking p.back() around
  // its face using ->next, so as to reach v; return true iff
  // successful (p.back() is updated along the way, regardless).
  // PRE: p.size() >= 2 and p is a creation path stub.
  bool find_target(HdsPath& p, std::size_t v) {
    assert(p.size() >= 2);
    if (vertices[v].halfedge == 0) {
      // new vertex, we can just insert it into this face
      p.back() = 0;
      return true;
    }
    // walk p.back() around its face
    std::size_t i = p.size()-1;
    HdsHalfedge* end = (i == 1 ? p[0] : p[i-1]->twin); 
    for (;;) {
      p[i] = p[i]->next;
      // must not self-cross
      if (i == 2) {
	// The following test is relevant (only) if p[2] is on the
	// same face as p[0], which can happen if the current
	// graph/drawing has cut edges. It is not in the current
	// drawing, but if we follow through with the creation path p,
	// then p[0] and p[1] will be connected by an edge. Hence, for
	// the remainder of p, we have to consider this connection as
	// a face boundary.
	if (p[2] == p[1])
	  p[2] = p[0];
	else if (p[2] == p[0])
	  p[2] = p[1];
      } else if (i == 3) {
	// Similar to above, we have two connections that are not
	// present in the drawing but have to be considered as face
	// boundaries for the purposes of finishing p:
	// (1) a connection between p[0] and p[1] (same as above,
	//     corresponds to the edge between the starting vertex and
	//     the first crossing along p) and
	// (2) a connection between p[1]->twin and p[2] (which
	//     corresponds to the edge between the first and the
	//     second crossing along p)
	if (p[3] == p[2])
	  p[3] = p[1]->twin;
	else if (p[3] == p[1]->twin)
	  p[3] = p[2];
	else if (p[3] == p[1])
	  p[3] = p[0];
	else if (p[3] == p[0])
	  p[3] = p[1];
      }
      if (p[i] == end) return false;
      if (p[i]->vertex->label == v) return true;
    }
    return false; // avoid compiler warning
  }
  
private:
  // try to build a creation path that extends the creation path stub
  // p[0],p[1],p[2] to reach v, using a crossing at p[2]. To find a
  // suitable crossing, advance p[2] using ->next and stop as soon as
  // it reaches p[1]->twin. In particular, this will change p[2]
  // accordingly, and may change p[3] as well, but not p[0] and p[1].
  // Return true iff a creation path was found. 
  //
  // PRE: (1) p.size() >= 4 and (2) p[0],p[1],[2] is a creation path
  //   stub. (no assumption concerning p[3], it will be appropriately
  //   initialized inside this function)
  bool find_second_crossing(HdsPath& p, std::size_t v) {
    assert(p.size() >= 4);
    std::size_t u = p[0]->vertex->label;
    for (;;) {
      p[2] = p[2]->next;
      // must not self-cross; see comment in find_target(...) above
      if (p[2] == p[1])
	p[2] = p[0];
      else if (p[2] == p[0])
	p[2] = p[1];
      if (p[2] == p[1]->twin) break;
      // simple 2-plane drawings only
      if (p[2]->edge->u == u || p[2]->edge->u == v ||
	  p[2]->edge->v == u || p[2]->edge->v == v ||
	  p[2]->edge->label == p[1]->edge->label ||
	  p[2]->edge->ncr >= kplane)
	continue;
      p[3] = p[2]->twin;
      if (find_target(p, v)) return true;
    }
    return false;
  }

  // try to build a creation path that extends the creation path stub
  // p[0],p[1] to reach v, using a crossing at p[1] (and possibly
  // another crossing at p[2], in case that p.size() == 4). To find a
  // suitable crossing, advance p[1] using ->next and stop as soon as
  // it reaches p[0]. In particular, this will change p[1]
  // accordingly, and may change p[2] as well, but not p[0].
  // Return true iff a creation path was found. 
  //
  // PRE: (1) p.size() >= 3 and (2) p[0],p[1] is a creation path stub.
  //   (no assumption concerning p[2], it will be appropriately
  //   initialized inside this function)
  bool find_first_crossing(HdsPath& p, std::size_t v) {
    assert(p.size() >= 3);
    std::size_t u = p[0]->vertex->label;
    for (;;) {
      p[1] = p[1]->next;
      if (p[1] == p[0]) break;
      // simple k-plane drawings only
      if (p[1]->edge->u == u || p[1]->edge->u == v ||
	  p[1]->edge->v == u || p[1]->edge->v == v ||
	  p[1]->edge->ncr >= kplane)
	continue;
      p[2] = p[1]->twin;
      if (p.size() == 3 && find_target(p, v)) return true;
      if (p.size() == 4 && find_second_crossing(p, v)) return true;
    }
    return false;
  }

public:
  // try to build a first creation path for an edge from u to v with c
  // prescribed (artificial) crossings; return the path, or an empty
  // path if no such creation path exists.
  // PRE: u has degree at least one in the drawing
  HdsPath first_path(std::size_t u, std::size_t v, std::size_t pcr = 0) {
    HdsHalfedge* end = vertices[u].halfedge;
    if (end == 0)
      throw std::runtime_error("don't start an edge from an isolated vertex");
    HdsPath p;
    p.emplace_back(end);
    p.emplace_back(end);
    // can we get there without any crossing? Try all positions in the
    // rotation of u.
    do {
      if (find_target(p, v)) return p;
      p[1] = p[0] = p[0]->twin->prev;
    } while (p[0] != end);
    // no? -> try the same with exactly one crossing
    if (pcr < 2 && kplane >= 1) {
      p.emplace_back();
      do {
	if (find_first_crossing(p, v)) return p;
	p[1] = p[0] = p[0]->twin->prev;
      } while (p[0] != end);
      // no? -> try the same with two crossings
      if (pcr < 1 && kplane >= 2) {
	p.emplace_back();
	do {
	  if (find_first_crossing(p, v)) return p;
	  p[1] = p[0] = p[0]->twin->prev;
	} while (p[0] != end);
      }
    }
    // no suitable creation path exists
    return HdsPath();
  }

  // try to find another creation path that extends the creation path
  // stub p to reach v, with pcr prescribed (artificial) crossings on
  // this edge. Return true iff another path was found.
  // PRE: (1) p.size() >= 2 and (2) p[0],...,p[p.size()-2] is a
  //      creation path stub.
  //
  // We enumerate creation paths the following way:
  // (1) shorter paths (i.e., with fewer crossings) before longer
  //     paths;
  // (2) explore the possibilities along each face boundary using
  //     ->next, until reaching the starting point again.
  bool next_path(HdsPath& p, std::size_t v, std::size_t pcr = 0) {
    assert(p.size() >= 2);
    HdsHalfedge* end = p[0]->vertex->halfedge;
    bool v_is_isolated = (vertices[v].halfedge == 0);
    if (p.size() == 4) {
      // first try to find another option to access v on the same face
      if (!v_is_isolated && find_target(p, v)) return true;
      // no such option -> take back last crossing
      if (find_second_crossing(p, v)) return true;
      // still no option -> take back both crossings
      do {
	if (find_first_crossing(p, v)) return true;
	// as a last resort, try a new starting edge
	p[1] = p[0] = p[0]->twin->prev;
      } while (p[0] != end);
      return false;
    } // if (p.size() == 4)
    else if (p.size() == 3) {
      // first try to find another option to access v on the same face
      if (!v_is_isolated && find_target(p, v)) return true;
      // no such option -> take back the crossing
      do {
	if (find_first_crossing(p, v)) return true;
	// no luck -> as a last resort, try a new starting edge
	p[1] = p[0] = p[0]->twin->prev;
      } while (p[0] != end);
      // no further option with one crossing -> start to explore options with two
      if (pcr + 1 < kplane) {
	p.emplace_back();
	do {
	  if (find_first_crossing(p, v)) return true;
	  p[1] = p[0] = p[0]->twin->prev;
	} while (p[0] != end);
      }
    } // if (p.size() == 3)
    else { // p.size() == 2
      // first try to find another option to access v on the same face
      if (!v_is_isolated && find_target(p, v)) return true;
      // no such option -> try a new starting edge
      for (;;) {
	p[1] = p[0] = p[0]->twin->prev;
	if (p[0] == end) break;
	if (find_target(p, v)) return true;
      }
      // no such option -> try with one crossing
      if (pcr < kplane) {
	p.emplace_back();
	do {
	  if (find_first_crossing(p, v)) return true;
	  p[1] = p[0] = p[0]->twin->prev;
	} while (p[0] != end);
	// no such option -> try with two crossings
	if (pcr + 1 < kplane) {
	  p.emplace_back();
	  do {
	    if (find_first_crossing(p, v)) return true;
	    p[1] = p[0] = p[0]->twin->prev;
	  } while (p[0] != end);
	}
      }
    } // (p.size() == 2)
    // out of options
    return false;
  }

  // debug output of all halfedges
  void dump() const {
    std::cerr << "{ #vertices = " << vertices.size()
	      << ", #edges = " << edges.size()
	      << ", #halfedges = " << halfedges.size()
	      << ", #crossings = " << crossings.size() << "\n";
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
      std::cerr << i->label << "\n"
		<< "\tprev = " << i->prev->label << "\n"
		<< "\tnext = " << i->next->label << "\n"
		<< "\ttwin = " << i->twin->label << "\n"
		<< "\tvert = " << i->vertex->label << "\n"
		<< "\tedge = " << i->edge->u << "-" << i->edge->v << "\n";
    std::size_t c = 0;
    for (auto i = edges.begin(); i != edges.end(); ++i, ++c)
      std::cerr << "e" << c << " " << i->u << " - " << i->v << "\n";
    std::cerr << "}\n";
  }

  // output a human readable recipe to draw the drawing, edge by edge
  void construction() {
    for (auto i = edges.begin(); i != edges.end(); ++i) {
      HdsPath& p = i->built;
      std::cerr << "draw an edge leaving " << i->u
		<< " cw after the edge to "
		<< (p[0]->edge->u == i->u ? p[0]->edge->v : p[0]->edge->u);
      for (std::size_t j = 1; j < p.size()-1; ++j) {
	// determine crossing orientation
	HdsHalfedge* src = p[j];
	while (src->vertex->label >= vertices.size())
	  src = src->next->twin->next;
	HdsHalfedge* tgt = p[j]->twin;
	while (tgt->vertex->label >= vertices.size())
	  tgt = tgt->next->twin->next;
	std::cerr << "\n\tacross the edge " << src->vertex->label
		  << "-" << tgt->vertex->label;
      }
      std::cerr << "\n\tto " << i->v << "\n";
    }
  }

  // graphml output of the drawing, to be viewed with yEd, for instance
  std::ostream& graphml_output(std::ostream& o) const {
    std::size_t total_nv = vertices.size() + crossings.size();
    std::vector<std::size_t> vcol(total_nv, total_nv);
    for (std::size_t i = 0; i < vertices.size(); ++i) vcol[i] = 0;
    return graphml_output(o, vcol);
  }
  
  std::ostream& graphml_output(std::ostream& o,
			       const std::vector<std::size_t>& vcol) const
  { return graphml_output_footer(graphml_output_header(o, vcol)); }
  
  std::ostream& graphml_output_header(std::ostream& o,
				      const std::vector<std::size_t>& vcol) const
  {
    o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
      << "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
      << "    xmlns:y=\"http://www.yworks.com/xml/graphml\"\n"
      << "    xsi:schemaLocation"
      << "=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
      << "  <key for=\"node\" id=\"d5\" yfiles.type=\"nodegraphics\"/>\n"
      << "  <key for=\"edge\" id=\"d9\" yfiles.type=\"edgegraphics\"/>\n"
      << "  <graph id=\"G\" edgedefault=\"undirected\">\n";
    for (std::size_t i = 0; i < vertices.size(); ++i)
      o << "    <node id=\"v" << i << "\">\n"
	<< "      <data key=\"d5\">\n"
        << "        <y:ShapeNode>\n"
	<< "          <y:NodeLabel>" << i << "</y:NodeLabel>\n"
	<< "          <y:Geometry height=\"30.0\" width=\"30.0\"/>\n"
	<< "          <y:Fill color=\"#"
	<< (vcol[i] < vcol.size() ? "99CCFF" : "C0C0C0")
	<< "\"/>\n"
	<< "          <y:Shape type=\"ellipse\"/>\n"
        << "        </y:ShapeNode>\n"
	<< "      </data>\n"
	<< "    </node>\n";
    for (std::size_t i = 0; i < crossings.size(); ++i)
      o << "    <node id=\"c" << i << "\">\n"
	<< "      <data key=\"d5\">\n"
        << "        <y:ShapeNode>\n"
	<< "          <y:NodeLabel>" << i << "</y:NodeLabel>\n"
	<< "          <y:Geometry height=\"15.0\" width=\"15.0\"/>\n"
	<< "          <y:Fill color=\"#"
	<< (vcol[i+vertices.size()] < vcol.size() ? "FF6600" : "C0C0C0")
	<< "\"/>\n"
	<< "          <y:Shape type=\"rectangle\"/>\n"
        << "        </y:ShapeNode>\n"
	<< "      </data>\n"
	<< "    </node>\n";
    std::vector<std::size_t> edgeparts(edges.size(), 0);
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
      if (i->label > i->twin->label) ++edgeparts[i->edge->label];
    
    std::size_t ec = 0;
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
      if (i->label > i->twin->label) continue;
      std::size_t u = i->vertex->label;
      std::size_t v = i->twin->vertex->label;
      std::size_t cr = edgeparts[i->edge->label] - 1;
      o << "    <edge id=\"e" << ec++ << "\" source=\"";
      if (u < vertices.size())
	o << "v" << u; else o << "c" << u - vertices.size();
      o << "\" target=\"";
      if (v < vertices.size())
	o << "v" << v; else o << "c" << v - vertices.size();
      o << "\"";
      std::string color = "3366FF";
      std::string width = "1.5";
      if (cr < 2) {
	width = "2";
	if (cr == 0) {
	  if (i->edge->ncr == 2) {
	    color = "339966";
	    width = "3";
	  } else if (i->edge->ncr == 1) {
	    color = "993366";
	  } else
	    color = "99CC00";
	} else 
	  color = "993366";
      }
      o << ">\n"
	<< "      <data key=\"d9\">\n"
	<< "        <y:PolyLineEdge>\n"
	<< "          <y:LineStyle color=\"#" << color << "\" "
	<<                "type=\"line\" width=\"" << width << "\"/>\n"
	<< "          <y:BendStyle smoothed=\"true\"/>\n"
	<< "        </y:PolyLineEdge>\n"
	<< "      </data>\n"
	<< "    </edge>\n";
    }
    return o;
  }
  
  std::ostream& graphml_output_footer(std::ostream& o) const {
    return o << "  </graph>\n</graphml>\n";
  }
  
  bool is_valid() const {
    // build index
    std::vector<const HdsHalfedge*> ind(halfedges.size());
    std::size_t count = 0;
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
      ind[i->label] = &*i;
      if (i->label != count++) return false;
    }
    std::vector<const HdsEdge*> ind_e(edges.size());
    count = 0;
    for (auto i = edges.begin(); i != edges.end(); ++i) {
      ind_e[i->label] = &*i;
      if (i->label != count++) return false;
    }
    std::vector<const HdsVertex*> ind_c(crossings.size());
    count = 0;
    for (auto i = crossings.begin(); i != crossings.end(); ++i) {
      ind_c[i->label - vertices.size()] = &*i;
      if (i->label - vertices.size() != count++) return false;
    }
    count = 0;
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
      if (i->vertex->label < vertices.size() &&
	  i->vertex != &(vertices[i->vertex->label]))
	return false;
      if (i->vertex->label >= vertices.size() &&
	  i->vertex != ind_c[i->vertex->label - vertices.size()])
	return false;
      if (i->prev != ind[i->prev->label]) return false;
      if (i->next != ind[i->next->label]) return false;
      if (i->twin != ind[i->twin->label]) return false;
      if (i->edge != ind_e[i->edge->label]) return false;
      if (i->next->prev != &*i) return false;
      if (i->prev->next != &*i) return false;
      if (i->twin->twin != &*i) return false;
    }
    count = 0;
    for (auto i = vertices.begin(); i != vertices.end(); ++i,++count) {
      if (i->label != count) return false;
      //if (i->halfedge == 0) return false;
      if (i->halfedge != 0) {
	if (i->halfedge->vertex->label != i->label) return false;
	if (i->halfedge != ind[i->halfedge->label]) return false;
      }
    }
    for (auto i = crossings.begin(); i != crossings.end(); ++i) {
      if (i->halfedge == 0 || i->halfedge->vertex->label != i->label)
	return false;
      if (i->halfedge != ind[i->halfedge->label]) return false;
    }

    // check that halfedges correctly implement their underlying edges
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
      std::size_t u = i->edge->u;
      if (i->vertex->label != u) continue;
      std::size_t v = i->edge->v;
      std::size_t iel = i->edge->label;
      auto e = i->twin;
      if (e->edge->label != iel) return false;
      if (e->vertex->label == v) continue;
      if (e->vertex->label < vertices.size()) return false;
      e = e->next->twin->next;
      if (e->edge->label != iel) return false;
      if (e->vertex->label == v) continue;
      if (e->vertex->label < vertices.size()) return false;
      e = e->next->twin->next;
      if (e->edge->label != iel) return false;
      if (e->vertex->label != v) return false;
    }
    return true;
  }
  
  std::vector<HdsVertex> vertices;
  std::list<HdsVertex> crossings;
  std::list<HdsHalfedge> halfedges;
  std::list<HdsEdge> edges;
};

// text output of rotation system
template < int kplane >
std::ostream& operator<<(std::ostream& o, const Drawing<kplane>& d) {
  o << "{ #vertices = " << d.vertices.size()
    << ", #edges = " << d.edges.size()
    << ", #halfedges = " << d.halfedges.size()
    << ", #crossings = " << d.crossings.size() << "\n";
  for (auto i = d.vertices.begin(); i != d.vertices.end(); ++i) {
    o << "\t" << i->label << " : ";
    HdsHalfedge* e = i->halfedge;
    if (e != 0)
      do {
	std::size_t target = (e->edge->u == i->label ? e->edge->v : e->edge->u);
	o << target;
	HdsHalfedge* f = e->twin;
	if (f->vertex->label != target) {
	  o << "(";
	  for (; f->vertex->label != target; f = f->next->twin->next) {
	    // determine crossing orientation
	    HdsHalfedge* g = f->next;
	    std::size_t c1 = g->edge->u;
	    std::size_t c2 = g->edge->v;
	    for (;;) {
	      if (g->vertex->label == c1) break;
	      if (g->vertex->label == c2) {
		std::swap(c1, c2);
		break;
	      }
	      g = g->next->twin->next;
	    }
	    o << "x" << c1 << "-" << c2;
	  }
	  o << ")";
	}
	o << " ";
	e = e->twin->prev;
      } while (e != i->halfedge);
    o << "\n";
  }
  return o << "}\n";
}

#endif
