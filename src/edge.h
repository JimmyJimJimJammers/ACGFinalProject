#ifndef _EDGE_H_
#define _EDGE_H_

#include <cassert>
#include <cstdlib>
#include "triangle.h"

class Vertex;
class Face;
class Triangle;

// ===================================================================
// half-edge data structure

class Edge { 

public:

  // ========================
  // CONSTRUCTORS & DESTRUCTOR
  Edge(Vertex *vs, Vertex *ve, Face *f);
    Edge(Vertex *vs, Vertex *ve, Triangle *t);//NEW*****
  ~Edge();

  // =========
  // ACCESSORS
  Vertex* getStartVertex() const { assert (start_vertex != NULL); return start_vertex; }
  Vertex* getEndVertex() const
  {
      //printf("SHERE8.111\n");
      assert (end_vertex != NULL);
      //printf("SHERE8.112\n");
      return end_vertex;
  }
  Edge* getNext() const { assert (next != NULL); return next; }
  Face* getFace() const { assert (face != NULL); return face; }
    Triangle* getTriangle() const { assert (triangle != NULL); return triangle; }//NEW*****
  Edge* getOpposite() const
  {
    // warning!  the opposite edge might be NULL!
    return opposite;
  }
  float Length() const;
    Edge* getTriangleEdge()
    {
        return triangleEdge;
    }
    Edge* getTriangleEdgeSubdiv()
    {
        return triangleEdgeSubdiv;
    }

  // =========
  // MODIFIERS
  void setOpposite(Edge *e)
  {
    assert (opposite == NULL); 
    assert (e != NULL);
    assert (e->opposite == NULL);
    opposite = e; 
    e->opposite = this; 
  }
  void clearOpposite() { 
    if (opposite == NULL) return; 
    assert (opposite->opposite == this); 
    opposite->opposite = NULL;
    opposite = NULL; 
  }
  void setNext(Edge *e) {
    assert (next == NULL);
    assert (e != NULL);
    assert (face == e->face);
    next = e;
  }
    void setTriangleEdge(Edge *e)
    {
        triangleEdge = e;
    }
    void setTriangleEdgeSubdiv(Edge *e)
    {
        triangleEdgeSubdiv = e;
    }
    
    void setStartVertex(Vertex* start)
    {
        assert(start != NULL);
        start_vertex = start;
    }
    void setEndVertex(Vertex* end)
    {
        assert(end != NULL);
        end_vertex = end;
    }

private:

  Edge(const Edge&) { assert(0); }
  Edge& operator=(const Edge&) { assert(0); exit(0); }

  // ==============
  // REPRESENTATION
  // in the half edge data adjacency data structure, the edge stores everything!
  // note: it's technically not necessary to store both vertices, but it makes
  //   dealing with non-closed meshes easier
  Vertex *start_vertex;
  Vertex *end_vertex;
    Triangle *triangle;//NEW**********
    Edge *triangleEdge;//used for converting quads to tris, and seeing if verts have already been created
    Edge *triangleEdgeSubdiv;//used for converting quads to tris, and seeing if verts have already been created
  Face *face;
  Edge *opposite;
  Edge *next;
};

// ===================================================================

#endif

