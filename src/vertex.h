#ifndef _VERTEX_H
#define _VERTEX_H

#include <glm/glm.hpp>
#include <stdio.h>

// ==========================================================

class Vertex {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Vertex(int i, const glm::vec3 &pos)
    {
        position = pos;
        index = i;
        s = 0;
        t = 0;
    }
    ~Vertex()
    {
        //printf("Deleted vert: %d\n", index);
    }
  
  // =========
  // ACCESSORS
  int getIndex() const { return index; }
  const glm::vec3& get() const { return position; }
  float get_s() const
  {
      assert(s >= 0);
      assert(s <= 1);
      //printf("About to return s(%f)\n", s);
      //checks for NaN
      assert(!(s != s));
      return s;
  }
  float get_t() const
  {
      assert(t >= 0);
      assert(t <= 1);
      //printf("About to return t(%f)\n", t);
      //checks for NaN
      assert(!(t != t));
      return t;
  }

  // =========
  // MODIFIERS
  void setTextureCoordinates(float _s, float _t)
    {
        assert(_s >= 0);
        assert(_s <= 1);
        assert(_t >= 0);
        assert(_t <= 1);
      s = _s; t = _t;
    }
    void set(glm::vec3 newPos){position = newPos;}

private:

  // ==============
  // REPRESENTATION
  glm::vec3 position;

  // texture coordinates
  // NOTE: arguably these should be stored at the faces of the mesh
  // rather than the vertices
  float s,t;

  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure
  int index;  

  // NOTE: the vertices don't know anything about adjacency.  In some
  // versions of this data structure they have a pointer to one of
  // their incoming edges.  However, this data is complicated to
  // maintain during mesh manipulation.
};

// ==========================================================

#endif

