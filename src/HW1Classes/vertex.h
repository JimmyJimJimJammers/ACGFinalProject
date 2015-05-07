#ifndef _VERTEX_H
#define _VERTEX_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdlib>

// ==========================================================
// Stores the vertex position, used by the Mesh class

class Vertex {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Vertex(int i, const glm::vec3 &pos) : position(pos) { index = i; }
  
  // =========
  // ACCESSORS
  int getIndex() const { return index; }
  double x() const { return position.x; }
  double y() const { return position.y; }
  double z() const { return position.z; }
  const glm::vec3& getPos() const { return position; }
    const glm::vec3& getNorm() const {return normal; }
    const glm::vec3& getColor() const {return color; }
    const float getCrease() const {return crease; }

  // =========
  // MODIFIERS
  void setPos(glm::vec3 v) { position = v; }
    void setNorm(glm::vec3 n) {normal = n; }
    void setColor(glm::vec3 c) {color = c; }
    void setCrease(float c) {crease = c;}

private:

  // don't use these constructors
  Vertex() { assert(0); exit(0); }
  Vertex(const Vertex&) { assert(0); exit(0); }
  Vertex& operator=(const Vertex&) { assert(0); exit(0); }
  
  // ==============
  // REPRESENTATION
  glm::vec3 position;
    glm::vec3 normal;
    glm::vec3 color;

  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure, 
  // but we use it for hashing
  int index;
  float crease;

  // NOTE: the vertices don't know anything about adjacency.  In some
  // versions of this data structure they have a pointer to one of
  // their incoming edges.  However, this data is very complicated to
  // maintain during mesh manipulation, so it has been omitted.

};

// ==========================================================

#endif

