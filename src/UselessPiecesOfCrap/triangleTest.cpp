//
//  triangleTest.cpp
//  
//
//  Created by James McCarthy on 5/3/15.
//
//
#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include <cstdlib>
#include "edge.h"
#include "material.h"
#include "triangleTest.h"
#include "vertex.h"


//class Edge;
class Triangle;

// ===========================================================
// Stores the indices to the 3 vertices of the triangles,
// used by the mesh class

// ========================
// CONSTRUCTOR & DESTRUCTOR
/*Triangle::Triangle() {
 edge = NULL;
 id = next_triangle_id;
 next_triangle_id++;
 material = NULL;
 }
 Triangle::Triangle(Material *m) {
 edge = NULL;
 id = next_triangle_id;
 next_triangle_id++;
 material = m;
 }*/

// =========
// ACCESSORS
/*Vertex* Triangle::operator[](int i) const {
 assert (edge != NULL);
 if (i==0) return edge->getStartVertex();
 if (i==1) return edge->getNext()->getStartVertex();
 if (i==2) return edge->getNext()->getNext()->getStartVertex();
 assert(0); exit(0);
 }*/
Triangle::Triangle(Material *m)
{
    edge = NULL;
    id = next_triangle_id;
    next_triangle_id++;
    material = m;
}

Vertex* Triangle::getVert(int v)
{
    if (v==0) return edge->getStartVertex();
    if (v==1) return edge->getNext()->getStartVertex();
    if (v==2) return edge->getNext()->getNext()->getStartVertex();
    else
        return NULL;
}
/*Edge* Triangle::getEdge() const
 {
 assert (edge != NULL);
 return edge;
 }//*/
/*void Triangle::setEdge(Edge *e)
 {
 assert (edge == NULL);
 edge = e;
 }*/
//int Triangle::getID() { return id; }


// don't use these constructors
//Triangle::Triangle(const Triangle &/*t*/) { assert(0); exit(0); }
//Triangle::Triangle& operator= (const Triangle &/*t*/) { assert(0); exit(0); }

#endif
