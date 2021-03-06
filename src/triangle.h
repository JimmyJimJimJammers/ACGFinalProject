//
//  triangle.h
//  
//
//  Created by James McCarthy on 5/3/15.
//
//

#ifndef ____triangle__
#define ____triangle__

#include <stdio.h>
#include "edge.h"
#include "material.h"
#include "vertex.h"

class Edge;


class Triangle
{

public:
    // ========================
    // CONSTRUCTOR & DESTRUCTOR
    Triangle()
    {
        edge = NULL;
        id = next_triangle_id;
        next_triangle_id++;
        material = NULL;
        generation = 0;
        mate = NULL;
        leadsToBorder = false;
        subdivIndex = 0;
        //centerDisp = 0;
        
        //mateEdge = NULL;
        //mateVert1 = NULL;
        //mateVert2 = NULL;
    }//*/
    Triangle(Material *m)
    {
        edge = NULL;
        id = next_triangle_id;
        next_triangle_id++;
        material = m;
        generation = 0;
        mate = NULL;
        leadsToBorder = false;
        subdivIndex = 0;
        //centerDisp = 0;
        //mateEdge = NULL;
        //mateVert1 = NULL;
        //mateVert2 = NULL;
    }//*/

    // =========
    // ACCESSORS
    Vertex* getVert(int v);

    Vertex* operator[](int i) const;/* {
        assert (edge != NULL);
        if (i==0) return edge->getStartVertex();
        if (i==1) return edge->getNext()->getStartVertex();
        if (i==2) return edge->getNext()->getNext()->getStartVertex();
        assert(0); exit(0);
    }//*/
    Edge* getEdge() const
    {
        assert (edge != NULL);
        return edge;
    }//*/
    void setEdge(Edge *e)
    {
        assert (edge == NULL);
        edge = e;
    }//*/
    int getID() { return id; }
    
    Material* getMaterial()
    {
        return material;
    }

    void setMaterial(Material *m)
    {
        material = m;
    }
    
    int getGeneration()
    {
        return generation;
    }
    void setGeneration(int g)
    {
        generation = g;
    }
    void setGeneration(Triangle* m)
    {
        mate = m;
    }
    
    int getSubdivIndex()
    {
        return subdivIndex;
    }
    void setSubdivIndex(int index)
    {
        subdivIndex = index;
    }
    
    Triangle* getMate()
    {
        return mate;
    }
    void setMate(Triangle *sMate)
    {
        mate = sMate;
    }
    
    bool isBorderEdge()
    {
        return (mate == NULL && generation%2 == 1);
    }
    bool leadsToBorderTri()
    {
        if (isBorderEdge()){ leadsToBorder = true;}
        return leadsToBorder;
    }
    void setLeadsToBorder(bool b)
    {
        leadsToBorder = b;
    }
    
    float getCenterDisp();
    float getAvgDisp();
    
    //void setVertNull(int i);
    
    

    /*Vertex* getMateVert1()
    {
        return mateVert1;
    }
    Vertex* getMateVert2()
    {
        return mateVert2;
    }
    
    void setMateEdge(Vertex* a, Vertex* b)
    {
        mateVert1 = a;
        mateVert2 = b;
    }*/
    
    float getArea();
    glm::vec3 computeNormal() const;
    glm::vec3 computeCentroid() const;

    // NOTE: If you want to modify a triangle, it is recommended that
    // you remove it from the mesh, delete it, create a triangle object
    // with the changes, and re-add it.  This will ensure the edges get
    // updated appropriately.

protected:

    // don't use these constructors
    // Triangle(const Triangle &/*t*/) { assert(0); exit(0); }
    // Triangle& operator= (const Triangle &/*t*/) { assert(0); exit(0); }

    // ==============
    // REPRESENTATION
    Material *material;
    Edge *edge;
    int id;
    int generation;
    Triangle* mate;
    int subdivIndex;
    bool leadsToBorder;
    //float centerDisp;
    //Vertex* mateVert1;
    //Vertex* mateVert2;
    
    // triangles are indexed starting at 0
    static int next_triangle_id;
};

#endif /* defined(____triangle__) */
