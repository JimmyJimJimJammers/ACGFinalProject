//
//  triangle.cpp
//  
//
//  Created by James McCarthy on 5/3/15.
//
//

#include "triangle.h"



Vertex* Triangle::operator[](int i) const
{
    assert (edge != NULL);
    if (i==0) return edge->getStartVertex();
    if (i==1) return edge->getNext()->getStartVertex();
    if (i==2) return edge->getNext()->getNext()->getStartVertex();
    assert(0); exit(0);
}


Vertex* Triangle::getVert(int v)
{
    if (v==0) return edge->getStartVertex();
    if (v==1) return edge->getNext()->getStartVertex();
    if (v==2) return edge->getNext()->getNext()->getStartVertex();
    else
        return NULL;
}

glm::vec3 Triangle::computeNormal() const
{
    // note: this face might be non-planar, so average the two triangle normals
    glm::vec3 a = (*this)[0]->get();
    glm::vec3 b = (*this)[1]->get();
    glm::vec3 c = (*this)[2]->get();
    //printf("Set the position vectors to compute normal\n");
    glm::vec3 v12 = b - a;
    glm::vec3 v23 = c - b;
    //printf("Did some math\n");
    glm::vec3 normal = glm::cross(v12,v23);
    //std::cout << "\t\t\tFaceNormal: " << normal.x << ", " << normal.y << ", " << normal.z << "\n";
    normal = glm::normalize(normal);
    //printf("Finished computing\n");
    //std::cout << "\t\t\tFaceNormal: " << normal.x << ", " << normal.y << ", " << normal.z << "\n";
    return normal;
}

glm::vec3 Triangle::computeCentroid() const {
    return (1.0f/3.0f) * ((*this)[0]->get() +
                          (*this)[1]->get() +
                          (*this)[2]->get());
}

float Triangle::getArea()
{
    // AB Dot AC = |AB||AC|cosθ
    //θ = |AB|*|AC| / AB Dot AC
    glm::vec3 ab = getVert(1)->get() - getVert(0)->get();
    glm::vec3 ac = getVert(2)->get() - getVert(0)->get();
    float angle = (glm::length(ab)*glm::length(ac))/glm::dot(ab, ac);
    // A = 1/2 |AB||AC|sinθ
    return .5f * glm::length(ab)*glm::length(ac)*sin(angle);
}


/*glm::vec3 ComputeNormal(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
    glm::vec3 v12 = p2;
    v12 -= p1;
    glm::vec3 v23 = p3;
    v23 -= p2;
    glm::vec3 normal = glm::cross(v12,v23);
    //std::cout << "\t\t\tFaceNormal: " << normal.x << ", " << normal.y << ", " << normal.z << "\n";
    normal = glm::normalize(normal);
    //std::cout << "\t\t\tFaceNormal: " << normal.x << ", " << normal.y << ", " << normal.z << "\n";
    return normal;
}*/

int Triangle::next_triangle_id = 0;