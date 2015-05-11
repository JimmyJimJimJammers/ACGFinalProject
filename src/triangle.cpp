//
//  triangle.cpp
//  
//
//  Created by James McCarthy on 5/3/15.
//
//

#include "triangle.h"
#include "math.h"



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

float Triangle::getCenterDisp()
{
    //printf("In compare disp\n");
    //the current distance from the original surface to the
    float currDisp = 0;
    float tempS = 0;
    float tempT = 0;
    
    assert(material != NULL);
    for (int i = 0; i < 3; i++)
    {
        assert((*this)[i] != NULL);
        //printf("Looking at vertex %d\n", (i+1));
        
        
        //printf("S: %f \n", (*this)[i]->get_s());
        //printf("T: %f \n", (*this)[i]->get_t());
        
        //calculate the S and T for the center point
        tempS += (1.0f/3.0f)*((*this)[i]->get_s());
        //printf("Got S\n");
        tempT += (1.0f/3.0f)*((*this)[i]->get_t());
        //printf("Got T\n");
        //printf("Disp for vert %d is %f\n", i, material->getDisplacementValue((*this)[i]->get_s(), (*this)[i]->get_t()).x);
        
        //calculate the interpolated disp (this is the value that the center point will be sitting at right now)
        currDisp += (1.0f/3.0f)*material->getDisplacementValue((*this)[i]->get_s(), (*this)[i]->get_t()).x;
    }
    //printf("Got the disp values and the S and T values\n");
    //the actual displacement value at the new point
    float realDisp = material->getDisplacementValue(tempS, tempT).x; //compare the difference between the interpolated displacement and the actual dispacement
    
    //printf("Got the real disp value, time to return\n");
    
    //the actual value we want to base out choice off of is how much it will change from where it is now to where it would be
    return fabs(realDisp - currDisp);
}

float Triangle::getAvgDisp()
{
    //printf("In compare disp\n");
    //the current distance from the original surface to the
    float currDisp = 0;
    float tempS = 0;
    float tempT = 0;
    
    assert(material != NULL);
    for (int i = 0; i < 3; i++)
    {
        assert((*this)[i] != NULL);
        //printf("Looking at vertex %d\n", (i+1));
        
        
        //printf("S: %f \n", (*this)[i]->get_s());
        //printf("T: %f \n", (*this)[i]->get_t());
        
        //calculate the S and T for the center point
        tempS += (1.0f/3.0f)*((*this)[i]->get_s());
        //printf("Got S\n");
        tempT += (1.0f/3.0f)*((*this)[i]->get_t());
        //printf("Got T\n");
        //printf("Disp for vert %d is %f\n", i, material->getDisplacementValue((*this)[i]->get_s(), (*this)[i]->get_t()).x);
        
        //calculate the interpolated disp (this is the value that the center point will be sitting at right now)
        currDisp += (1.0f/3.0f)*material->getDisplacementValue((*this)[i]->get_s(), (*this)[i]->get_t()).x;
    }
    //printf("Got the disp values and the S and T values\n");
    //the actual displacement value at the new point
    float realDisp = material->getDisplacementValue(tempS, tempT).x; //compare the difference between the interpolated displacement and the actual dispacement
    
    //lets also factor in the 3 points between the center point and the outer vertices
    float caDisp = .5f*material->getDisplacementValue((*this)[0]->get_s(), (*this)[0]->get_t()).x;
    caDisp += .5f*currDisp;
    float caReal = material->getDisplacementValue((tempS+(*this)[0]->get_s())/2.0f, (tempT+(*this)[0]->get_t())/2.0f).x;
    
    float cbDisp = .5f*material->getDisplacementValue((*this)[1]->get_s(), (*this)[1]->get_t()).x;
    cbDisp += .5f*currDisp;
    float cbReal = material->getDisplacementValue((tempS+(*this)[1]->get_s())/2.0f, (tempT+(*this)[1]->get_t())/2.0f).x;
    
    float ccDisp = .5f*material->getDisplacementValue((*this)[2]->get_s(), (*this)[2]->get_t()).x;
    ccDisp += .5f*currDisp;
    float ccReal = material->getDisplacementValue((tempS+(*this)[2]->get_s())/2.0f, (tempT+(*this)[2]->get_t())/2.0f).x;
    
    
    
    //printf("Got the real disp value, time to return\n");
    
    //the actual value we want to base out choice off of is how much it will change from where it is now to where it would be
    return .5f*fabs(realDisp - currDisp) + (1.0f/6.0f)*fabs(caReal - caDisp) + (1.0f/6.0f)*fabs(cbReal - cbDisp) + (1.0f/6.0f)*fabs(ccReal - ccDisp);
}

glm::vec3 Triangle::computeNormal() const
{
    // note: this face might be non-planar, so average the two triangle normals
    //printf("Problem?1\n");
    assert(this->getEdge() != NULL);
    //printf("Problem?1.25\n");
    assert(this->getEdge()->getNext() != NULL);
    //printf("Problem?1.5\n");
    assert(this->getEdge()->getNext()->getNext() != NULL);
    //printf("Problem?1.75\n");
    
    glm::vec3 a = (*this)[0]->get();
    //printf("Problem?2\n");
    glm::vec3 b = (*this)[1]->get();
    //printf("Problem?3\n");
    glm::vec3 c = (*this)[2]->get();
    //printf("Problem?4\n");
    //printf("Set the position vectors to compute normal\n");
    glm::vec3 v12 = b - a;
    glm::vec3 v23 = c - b;
    //printf("Did some math\n");
    //printf("Problem?5\n");
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

/*void Triangle::setVertNull(int i)
{
    if (i==0) edge->setStartVertex(NULL);
    if (i==1) edge->getNext()->setStartVertex(NULL);
    if (i==2) edge->getNext()->getNext()->setStartVertex(NULL);
}*/


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