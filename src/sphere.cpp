#include "utils.h"
#include "material.h"
#include "argparser.h"
#include "sphere.h"
#include "vertex.h"
#include "mesh.h"
#include "ray.h"
#include "hit.h"
#include "math.h"

// ====================================================================
// ====================================================================

bool Sphere::intersect(const Ray &r, Hit &h) const {




  // ==========================================
  // ASSIGNMENT:  IMPLEMENT SPHERE INTERSECTION
  // ==========================================

  // plug the explicit ray equation into the implict sphere equation and solve
    
    //Ray Equation
    //P(t) = origin + t*direction
    //Sphere Equation
    //H(P) = P*P - r^2 = 0 (P is only if at origin, otherwise it need to be (x+i)^2 + (y+j)^2 + (z+k)^2)
    
    //Combined Equation
    //((Ro + t*Rd) - SphereOrigin)^2 - r^2
    //((X  + T*Y)  -      Z)^2       - R^2
    //T^2*Y^2 + 2*T*X*Y - 2*T*Y*Z + X^2 - 2*X*Z + Z^2 - R^2
    //(Y^2)(T^2) + (2*x*y - 2*Y*Z)(T) + (X^2 - 2*X*Z + Z^2 - R^2)
    //  A                 B                         C
    
    //A = y^2
    //A = Rd^2
    float a = 1;
    
    //B = 2y*(x - z)
    //B = 2Rd*(Ro - SphereOrigin)
    float b = (float)2*glm::dot(r.getDirection(),(r.getOrigin() - center) );
    
    //C = x^2 - 2*x*z + z^2 - r^2
    //C = Ro^2 - 2*Ro*SphereOrigin + SphereOrigin^2 - r^2
    float c = glm::dot(r.getOrigin(),r.getOrigin()) - 2*glm::dot(r.getOrigin(), center) + glm::dot(center, center) - (float)radius*(float)radius;
    
    //float epsilon = .001;
    
    if (b*b - 4*a*c < 0 || (-b + sqrt(b*b - 4*a*c) < 0 && -b - sqrt(b*b - 4*a*c) < 0)) // if the ray doesn't intersect the sphere
    {
        return false;
    }
    
    
    
    //printf("Origin: (%f, %f, %f), %p\n", r.getOrigin().x, r.getOrigin().y, r.getOrigin().z, material);
    
    if (-b + sqrt(b*b - 4*a*c) > 0)   //if outside sphere
    {
        float sol1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
        float sol2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
        
        glm::vec3 normal;
        if(sol1 < sol2) // sol 1 is closer use sol 1
        {
            if(sol1 < 0)              //if smaller one is behind me dont use it
            {
                //if hit.getT is smaller than root or if it's zero
                if (h.getT() == 0 || h.getT() > sol2)
                {
                    //printf("1\n");
                    //get normal using               point of intersection    -    sphere center
                    normal = glm::normalize((r.getOrigin() + sol2*r.getDirection()) - center);
                    
                    //update hit data to be sol2
                    h.set(sol2, material, normal);
                    return true;
                }
            }
            else                            //if it's outside the sphere and fine
            {
                if (h.getT() == 0 || h.getT() > sol1)
                {
                    //printf("2\n");
                    //get normal using               point of intersection    -    sphere center
                    normal = glm::normalize((r.getOrigin() + sol1*r.getDirection()) - center);
                    
                    //update hit data to be sol1
                    h.set(sol1, material, normal);
                    return true;
                }
            }
        }
        else //sol 2 is closer than sol 1
        {
            if(sol2 < 0) //use sol1 instead
            {
                if (h.getT() == 0 || h.getT() > sol1)
                {
                    //printf("3\n");
                    //get normal using               point of intersection    -    sphere center
                    normal = glm::normalize((r.getOrigin() + sol1*r.getDirection()) - center);
                    
                    //update hit data to be sol1
                    h.set(sol1, material, normal);
                    return true;
                }
            }
            else //use sol2
            {
                //printf("4\n");
                //get normal using               point of intersection    -    sphere center
                normal = glm::normalize((r.getOrigin() + sol2*r.getDirection()) - center);
                
                //update hit data to be sol2
                h.set(sol2, material, normal);
                return true;
            }
        }
        
    }
    
  // return true if the sphere was intersected, and update the hit
  // data structure to contain the value of t for the ray at the
  // intersection point, the material, and the normal


  return false;

} 

// ====================================================================
// ====================================================================

// helper function to place a grid of points on the sphere
glm::vec3 ComputeSpherePoint(float s, float t, const glm::vec3 center, float radius) {
  float angle = 2*M_PI*s;
  float y = -cos(M_PI*t);
  float factor = sqrt(1-y*y);
  float x = factor*cos(angle);
  float z = factor*-sin(angle);
  glm::vec3 answer = glm::vec3(x,y,z);
  answer *= radius;
  answer += center;
  return answer;
}

void Sphere::addRasterizedFaces(Mesh *m, ArgParser *args) {
  
  // and convert it into quad patches for radiosity
  int h = args->sphere_horiz;
  int v = args->sphere_vert;
  assert (h % 2 == 0);
  int i,j;
  int va,vb,vc,vd;
  Vertex *a,*b,*c,*d;
  int offset = m->numVertices(); //vertices.size();

  // place vertices
  m->addVertex(center+radius*glm::vec3(0,-1,0));  // bottom
  for (j = 1; j < v; j++) {  // middle
    for (i = 0; i < h; i++) {
      float s = i / float(h);
      float t = j / float(v);
      m->addVertex(ComputeSpherePoint(s,t,center,radius));
    }
  }
  m->addVertex(center+radius*glm::vec3(0,1,0));  // top

  // the middle patches
  for (j = 1; j < v-1; j++) {
    for (i = 0; i < h; i++) {
      va = 1 +  i      + h*(j-1);
      vb = 1 + (i+1)%h + h*(j-1);
      vc = 1 +  i      + h*(j);
      vd = 1 + (i+1)%h + h*(j);
      a = m->getVertex(offset + va);
      b = m->getVertex(offset + vb);
      c = m->getVertex(offset + vc);
      d = m->getVertex(offset + vd);
      m->addRasterizedPrimitiveFace(a,b,d,c,material);
    }
  }

  for (i = 0; i < h; i+=2) {
    // the bottom patches
    va = 0;
    vb = 1 +  i;
    vc = 1 + (i+1)%h;
    vd = 1 + (i+2)%h;
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(d,c,b,a,material);
    // the top patches
    va = 1 + h*(v-1);
    vb = 1 +  i      + h*(v-2);
    vc = 1 + (i+1)%h + h*(v-2);
    vd = 1 + (i+2)%h + h*(v-2);
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(b,c,d,a,material);
  }
}
