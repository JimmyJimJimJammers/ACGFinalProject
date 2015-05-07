#include "glCanvas.h"

#include <iostream>
#include <algorithm>

#include "argparser.h"
#include "photon_mapping.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "kdtree.h"
#include "utils.h"
#include "raytracer.h"
#include "math.h"
#include "unistd.h"
#include <climits>
#include "math.h"


// ==========
// DESTRUCTOR
PhotonMapping::~PhotonMapping() {
  // cleanup all the photons
  delete kdtree;
}

glm::vec3 generateRandomDirection(glm::vec3 normal, ArgParser* args)
{
    //generate random point in a 1x1x1 cube
    //Then throw away values from outside the unit sphere
    //also throw out values that create an angle greater than 90 degrees from the normal
    glm::vec3 retDir(2, 2, 2);
    ////printf("\nRadius > 1: %f\tAngle: %f\n", retDir.x*retDir.x + retDir.y*retDir.y + retDir.z*retDir.z, glm::dot(normal, retDir));
    while(retDir.x*retDir.x + retDir.y*retDir.y + retDir.z*retDir.z > 1 || acos(glm::dot(normal, retDir)/(glm::length(retDir)*glm::length(normal))) > M_PI/2.0f)
    {
        retDir = glm::vec3(args->getRandf()*2 - 1, args->getRandf()*2 - 1, args->getRandf()*2 - 1);
        //retDir = glm::normalize(retDir);
        ////printf("RetDir: (%f %f %f)\tRadius: %f\tAngle: %f\n", retDir.x, retDir.y, retDir.z, retDir.x*retDir.x + retDir.y*retDir.y + retDir.z*retDir.z, glm::dot(normal, retDir));
    }
    
    return glm::normalize(retDir);
}

// ========================================================================
// Recursively trace a single photon

void PhotonMapping::TracePhoton(const glm::vec3 &position, const glm::vec3 &direction, 
				const glm::vec3 &energy, int iter) {



    ////printf("Started Tracing\n");
  // ==============================================
  // ASSIGNMENT: IMPLEMENT RECURSIVE PHOTON TRACING
  // ==============================================

  // Trace the photon through the scene.  At each diffuse or
  // reflective bounce, store the photon in the kd tree.
    
    
    ////printf("Segfault!\n");
    //store photon
    ////printf("Pos: (%f %f %f)\tDirection: (%f %f %f)\tEnergy: (%f %f %f)\tIteration: %d\n", position.x, position.y, position.z, direction.x, direction.y, direction.z, energy.x, energy.y, energy.z, iter);
    /*if (iter > 0)
    {
        Photon p(position, direction, energy, iter);
        kdtree->AddPhoton(p);
    }//*/
    
    glm::vec3 newPosition;
    glm::vec3 newDirection;
    glm::vec3 newEnergy;
    
    //if (glm::length(energy) > .001)/**************************CHANGE THIS TO A FRACTION OF THE NUMBER OF PHOTONS BEING SHOT************/
    if (iter < args->num_bounces)
    {
        ////printf("\tIN THE IF!\n");
        //to avoid self collisions I'm gonna move the position slightly forward in "direction"
        newPosition = position + .001f*direction;
        
        //cast the ray
        Ray r(newPosition, direction);
        Hit h;
        ////printf("\tBefore Trace Ray!\n");
        raytracer->CastRay(r, h, false);
        ////printf("\tAfter Trace Ray!\n");
        
        //adjust position
        newPosition = r.pointAtParameter(h.getT());
        
        ////printf("\tAfter New Pos Set!\n");
        
        if (h.getMaterial() != NULL)
        {
            //adjust energy and direction according to material
            if (glm::length(h.getMaterial()->getDiffuseColor()) > 0)
            {
                ////printf("\t\tDiffuse\n\n");
                //adjust energy
                newEnergy = energy*h.getMaterial()->getDiffuseColor();
                //adjust direction
                newDirection = RandomDiffuseDirection(h.getNormal());//generateRandomDirection(h.getNormal(), args);
            }
            if (glm::length(h.getMaterial()->getReflectiveColor()) > 0)
            {
                ////printf("\t\tReflective\n\n");
                //adjust energy
                newEnergy = energy*h.getMaterial()->getReflectiveColor();
                //adjust direction
                newDirection = direction - ((float)2*glm::dot(direction, h.getNormal()))*h.getNormal();
            }
            
            ////printf("\tRight before trace photon!\n");
            if (iter > 0)
            {
                Photon p(newPosition, direction, newEnergy, iter+1);
                kdtree->AddPhoton(p);
            }

            TracePhoton(newPosition, newDirection, newEnergy, iter+1);
            
        }
    }

  // One optimization is to *not* store the first bounce, since that
  // direct light can be efficiently computed using classic ray
  // tracing.

    ////printf("Finished Tracing\n");


}


// ========================================================================
// Trace the specified number of photons through the scene

void PhotonMapping::TracePhotons() {
    ////printf("trace photons\n");
    
  // first, throw away any existing photons
  delete kdtree;
    

  // consruct a kdtree to store the photons
  BoundingBox *bb = mesh->getBoundingBox();
  glm::vec3 min = bb->getMin();
  glm::vec3 max = bb->getMax();
  glm::vec3 diff = max-min;
  min -= 0.001f*diff;
  max += 0.001f*diff;
  kdtree = new KDTree(BoundingBox(min,max));

    
  // photons emanate from the light sources
  const std::vector<Face*>& lights = mesh->getLights();

  // compute the total area of the lights
  float total_lights_area = 0;
  for (unsigned int i = 0; i < lights.size(); i++) {
    total_lights_area += lights[i]->getArea();
  }

  // shoot a constant number of photons per unit area of light source
  // (alternatively, this could be based on the total energy of each light)
  for (unsigned int i = 0; i < lights.size(); i++) {
    float my_area = lights[i]->getArea();
    int num = args->num_photons_to_shoot * my_area / total_lights_area;
    // the initial energy for this photon
    glm::vec3 energy = my_area/float(num) * lights[i]->getMaterial()->getEmittedColor();
    glm::vec3 normal = lights[i]->computeNormal();
      
    for (int j = 0; j < num; j++) {
      glm::vec3 start = lights[i]->RandomPoint();
      // the initial direction for this photon (for diffuse light sources)
      glm::vec3 direction = RandomDiffuseDirection(normal);
        
      TracePhoton(start,direction,energy,0);
    }
  }
}


// ======================================================================



struct photonDist
{
    float distToPoint;
    Photon* photon;
};


// helper function
bool closest_photon(const std::pair<Photon,float> &a, const std::pair<Photon,float> &b) {
  return (a.second < b.second);
}


struct less_than_key
{
    inline bool operator() (const photonDist &photonDist1, const photonDist &photonDist2)
    {
        return (photonDist1.distToPoint < photonDist2.distToPoint);
    }
};


//gets rid of the farthest photons leaving only the desired # of photons within the radius
bool reduceNumPhotons(int targetNum, std::vector<Photon> &photons, const glm::vec3 &point,float &radius)
{
    std::vector<photonDist> listForSort;
    
    
    
    //step1 go through and eliminate photons not in radius, if too few are left then return false
    for (int i = 0; i < photons.size(); i++)
    {
        float dist = glm::distance(point, photons[i].getPosition());
        if (dist <= radius)
        {
            photonDist temp;
            temp.distToPoint = dist;
            temp.photon = &(photons[i]);
            listForSort.push_back(temp);
        }
    }
    if (listForSort.size() < targetNum)
    {
        //printf("REJECTED!\n");
        radius *= 2;
        return false;
    }
    
    std::sort(listForSort.begin(), listForSort.end(), less_than_key());
    
    photons.clear();
    
    
    float farthestDist = 0;
    for (int i = 0; i<targetNum; i++)
    {
        if (listForSort[i].distToPoint > farthestDist)
        {
            farthestDist = listForSort[i].distToPoint;
        }
        ////printf("Dist: %f\n", listForSort[i].distToPoint);
        photons.push_back(*listForSort[i].photon);
    }
    
    //change radius to be the distance to the farthest point in our list
    radius = farthestDist;
    
    return true;
}


// ======================================================================
glm::vec3 PhotonMapping::GatherIndirect(const glm::vec3 &point, const glm::vec3 &normal, const glm::vec3 &direction_from) const {


  if (kdtree == NULL) { 
    std::cout << "WARNING: Photons have not been traced throughout the scene." << std::endl;
    return glm::vec3(0,0,0); 
  }



  // ================================================================
  // ASSIGNMENT: GATHER THE INDIRECT ILLUMINATION FROM THE PHOTON MAP
  // ================================================================

  // collect the closest args->num_photons_to_collect photons
  // determine the radius that was necessary to collect that many photons
  // average the energy of those photons over that radius
    float radius = args->avgRadius;
    //float radius = args->num_photons_to_collect/(float)args->num_photons_to_shoot;
    BoundingBox box(point-radius, point+radius);
    std::vector<Photon> photons;
    //int previousNums[4] = {-1, -1, -1, -1};
    //int i = 0;
    //bool lastUp = false;
    //bool lastDown = false;
    glm::vec3 answer(0, 0, 0);
    
    do
    {
        while(photons.size() < args->num_photons_to_collect)
        {
            photons.clear();
            //collect photons
            kdtree->CollectPhotonsInBox(box, photons);
            //printf("\tPhotons: %d\tRadius: %f\tChangeAmount: %f\n", (int)photons.size(), radius, changeAmount);
            
            //if we have too many
            if (photons.size() < args->num_photons_to_collect)
            {
                //printf("FAILURE!\n");
                radius *= 2;
                box = BoundingBox(point-radius, point+radius);
            }
        }
    }//while the number of reduced photons in the radius does not match the target keep expanding the radius
    while (!reduceNumPhotons(args->num_photons_to_collect, photons, point, radius));

    /*float changeAmount = .5;
    float radius = 1;
    BoundingBox box(point-radius, point+radius);
    std::vector<Photon> photons;
    int previousNums[4] = {-1, -1, -1, -1};
    int i = 0;
    bool lastUp = false;
    bool lastDown = false;
    glm::vec3 answer(0, 0, 0);
    
    while(photons.size() > args->num_photons_to_collect + args->num_photons_to_collect || photons.size() < args->num_photons_to_collect - args->num_photons_to_collect*.5f )
    {
        ////printf("Goal: %d\n", args->num_photons_to_collect);
        ////printf("Photons: %d\tRadius: %f\tChangeAmount: %f\n", (int)photons.size(), radius, changeAmount);
        if(changeAmount >= radius)
        {
            changeAmount = .5*radius;
        }
        
        photons.clear();
        //collect photons
        kdtree->CollectPhotonsInBox(box, photons);
        ////printf("\tPhotons: %d\tRadius: %f\tChangeAmount: %f\n", (int)photons.size(), radius, changeAmount);
        
        previousNums[i] = photons.size();
        i++;
        i = i%4;
        
        //if we have a flickering between two numbers then break out.
        if (previousNums[0] == previousNums[2] && previousNums[1] == previousNums[3])
        {
            break;
        }
        
        
        //sleep(2);
        
        //if we have too many
        if (photons.size() > args->num_photons_to_collect)
        {
            if(lastUp)
            {
                changeAmount *= .5;
                lastUp = false;
            }
            radius -= changeAmount;
            box = BoundingBox(point-radius, point+radius);
            lastDown = true;
        }
        else if(photons.size() < args->num_photons_to_collect) //if we have too few
        {
            if(lastDown)
            {
                changeAmount *= .5;
                lastDown = false;
            }
            
            radius += changeAmount;
            box = BoundingBox(point-radius, point+radius);
            lastUp = true;
        }
    }*/
    
    
    
    /*if (args->lastRadius == radius) //need to decrease radius if consistently overshooting
    {
        args->avgRadius -= .25f*args->avgRadius;
    }
    else //otherwise take the midpoint between this radius and the average radius
    {*/
        args->avgRadius += radius;
        args->avgRadius /= 2;
    //}
    args->numRadiuses++;
    args->lastRadius = radius;
    
    
    for (int i = 0; i < photons.size(); i++)
    {
        answer += photons[i].getEnergy();
    }
    
    answer /= M_PI*radius*radius;
    
    answer *= 3;
    
    //printf("Added photon energy: (%f %f %f)\tradius: %f\tavgRadius: %f\n\n", answer.x, answer.y, answer.z, radius, args->avgRadius);

  // return the color
  return answer;


}


// ======================================================================
// PHOTON VISUALIZATION FOR DEBUGGING
// ======================================================================

void PhotonMapping::initializeVBOs() {
  HandleGLError("enter photonmapping initializevbos()");
  glGenBuffers(1, &photon_direction_verts_VBO);
  glGenBuffers(1, &photon_direction_indices_VBO);
  glGenBuffers(1, &kdtree_verts_VBO);
  glGenBuffers(1, &kdtree_edge_indices_VBO);
  HandleGLError("leave photonmapping initializevbos()");
}

void PhotonMapping::setupVBOs() {
  HandleGLError("enter photonmapping setupvbos()");

  photon_direction_verts.clear();
  photon_direction_indices.clear();
  kdtree_verts.clear();
  kdtree_edge_indices.clear();

  // initialize the data
  BoundingBox *bb = mesh->getBoundingBox();
  float max_dim = bb->maxDim();

  if (kdtree == NULL) return;
  std::vector<const KDTree*> todo;  
  todo.push_back(kdtree);
  while (!todo.empty()) {
    const KDTree *node = todo.back();
    todo.pop_back(); 
    if (node->isLeaf()) {

      // initialize photon direction vbo
      const std::vector<Photon> &photons = node->getPhotons();
      int num_photons = photons.size();
      for (int i = 0; i < num_photons; i++) {
	const Photon &p = photons[i];
	glm::vec3 energy = p.getEnergy()*float(args->num_photons_to_shoot);
        glm::vec4 color(energy.x,energy.y,energy.z,1);
	const glm::vec3 &position = p.getPosition();
	glm::vec3 other = position - p.getDirectionFrom()*0.02f*max_dim;
        addEdgeGeometry(photon_direction_verts,photon_direction_indices,
                        position,other,color,color,max_dim*0.0005f,0);
      }

      // initialize kdtree vbo
      float thickness = 0.001*max_dim;
      glm::vec3 A = node->getMin();
      glm::vec3 B = node->getMax();
      glm::vec4 black(1,0,0,1);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,A.z),glm::vec3(A.x,A.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,B.z),glm::vec3(A.x,B.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,B.z),glm::vec3(A.x,B.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,A.z),glm::vec3(A.x,A.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,A.y,A.z),glm::vec3(B.x,A.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,A.y,B.z),glm::vec3(B.x,B.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,B.y,B.z),glm::vec3(B.x,B.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,B.y,A.z),glm::vec3(B.x,A.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,A.z),glm::vec3(B.x,A.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,B.z),glm::vec3(B.x,A.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,B.z),glm::vec3(B.x,B.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,A.z),glm::vec3(B.x,B.y,A.z),black,black,thickness,thickness);

    } else {
      todo.push_back(node->getChild1());
      todo.push_back(node->getChild2());
    } 
  }



  // copy the data to each VBO
  if (photon_direction_verts.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,photon_direction_verts_VBO); 
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(VBOPosNormalColor) * photon_direction_verts.size(),
                 &photon_direction_verts[0],
                 GL_STATIC_DRAW); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,photon_direction_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(VBOIndexedTri) * photon_direction_indices.size(),
                 &photon_direction_indices[0], GL_STATIC_DRAW);
  }
  if (kdtree_verts.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,kdtree_verts_VBO); 
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(VBOPosNormalColor) * kdtree_verts.size(),
                 &kdtree_verts[0],
                 GL_STATIC_DRAW); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,kdtree_edge_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(VBOIndexedTri) * kdtree_edge_indices.size(),
                 &kdtree_edge_indices[0], GL_STATIC_DRAW);
  }

  HandleGLError("leave photonmapping setupvbos()");
}

void PhotonMapping::drawVBOs() {
  HandleGLError("enter photonmapping drawvbos()");

  glUniform1i(GLCanvas::colormodeID, 1);
  if (args->render_photons && photon_direction_verts.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,photon_direction_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,photon_direction_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
    glDrawElements(GL_TRIANGLES,
                   photon_direction_indices.size()*3,
                   GL_UNSIGNED_INT, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
  }

  if (args->render_kdtree && kdtree_edge_indices.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,kdtree_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,kdtree_edge_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
    glDrawElements(GL_TRIANGLES,
                   kdtree_edge_indices.size()*3,
                   GL_UNSIGNED_INT, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
  }

  HandleGLError("leave photonmapping drawvbos()");
}

void PhotonMapping::cleanupVBOs() {
  glDeleteBuffers(1, &photon_direction_verts_VBO);
  glDeleteBuffers(1, &photon_direction_indices_VBO);
  glDeleteBuffers(1, &kdtree_verts_VBO);
  glDeleteBuffers(1, &kdtree_edge_indices_VBO);
}

