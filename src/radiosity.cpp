#include "glCanvas.h"
#include "radiosity.h"
#include "mesh.h"
#include "face.h"
#include "sphere.h"
#include "raytree.h"
#include "raytracer.h"
#include "utils.h"
#include "math.h"

// ================================================================
// CONSTRUCTOR & DESTRUCTOR
// ================================================================
Radiosity::Radiosity(Mesh *m, ArgParser *a) {
  mesh = m;
  args = a;
  num_faces = -1;
    totalFormFactori = NULL;
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
  Reset();
}

Radiosity::~Radiosity() {
  Cleanup();
  cleanupVBOs();
}

void Radiosity::Cleanup() {
  delete [] formfactors;
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;
    delete [] totalFormFactori;
  num_faces = -1;
    totalFormFactori = NULL;
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
}

void Radiosity::Reset() {
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;
    
  // create and fill the data structures
  num_faces = mesh->numFaces();
  area = new float[num_faces];
  undistributed = new glm::vec3[num_faces];
  absorbed = new glm::vec3[num_faces];
  radiance = new glm::vec3[num_faces];
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    f->setRadiosityPatchIndex(i);
    setArea(i,f->getArea());
    glm::vec3 emit = f->getMaterial()->getEmittedColor();
    setUndistributed(i,emit);
    setAbsorbed(i,glm::vec3(0,0,0));
    setRadiance(i,emit);
  }

  // find the patch with the most undistributed energy
  findMaxUndistributed();
}


// =======================================================================================
// =======================================================================================

void Radiosity::findMaxUndistributed() {
  // find the patch with the most undistributed energy 
  // don't forget that the patches may have different sizes!
  max_undistributed_patch = -1;
  total_undistributed = 0;
  total_area = 0;
  float max = -1;
  for (int i = 0; i < num_faces; i++) {
    float m = glm::length(getUndistributed(i)) * getArea(i);
    total_undistributed += m;
    total_area += getArea(i);
    if (max < m) {
      max = m;
      max_undistributed_patch = i;
    }
  }
  assert (max_undistributed_patch >= 0 && max_undistributed_patch < num_faces);
}

//quick helper function to determine if two points are roughly the same
bool comparePoints(glm::vec3 point1, glm::vec3 point2)
{
    if(point1.x < point2.x + .001 && point1.x > point2.x - .001)
        if(point1.y < point2.y + .001 && point1.y > point2.y - .001)
            if(point1.z < point2.z + .001 && point1.z > point2.z - .001)
                return true;
    return false;
}


void Radiosity::ComputeFormFactors() {
  assert (formfactors == NULL);
  assert (num_faces > 0);
  formfactors = new float[num_faces*num_faces];
    totalFormFactori = new float[num_faces];


  // =====================================
  // ASSIGNMENT:  COMPUTE THE FORM FACTORS
  // =====================================

    //use getRasterizedPrimitiveFace(int i) for the forms (maybe, apparently there are 0 of these)
    //use the original geometry for checking occlusion (possibly Face* getOriginalQuad(int i))
    /*printf("NumFormFactors: %d\tNumOriginalQuads: %d\tNumRasterizedPrimitiveFaces: %d\tNumPrimitiveFaces: %d\n", num_faces, mesh->numOriginalQuads(), mesh->numRasterizedPrimitiveFaces(), mesh->numPrimitives()); //14, 14, 0, 0  to  //56, 14, 0, 0
    
    Ray ray1(glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
    Ray ray2(glm::vec3(0, 0, 0), glm::vec3(1, 0, 0));
    
    printf("Angle: %f\n", ray1.getAngle(ray2));*/
    
    
    /*START SMALL (Like .5 or smaller and get larger) double until makes too large, then halve the distance until going too small, take the number before going too small, toss out the extras by directly going through vector and eliminating the n smallest*/
    
    for (int i = 0; i < num_faces; i++)
    {
        totalFormFactori[i] = 0;
        glm::vec3 point1 = mesh->getFace(i)->RandomPoint();
        glm::vec3 faceNormali = mesh->getFace(i)->computeNormal();
        
        //move the start point slightly off the surface to avoid self collision
        point1 += (float).001*glm::normalize(faceNormali);
        
        //float totalFormFactori = 0;
        for (int j = 0; j < num_faces; j++)
        {
            if (i == j)
            {
                setFormFactor(i, j, 0);
                continue;
            }
            glm::vec3 faceNormalj = mesh->getFace(j)->computeNormal();
            float formFactorTemp = 0;
            
            for (int k = 0; k<args->num_form_factor_samples; k++)
            {
                
                glm::vec3 point2 = mesh->getFace(j)->RandomPoint();
                
                //move the end point slightly off the surface to avoid self collision
                //point2 += (float).001*glm::normalize(faceNormalj);
                
                Ray r(point1, glm::normalize(point2-point1));
                
                Hit h;
                raytracer->CastRay(r,h,true);
                glm::vec3 hitPoint = r.pointAtParameter(h.getT());
                if (comparePoints(point2, hitPoint) == 0) //if the points are not roughly the same then we have occlusion
                {
                    continue;
                }
                
                                  //   cos(       ThetaI            )  *cos(        ThetaJ                 ) /
                formFactorTemp +=   ( (cos( r.getAngle( faceNormali ) )*cos( r.getAngle( -1.0f*faceNormalj ) ) ) /
                                  //  PI*R^2 /
                                    ( M_PI*glm::dot( point2-point1, point2-point1 ) ) )/
                                  // numFormFactors
                                    args->num_form_factor_samples;
                
                //printf("(cos(%f)*cos(%f)/3.1415*%f)/%d = %f\n", r.getAngle( faceNormali ), r.getAngle( -1.0f*faceNormalj ), glm::dot( point2-point1, point2-point1 ), args->num_form_factor_samples, ( (cos( r.getAngle( faceNormali ) )*cos( r.getAngle( -1.0f*faceNormalj ) ) ) / ( M_PI*glm::dot( point2-point1, point2-point1 ) ) )/args->num_form_factor_samples);
                //printf("Total Form Factor %d: %f\n", i, totalFormFactori[i]);
                //printf("JPoint: (%f %f %f)\tHitPoint: (%f %f %f)\tSame: %d\n", point2.x, point2.y, point2.z, hitPoint.x, hitPoint.y, hitPoint.z, comparePoints(point2, hitPoint));
                
            }
            
            //printf("FormFactor(%d, %d): %f\n", i, j, formFactorTemp);
            setFormFactor(i, j, formFactorTemp);
            totalFormFactori[i]+= formFactorTemp;
            
            //cast n rays from patch 1 to patch 2 and see how many are occluded
            //integrate point to point form factor Cos(thetaI)*Cos(thetaJ)/PI*|R|^2
            //this will use getFace[i].getRandomPoint()
            
            
            //for testing purposes we are gonna use the centers of the faces to start
            
        }
        //printf("Total Form Factor %d: %f\n", i, totalFormFactori[i]);
    }
    
    //SUM ALL UNDISTRIBUTED EVERY ITERATION***********************************************************************
    
    


}





// ================================================================
// ================================================================

float Radiosity::Iterate() {
  if (formfactors == NULL) 
    ComputeFormFactors();
  assert (formfactors != NULL);




  // ==========================================
  // ASSIGNMENT:  IMPLEMENT RADIOSITY ALGORITHM
  // ==========================================

    //printf("\n\n");
    
    //take the thing with the most undistributed radiosity
    findMaxUndistributed();
    
    float totalUnd = 0;
    
    for(int j = 0; j<num_faces; j++)
    {
        if (j == max_undistributed_patch)
        {
            continue;
        }
        //use max_undistributed_patch for the thing giving away radience
        
        //the j panel absorbs some of the light and reflects the rest
        
        /*printf("MaxUndist: (%f %f %f)\tJDiffuseColor: (%f %f %f)\tAbsorbed: (%f %f %f)\tRadiance: (%f %f %f)\tUndistributed: (%f %f %f)\t\n\n", getUndistributed(max_undistributed_patch).x, getUndistributed(max_undistributed_patch).y, getUndistributed(max_undistributed_patch).z,
               mesh->getFace(j)->getMaterial()->getDiffuseColor().x, mesh->getFace(j)->getMaterial()->getDiffuseColor().y, mesh->getFace(j)->getMaterial()->getDiffuseColor().z,
               getAbsorbed(j).x, getAbsorbed(j).y, getAbsorbed(j).z,
               getRadiance(j).x, getRadiance(j).y, getRadiance(j).x,
               getUndistributed(j).x, getUndistributed(j).y, getUndistributed(j).z);//*/
        
        
        //give it to everyone else based on form factor
        setRadiance(j, getUndistributed(max_undistributed_patch)*(getFormFactor(max_undistributed_patch, j)/
                       totalFormFactori[max_undistributed_patch])*mesh->getFace(j)->getMaterial()->getDiffuseColor() +//- getAbsorbed(j) + *******INSTEAD OF GET ABSORBED THIS SHOULD BE PANEL COLOR)
                        getRadiance(j));//*/
        
        setAbsorbed(j,  (getUndistributed(max_undistributed_patch)*(1.0f-mesh->getFace(j)->getMaterial()->getDiffuseColor()))*
                    (getFormFactor(max_undistributed_patch, j)/totalFormFactori[max_undistributed_patch]) +
                    (getAbsorbed(j)));//*/
        
        setUndistributed(j, getUndistributed(max_undistributed_patch)*(getFormFactor(max_undistributed_patch, j)/
                            totalFormFactori[max_undistributed_patch])*mesh->getFace(j)->getMaterial()->getDiffuseColor() +
                         getUndistributed(j));//*/
        
        totalUnd += glm::length(getUndistributed(j));
        
        /*printf("MaxUndist: (%f %f %f)\tJDiffuseColor: (%f %f %f)\tAbsorbed: (%f %f %f)\tRadiance: (%f %f %f)\tUndistributed: (%f %f %f)\t\n", getUndistributed(max_undistributed_patch).x, getUndistributed(max_undistributed_patch).y, getUndistributed(max_undistributed_patch).z,
               mesh->getFace(j)->getMaterial()->getDiffuseColor().x, mesh->getFace(j)->getMaterial()->getDiffuseColor().y, mesh->getFace(j)->getMaterial()->getDiffuseColor().z,
               getAbsorbed(j).x, getAbsorbed(j).y, getAbsorbed(j).z,
               getRadiance(j).x, getRadiance(j).y, getRadiance(j).x,
               getUndistributed(j).x, getUndistributed(j).y, getUndistributed(j).z);
        
        printf("***************************************************************************************\n");//*/
        
        
        
    }
    //setRadiance(max_undistributed_patch, glm::vec3(0,0,0));
    setUndistributed(max_undistributed_patch, glm::vec3(0, 0, 0));
    
    
    
  //==============================================
  // return the total light yet undistributed
  //==============================================
    
  // (so we can decide when the solution has sufficiently converged)

    //printf("TotalUndistributed: %f", totalUnd);
    
  return totalUnd;




}


// =======================================================================================
// VBO & DISPLAY FUNCTIONS
// =======================================================================================

// for interpolation
void CollectFacesWithVertex(Vertex *have, Face *f, std::vector<Face*> &faces) {
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == f) return;
  }
  if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2] && have != (*f)[3]) return;
  faces.push_back(f);
  for (int i = 0; i < 4; i++) {
    Edge *ea = f->getEdge()->getOpposite();
    Edge *eb = f->getEdge()->getNext()->getOpposite();
    Edge *ec = f->getEdge()->getNext()->getNext()->getOpposite();
    Edge *ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
    if (ea != NULL) CollectFacesWithVertex(have,ea->getFace(),faces);
    if (eb != NULL) CollectFacesWithVertex(have,eb->getFace(),faces);
    if (ec != NULL) CollectFacesWithVertex(have,ec->getFace(),faces);
    if (ed != NULL) CollectFacesWithVertex(have,ed->getFace(),faces);
  }
}

// different visualization modes
glm::vec3 Radiosity::setupHelperForColor(Face *f, int i, int j) {
  assert (mesh->getFace(i) == f);
  assert (j >= 0 && j < 4);
  if (args->render_mode == RENDER_MATERIALS) {
    return f->getMaterial()->getDiffuseColor();
  } else if (args->render_mode == RENDER_RADIANCE && args->interpolate == true) {
    std::vector<Face*> faces;
    CollectFacesWithVertex((*f)[j],f,faces);
    float total = 0;
    glm::vec3 color = glm::vec3(0,0,0);
    glm::vec3 normal = f->computeNormal();
    for (unsigned int i = 0; i < faces.size(); i++) {
      glm::vec3 normal2 = faces[i]->computeNormal();
      float area = faces[i]->getArea();
      if (glm::dot(normal,normal2) < 0.5) continue;
      assert (area > 0);
      total += area;
      color += float(area) * getRadiance(faces[i]->getRadiosityPatchIndex());
    }
    assert (total > 0);
    color /= total;
    return color;
  } else if (args->render_mode == RENDER_LIGHTS) {
    return f->getMaterial()->getEmittedColor();
  } else if (args->render_mode == RENDER_UNDISTRIBUTED) { 
    return getUndistributed(i);
  } else if (args->render_mode == RENDER_ABSORBED) {
    return getAbsorbed(i);
  } else if (args->render_mode == RENDER_RADIANCE) {
    return getRadiance(i);
  } else if (args->render_mode == RENDER_FORM_FACTORS) {
    if (formfactors == NULL) ComputeFormFactors();
    float scale = 0.2 * total_area/getArea(i);
    float factor = scale * getFormFactor(max_undistributed_patch,i);
    return glm::vec3(factor,factor,factor);
  } else {
    assert(0);
  }
  exit(0);
}



void Radiosity::initializeVBOs() {
  // create a pointer for the vertex & index VBOs
  glGenBuffers(1, &mesh_tri_verts_VBO);
  glGenBuffers(1, &mesh_tri_indices_VBO);
  glGenBuffers(1, &mesh_textured_tri_indices_VBO);
}


void Radiosity::setupVBOs() {
  HandleGLError("enter radiosity setupVBOs()");
  mesh_tri_verts.clear();
  mesh_tri_indices.clear();
  mesh_textured_tri_indices.clear();

  // initialize the data in each vector
  int num_faces = mesh->numFaces();
  assert (num_faces > 0);
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    const Edge *e = f->getEdge();
    glm::vec3 normal = f->computeNormal();

    double avg_s = 0;
    double avg_t = 0;
    glm::vec3 avg_color(0,0,0);

    int start = mesh_tri_verts.size();

    // wireframe is normally black, except when it's the special
    // patch, then the wireframe is red
    glm::vec4 wireframe_color(0,0,0,0.5);
    if (args->render_mode == RENDER_FORM_FACTORS && i == max_undistributed_patch) {
      wireframe_color = glm::vec4(1,0,0,1);
    }

    // add the 4 corner vertices
    for (int j = 0; j < 4; j++) {
      glm::vec3 pos = ((*f)[j])->get();
      double s = (*f)[j]->get_s();
      double t = (*f)[j]->get_t();
      glm::vec3 color = setupHelperForColor(f,i,j);
      color = glm::vec3(linear_to_srgb(color.r),
                        linear_to_srgb(color.g),
                        linear_to_srgb(color.b));
      avg_color += 0.25f * color;
      mesh_tri_verts.push_back(VBOPosNormalColor(pos,normal,
                                                 glm::vec4(color.r,color.g,color.b,1.0),
                                                 wireframe_color,
                                                 s,t));
      avg_s += 0.25 * s;
      avg_t += 0.25 * t;
      e = e->getNext();
    }

    // the centroid (for wireframe rendering)
    glm::vec3 centroid = f->computeCentroid();
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,
                                               glm::vec4(avg_color.r,avg_color.g,avg_color.b,1),
                                               glm::vec4(1,1,1,1),
                                               avg_s,avg_t));

    if (f->getMaterial()->hasTextureMap()) {
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+0,start+1,start+4));
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+4));
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+2,start+3,start+4));
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+3,start+0,start+4));
    } else {
      mesh_tri_indices.push_back(VBOIndexedTri(start+0,start+1,start+4));
      mesh_tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+4));
      mesh_tri_indices.push_back(VBOIndexedTri(start+2,start+3,start+4));
      mesh_tri_indices.push_back(VBOIndexedTri(start+3,start+0,start+4));
    }
  }
  assert ((int)mesh_tri_verts.size() == num_faces*5);
  assert ((int)mesh_tri_indices.size() + (int)mesh_textured_tri_indices.size() == num_faces*4);
  
  // copy the data to each VBO
  glBindBuffer(GL_ARRAY_BUFFER,mesh_tri_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,
	       sizeof(VBOPosNormalColor) * num_faces * 5,
	       &mesh_tri_verts[0],
	       GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	       sizeof(VBOIndexedTri) * mesh_tri_indices.size(),
	       &mesh_tri_indices[0], GL_STATIC_DRAW);
  if (mesh_textured_tri_indices.size() > 0) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_textured_tri_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(VBOIndexedTri) * mesh_textured_tri_indices.size(),
                 &mesh_textured_tri_indices[0], GL_STATIC_DRAW);
    
  }

  HandleGLError("radiosity setupVBOs() just before texture");

  // WARNING: this naive VBO implementation only allows a single texture
  // FIXME: something still buggy about textures
  int num_textured_materials = 0;
  for (unsigned int mat = 0; mat < mesh->materials.size(); mat++) {
    Material *m = mesh->materials[mat];
    if (m->hasTextureMap()) {
      // FIXME: old gl...
      //glBindTexture(GL_TEXTURE_2D,m->getTextureID());
      num_textured_materials++;
    }
  }
  assert (num_textured_materials <= 1);

  HandleGLError("leave radiosity setupVBOs()");
}


void Radiosity::drawVBOs() {

  // =====================
  // DRAW ALL THE POLYGONS

    //printf("#MeshTri: %d, #TexMeshTri: %d, numFaces: %d (*4)\n", (int)mesh_tri_indices.size(), (int)mesh_textured_tri_indices.size(), num_faces);
    
  //assert ((int)mesh_tri_indices.size() + (int)mesh_textured_tri_indices.size() == num_faces*4); //**************************ISSUE HERE************************

  // render with Phong lighting?
  if (args->render_mode == RENDER_MATERIALS) {
    // yes
    glUniform1i(GLCanvas::colormodeID, 1);
  } else {
    // no
    glUniform1i(GLCanvas::colormodeID, 0);
  }

  // render untextured faces
  glBindBuffer(GL_ARRAY_BUFFER,mesh_tri_verts_VBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_tri_indices_VBO); 
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
  glEnableVertexAttribArray(4);
  glVertexAttribPointer(4, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)*2));
  glDrawElements(GL_TRIANGLES, mesh_tri_indices.size()*3,GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(3);
  glDisableVertexAttribArray(4);


  // render faces with textures
  if (mesh_textured_tri_indices.size() > 0) {

    // FIXME: there is something buggy with textures still
    //glUniform1i(GLCanvas::colormodeID, 2);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, GLCanvas::textureID);
    GLCanvas::mytexture = glGetUniformLocation(GLCanvas::programID, "mytexture");
    glUniform1i(GLCanvas::mytexture, /*GL_TEXTURE*/0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_textured_tri_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)*2));
    glDrawElements(GL_TRIANGLES, mesh_textured_tri_indices.size()*3, GL_UNSIGNED_INT, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);

    //glUniform1i(GLCanvas::colormodeID, 1);
  }

  HandleGLError(); 
}


void Radiosity::cleanupVBOs() {
  glDeleteBuffers(1, &mesh_tri_verts_VBO);
  glDeleteBuffers(1, &mesh_tri_indices_VBO);
  glDeleteBuffers(1, &mesh_textured_tri_indices_VBO);

  glDeleteTextures(1, &GLCanvas::textureID);
}

