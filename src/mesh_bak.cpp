#include "glCanvas.h"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <string>
#include <utility>

#include "argparser.h"
#include "vertex.h"
#include "boundingbox.h"
#include "mesh.h"
#include "edge.h"
#include "face.h"
#include "primitive.h"
#include "sphere.h"
#include "cylinder_ring.h"
#include "ray.h"
#include "hit.h"
#include "camera.h"

#include <cfloat>


// =======================================================================
// DESTRUCTOR
// =======================================================================

Mesh::~Mesh() {
  unsigned int i;
  for (i = 0; i < rasterized_primitive_faces.size(); i++) {
    Face *f = rasterized_primitive_faces[i];
    removeFaceEdges(f);
    delete f;
  }
  if (subdivided_quads.size() != original_quads.size()) {
    for (i = 0; i < subdivided_quads.size(); i++) {
      Face *f = subdivided_quads[i];
      removeFaceEdges(f);
      delete f;
    }
  }
  for (i = 0; i < original_quads.size(); i++) {
    Face *f = original_quads[i];
    removeFaceEdges(f);
    delete f;
  }
  for (i = 0; i < primitives.size(); i++) { delete primitives[i]; }
  for (i = 0; i < materials.size(); i++) { delete materials[i]; }
  for (i = 0; i < vertices.size(); i++) { delete vertices[i]; }
  delete bbox;
}

// =======================================================================
// MODIFIERS:   ADD & REMOVE
// =======================================================================

Vertex* Mesh::addVertex(const glm::vec3 &position) {
  int index = numVertices();
  vertices.push_back(new Vertex(index,position));
  // extend the bounding box to include this point
  if (bbox == NULL) 
    bbox = new BoundingBox(position,position);
  else 
    bbox->Extend(position);
  return vertices[index];
}

void Mesh::addPrimitive(Primitive* p) {
  primitives.push_back(p);
  p->addRasterizedFaces(this,args);
}

void Mesh::addFace(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material, enum FACE_TYPE face_type) {
  // create the face
  Face *f = new Face(material);
  // create the edges
  Edge *ea = new Edge(a,b,f);
  Edge *eb = new Edge(b,c,f);
  Edge *ec = new Edge(c,d,f);
  Edge *ed = new Edge(d,a,f);
  // point the face to one of its edges
  f->setEdge(ea);
  // connect the edges to each other
  ea->setNext(eb);
  eb->setNext(ec);
  ec->setNext(ed);
  ed->setNext(ea);
  // verify these edges aren't already in the mesh 
  // (which would be a bug, or a non-manifold mesh)
  assert (edges.find(std::make_pair(a,b)) == edges.end());
  assert (edges.find(std::make_pair(b,c)) == edges.end());
  assert (edges.find(std::make_pair(c,d)) == edges.end());
  assert (edges.find(std::make_pair(d,a)) == edges.end());
  // add the edges to the master list
  edges[std::make_pair(a,b)] = ea;
  edges[std::make_pair(b,c)] = eb;
  edges[std::make_pair(c,d)] = ec;
  edges[std::make_pair(d,a)] = ed;
  // connect up with opposite edges (if they exist)
  edgeshashtype::iterator ea_op = edges.find(std::make_pair(b,a)); 
  edgeshashtype::iterator eb_op = edges.find(std::make_pair(c,b)); 
  edgeshashtype::iterator ec_op = edges.find(std::make_pair(d,c)); 
  edgeshashtype::iterator ed_op = edges.find(std::make_pair(a,d)); 
  if (ea_op != edges.end()) { ea_op->second->setOpposite(ea); }
  if (eb_op != edges.end()) { eb_op->second->setOpposite(eb); }
  if (ec_op != edges.end()) { ec_op->second->setOpposite(ec); }
  if (ed_op != edges.end()) { ed_op->second->setOpposite(ed); }
  // add the face to the appropriate master list
  if (face_type == FACE_TYPE_ORIGINAL) {
    original_quads.push_back(f);
    subdivided_quads.push_back(f);
  } else if (face_type == FACE_TYPE_RASTERIZED) {
    rasterized_primitive_faces.push_back(f); 
  } else {
    assert (face_type == FACE_TYPE_SUBDIVIDED);
    subdivided_quads.push_back(f);
  }
  // if it's a light, add it to that list too
  if (glm::length(material->getEmittedColor()) > 0 && face_type == FACE_TYPE_ORIGINAL) {
    original_lights.push_back(f);
  }
}

void Mesh::removeFaceEdges(Face *f) {
  // helper function for face deletion
  Edge *ea = f->getEdge();
  Edge *eb = ea->getNext();
  Edge *ec = eb->getNext();
  Edge *ed = ec->getNext();
  assert (ed->getNext() == ea);
  Vertex *a = ea->getStartVertex();
  Vertex *b = eb->getStartVertex();
  Vertex *c = ec->getStartVertex();
  Vertex *d = ed->getStartVertex();
  // remove elements from master lists
  edges.erase(std::make_pair(a,b)); 
  edges.erase(std::make_pair(b,c)); 
  edges.erase(std::make_pair(c,d)); 
  edges.erase(std::make_pair(d,a)); 
  // clean up memory
  delete ea;
  delete eb;
  delete ec;
  delete ed;
}

// ==============================================================================
// EDGE HELPER FUNCTIONS

Edge* Mesh::getEdge(Vertex *a, Vertex *b) const {
  edgeshashtype::const_iterator iter = edges.find(std::make_pair(a,b));
  if (iter == edges.end()) return NULL;
  return iter->second;
}

Vertex* Mesh::getChildVertex(Vertex *p1, Vertex *p2) const {
  vphashtype::const_iterator iter = vertex_parents.find(std::make_pair(p1,p2)); 
  if (iter == vertex_parents.end()) return NULL;
  return iter->second; 
}

void Mesh::setParentsChild(Vertex *p1, Vertex *p2, Vertex *child) {
  assert (vertex_parents.find(std::make_pair(p1,p2)) == vertex_parents.end());
  vertex_parents[std::make_pair(p1,p2)] = child; 
}

//
// ===============================================================================
// the load function parses our (non-standard) extension of very simple .obj files
// ===============================================================================

void Mesh::Load(ArgParser *_args) {
  args = _args;

  std::string file = args->path+'/'+args->input_file;

  std::ifstream objfile(file.c_str());
  if (!objfile.good()) {
    std::cout << "ERROR! CANNOT OPEN " << file << std::endl;
    return;
  }

  std::string token;
  Material *active_material = NULL;
  camera = NULL;
  background_color = glm::vec3(1,1,1);
 
  while (objfile >> token) {
    if (token == "v") {
      float x,y,z;
      objfile >> x >> y >> z;
      addVertex(glm::vec3(x,y,z));
    } else if (token == "vt") {
      assert (numVertices() >= 1);
      float s,t;
      objfile >> s >> t;
      getVertex(numVertices()-1)->setTextureCoordinates(s,t);
    } else if (token == "f") {
      int a,b,c,d;
      objfile >> a >> b >> c >> d;
      a--;
      b--;
      c--;
      d--;
      assert (a >= 0 && a < numVertices());
      assert (b >= 0 && b < numVertices());
      assert (c >= 0 && c < numVertices());
      assert (d >= 0 && d < numVertices());
      assert (active_material != NULL);
      addOriginalQuad(getVertex(a),getVertex(b),getVertex(c),getVertex(d),active_material);
    }
    else if (token == "s") {
      float x,y,z,r;
      objfile >> x >> y >> z >> r;
      assert (active_material != NULL);
      addPrimitive(new Sphere(glm::vec3(x,y,z),r,active_material));
    }
    else if (token == "r") {
      float x,y,z,h,r,r2;
      objfile >> x >> y >> z >> h >> r >> r2;
      assert (active_material != NULL);
      addPrimitive(new CylinderRing(glm::vec3(x,y,z),h,r,r2,active_material));
    }
    else if (token == "background_color") {
      float r,g,b;
      objfile >> r >> g >> b;
      background_color = glm::vec3(r,g,b);
    }
    else if (token == "PerspectiveCamera") {
      camera = new PerspectiveCamera();
      objfile >> *(PerspectiveCamera*)camera;
    }
    else if (token == "OrthographicCamera") {
      camera = new OrthographicCamera();
      objfile >> *(OrthographicCamera*)camera;
    }
    else if (token == "m") {
      // this is not standard .obj format!!
      // materials
      int m;
      objfile >> m;
      assert (m >= 0 && m < (int)materials.size());
      active_material = materials[m];
    }
    else if (token == "material")
    {
      // this is not standard .obj format!!
      std::string texture_file = "";
        std::string displace_file = "";
        
      glm::vec3 diffuse(0,0,0);
        glm::vec3 displace(0, 0, 0);
      float r,g,b;
        float j, k, l;
        
      objfile >> token;
        
        ////printf("\n\nToken: %s\n", token.c_str());
        if (token == "displace_file")
        {
            objfile >> displace_file;
            // prepend the directory name
            displace_file = args->path + '/' + displace_file;
            objfile >> token;
        }
        else if(token == "displace")
        {
            objfile >> j >> k >> l;
            displace = glm::vec3(j, k, l);
            objfile >> token;
        }
        
        ////printf("DisplaceFile: %s\n", displace_file.c_str());
        
        if (token == "diffuse")
        {
            objfile >> r >> g >> b;
            diffuse = glm::vec3(r,g,b);
        }
        else
        {
            assert (token == "texture_file");
            objfile >> texture_file;
            // prepend the directory name
            texture_file = args->path + '/' + texture_file;
        }
        
      glm::vec3 reflective, emitted;
      objfile >> token >> r >> g >> b;
      assert (token == "reflective");
      reflective = glm::vec3(r,g,b);
      float roughness = 0;
      objfile >> token;
      if (token == "roughness") {
	objfile >> roughness;
	objfile >> token;
      } 
      assert (token == "emitted");
      objfile >> r >> g >> b;
      emitted = glm::vec3(r,g,b);
      materials.push_back(new Material(texture_file, displace_file, diffuse, displace, reflective,emitted,roughness));
    }
    else {
      std::cout << "UNKNOWN TOKEN " << token << std::endl;
      exit(0);
    }
  }
  std::cout << " mesh loaded: " << numFaces() << " faces and " << numEdges() << " edges." << std::endl;

    /*for (int i = 0; i<materials.size(); i++)
    {
        //printf("DisplaceFile Mesh: %d\n", materials[i]->hasDisplacementMap());
    }*/
    
    
    
  if (camera == NULL) {
    // if not initialized, position a perspective camera and scale it so it fits in the window
    assert (bbox != NULL);
    glm::vec3 point_of_interest; bbox->getCenter(point_of_interest);
    float max_dim = bbox->maxDim();
    glm::vec3 camera_position = point_of_interest + glm::vec3(0,0,4*max_dim);
    glm::vec3 up = glm::vec3(0,1,0);
    camera = new PerspectiveCamera(camera_position, point_of_interest, up, 20 * M_PI/180.0);    
  }
}

// =================================================================
// SUBDIVISION
// =================================================================

Vertex* Mesh::AddEdgeVertex(Vertex *a, Vertex *b) {
  Vertex *v = getChildVertex(a,b);
  if (v != NULL) return v;
  glm::vec3 pos = 0.5f*a->get() + 0.5f*b->get();
  float s = 0.5f*a->get_s() + 0.5f*b->get_s();
  float t = 0.5f*a->get_t() + 0.5f*b->get_t();
  v = addVertex(pos);
  v->setTextureCoordinates(s,t);
  setParentsChild(a,b,v);
  return v;
}

//gets the vert at the center of the face and set it's S and T to be the same proportion of the 4 original verts
Vertex* Mesh::AddMidVertex(Vertex *a, Vertex *b, Vertex *c, Vertex *d) {
  glm::vec3 pos = 0.25f*a->get() + 0.25f*b->get() + 0.25f*c->get() + 0.25f*d->get();
  float s = 0.25f*a->get_s() + 0.25f*b->get_s() + 0.25f*c->get_s() + 0.25f*d->get_s();
  float t = 0.25f*a->get_t() + 0.25f*b->get_t() + 0.25f*c->get_t() + 0.25f*d->get_t();
  Vertex *v = addVertex(pos);
  v->setTextureCoordinates(s,t);
  return v;
}

void Mesh::Subdivision()
{

  //if it's the firs time it's being subdivided
  bool first_subdivision = false;
  if (original_quads.size() == subdivided_quads.size()) {
    first_subdivision = true;
  }

    //get all of the subdivided faces
  std::vector<Face*> tmp = subdivided_quads;
    //delete the old list of them
  subdivided_quads.clear();
  
    //for all subdivided faces
  for (unsigned int i = 0; i < tmp.size(); i++)
  {
      //make a temp face
    Face *f = tmp[i];
    
      //make some temp verts for that face
    Vertex *a = (*f)[0];
    Vertex *b = (*f)[1];
    Vertex *c = (*f)[2];
    Vertex *d = (*f)[3];
    // add new vertices on the edges
    Vertex *ab = AddEdgeVertex(a,b);
    Vertex *bc = AddEdgeVertex(b,c);
    Vertex *cd = AddEdgeVertex(c,d);
    Vertex *da = AddEdgeVertex(d,a);
    // add new point in the middle of the patch
    Vertex *mid = AddMidVertex(a,b,c,d);

      //make sure all the old edges are not null
    assert (getEdge(a,b) != NULL);
    assert (getEdge(b,c) != NULL);
    assert (getEdge(c,d) != NULL);
    assert (getEdge(d,a) != NULL);

    // copy the color and emission from the old patch to the new
    Material *material = f->getMaterial();
    if (!first_subdivision) {
      removeFaceEdges(f);
      delete f;
    }

    // create the new faces
    addSubdividedQuad(a,ab,mid,da,material);
    addSubdividedQuad(b,bc,mid,ab,material);
    addSubdividedQuad(c,cd,mid,bc,material);
    addSubdividedQuad(d,da,mid,cd,material);

      //make sure all the new edges exist
    assert (getEdge(a,ab) != NULL);
    assert (getEdge(ab,b) != NULL);
    assert (getEdge(b,bc) != NULL);
    assert (getEdge(bc,c) != NULL);
    assert (getEdge(c,cd) != NULL);
    assert (getEdge(cd,d) != NULL);
    assert (getEdge(d,da) != NULL);
    assert (getEdge(da,a) != NULL);
  }
}

// =================================================================
// DISPLACEMENT SUBDIVISION
// =================================================================

Mesh::newVertReturn Mesh::AddEdgeVertexDisp(Vertex *a, Vertex *b, Material *m)
{
    Mesh::newVertReturn returnValue;
    
    
    Vertex *v = getChildVertex(a,b);
    if (v != NULL)
    {
        returnValue.newVert = v;
        returnValue.significant = false;
        
        return returnValue;
    }
    //glm::vec3 pos = 0.5f*a->get() + 0.5f*b->get();
    //float s = 0.5f*a->get_s() + 0.5f*b->get_s();
    //float t = 0.5f*a->get_t() + 0.5f*b->get_t();
    
    glm::vec3 pos(0, 0, 0);
    float s = 0;
    float t = 0;
    float epsilon = 0.00001;
    
    float maxHeight = FLT_MIN;
    
    for (int i = 0; i < args->disp_increments; i++)
    {
        float fraction = (float)(i+1)/(float)(args->disp_increments + 1);
        glm::vec3 tempPos = fraction*a->get() + (1-fraction)*b->get();
        float tempS = fraction*a->get_s() + (1-fraction)*b->get_s();
        float tempT = fraction*a->get_t() + (1-fraction)*b->get_t();
        
        //the current distance from the original surface to the interpolated point between the parameter points A and B
        float currDisp = fraction*glm::length(m->getDisplacementValue(a->get_s(), a->get_t())) + (1-fraction)*glm::length(m->getDisplacementValue(b->get_s(), b->get_t()));
        //the actual displacement value at the new point
        float tempDisp = glm::length(m->getDisplacementValue(tempS, tempT)); //compare the difference between the interpolated displacement and the actual dispacement
        
        //printf("\t(tempDisp[%f] - currDisp[%f]) = %f\n", tempDisp, currDisp, (tempDisp - currDisp));
        
        if (std::fabs(tempDisp - currDisp) > maxHeight)
        {
            maxHeight = std::fabs(tempDisp - currDisp);
            pos = tempPos;
            s = tempS;
            t = tempT;
        }
        
    }
    //printf("MaxHeight: %f\n", maxHeight);
    
    v = addVertex(pos);
    v->setTextureCoordinates(s,t);
    setParentsChild(a,b,v);
    
    returnValue.newVert = v;
    
    if (std::fabs(maxHeight) < epsilon)
    {
        printf("EdgeVert: Negligible height difference, no vert should be created.\n");
        returnValue.significant = false;
    }
    else
    {
        returnValue.significant = true;
    }
    
    return returnValue;
}

Mesh::newVertReturn Mesh::AddMidVertexDisp(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *m)
{
    Mesh::newVertReturn returnValue;
    //find max disp value between points A and C and between points B and D, and use the largest
    
    glm::vec3 pos(0, 0, 0);
    float s = 0;
    float t = 0;
    float epsilon = 0.00001;
    
    float maxHeight = FLT_MIN;
    
    //go through the the points between A and C
    for (int i = 0; i < args->disp_increments; i++)
    {
        float fraction = (float)(i+1)/(float)(args->disp_increments + 1);
        glm::vec3 tempPos = fraction*a->get() + (1-fraction)*c->get();
        float tempS = fraction*a->get_s() + (1-fraction)*c->get_s();
        float tempT = fraction*a->get_t() + (1-fraction)*c->get_t();
        
        //the current distance from the original surface to the interpolated point between the parameter points A and B
        float currDisp = fraction*glm::length(m->getDisplacementValue(a->get_s(), a->get_t())) + (1-fraction)*glm::length(m->getDisplacementValue(c->get_s(), c->get_t()));
        //the actual displacement value at the new point
        float tempDisp = glm::length(m->getDisplacementValue(tempS, tempT)); //compare the difference between the interpolated displacement and the actual dispacement
        
        //printf("\t(tempDisp[%f] - currDisp[%f]) = %f\n", tempDisp, currDisp, (tempDisp - currDisp));
        
        if (std::fabs(tempDisp - currDisp) > maxHeight)
        {
            maxHeight = std::fabs(tempDisp - currDisp);
            pos = tempPos;
            s = tempS;
            t = tempT;
        }
        
    }
    
    //go through the points between B and D
    for (int i = 0; i < args->disp_increments; i++)
    {
        float fraction = (float)(i+1)/(float)(args->disp_increments + 1);
        glm::vec3 tempPos = fraction*b->get() + (1-fraction)*d->get();
        float tempS = fraction*b->get_s() + (1-fraction)*d->get_s();
        float tempT = fraction*b->get_t() + (1-fraction)*d->get_t();
        
        //the current distance from the original surface to the interpolated point between the parameter points A and B
        float currDisp = fraction*glm::length(m->getDisplacementValue(b->get_s(), b->get_t())) + (1-fraction)*glm::length(m->getDisplacementValue(d->get_s(), d->get_t()));
        //the actual displacement value at the new point
        float tempDisp = glm::length(m->getDisplacementValue(tempS, tempT)); //compare the difference between the interpolated displacement and the actual dispacement
        
        //printf("\t(tempDisp[%f] - currDisp[%f]) = %f\n", tempDisp, currDisp, (tempDisp - currDisp));
        
        if (std::fabs(tempDisp - currDisp) > maxHeight)
        {
            maxHeight = std::fabs(tempDisp - currDisp);
            pos = tempPos;
            s = tempS;
            t = tempT;
        }
        
    }
    
    //printf("MaxHeight: %f\n", maxHeight);
    
    Vertex *v = addVertex(pos);
    v->setTextureCoordinates(s,t);
    
    returnValue.newVert = v;
    
    if (std::fabs(maxHeight) < epsilon)
    {
        printf("EdgeVert: Negligible height difference, no vert should be created.\n");
        returnValue.significant = false;
    }
    else
    {
        returnValue.significant = true;
    }
    
    return returnValue;
}

//does the same thing as Subdivision(), but moves the created vertices to the largest displacement value along the edges they are created on
void Mesh::SubdivisionAltered()
{
    
    //if it's the firs time it's being subdivided
    bool first_subdivision = false;
    if (original_quads.size() == subdivided_quads.size()) {
        first_subdivision = true;
    }
    
    //get all of the subdivided faces
    std::vector<Face*> tmp = subdivided_quads;
    //delete the old list of them
    subdivided_quads.clear();
    
    //for all subdivided faces
    for (unsigned int i = 0; i < tmp.size(); i++)
    {
        printf("Face: %d\n", i);
        //make a temp face
        Face *f = tmp[i];
        
        //make some temp verts for that face
        Vertex *a = (*f)[0];
        a->setTextureCoordinates((*f)[0]->get_s(), (*f)[0]->get_t());
        
        Vertex *b = (*f)[1];
        b->setTextureCoordinates((*f)[1]->get_s(), (*f)[1]->get_t());
        
        Vertex *c = (*f)[2];
        c->setTextureCoordinates((*f)[2]->get_s(), (*f)[2]->get_t());
        
        Vertex *d = (*f)[3];
        d->setTextureCoordinates((*f)[3]->get_s(), (*f)[3]->get_t());
        //printf("A_ST: (%f %f) B_ST: (%f %f) C_ST: (%f %f) D_ST: (%f %f)", a->get_s(), a->get_t(), b->get_s(), b->get_t(), c->get_s(), c->get_t(), d->get_s(), d->get_t());
        
        
        Mesh::newVertReturn abTemp = AddEdgeVertexDisp(a,b, f->getMaterial());
        Mesh::newVertReturn bcTemp = AddEdgeVertexDisp(a,b, f->getMaterial());
        Mesh::newVertReturn cdTemp = AddEdgeVertexDisp(a,b, f->getMaterial());
        Mesh::newVertReturn daTemp = AddEdgeVertexDisp(a,b, f->getMaterial());
        Mesh::newVertReturn midTemp = AddEdgeVertexDisp(a,b, f->getMaterial());
        
        // add new vertices on the edges
        Vertex *ab = abTemp.newVert;// = new Vertex();
        Vertex *bc = bcTemp.newVert;// = new Vertex();
        Vertex *cd = cdTemp.newVert;//= new Vertex();
        Vertex *da = daTemp.newVert;// = new Vertex();
        // add new point in the middle of the patch
        Vertex *mid = midTemp.newVert;// = new Vertex();
        
        // copy the color and emission from the old patch to the new
        Material *material = f->getMaterial();
        if (!first_subdivision) {
            removeFaceEdges(f);
            delete f;
        }
        printf("Made it here1\n");
        //create the new verts, and if none of them have significant displacement differences, then don't subdivide this face
        if(!abTemp.significant &&
           !bcTemp.significant &&
           !cdTemp.significant &&
           !daTemp.significant &&
           !midTemp.significant)
        {
            //addSubdividedQuad(a, b, c, d, material);
            //continue;
        }
        printf("Made it here2\n");
        //printf("AB_ST: (%f %f) BC_ST: (%f %f) CD_ST: (%f %f) DA_ST: (%f %f)", ab->get_s(), ab->get_t(), bc->get_s(), bc->get_t(), cd->get_s(), cd->get_t(), da->get_s(), da->get_t());
        
        //make sure all the old edges are not null
        assert (getEdge(a,b) != NULL);
        assert (getEdge(b,c) != NULL);
        assert (getEdge(c,d) != NULL);
        assert (getEdge(d,a) != NULL);
        
        printf("Made it here3\n");
        // create the new faces
        addSubdividedQuad(a,ab,mid,da,material);
        addSubdividedQuad(b,bc,mid,ab,material);
        addSubdividedQuad(c,cd,mid,bc,material);
        addSubdividedQuad(d,da,mid,cd,material);
        
        printf("Made it here4\n");
        //make sure all the new edges exist
        assert (getEdge(a,ab) != NULL);
        assert (getEdge(ab,b) != NULL);
        assert (getEdge(b,bc) != NULL);
        assert (getEdge(bc,c) != NULL);
        assert (getEdge(c,cd) != NULL);
        assert (getEdge(cd,d) != NULL);
        assert (getEdge(d,da) != NULL);
        assert (getEdge(da,a) != NULL);
        printf("Made it here5\n");
    }
    
    printf("\n");
}

void Mesh::DisplacementSubdivision()
{
    //subdivide the mesh based on current subdivided mesh
    SubdivisionAltered();
}

