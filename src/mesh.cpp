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
#include "vbo_structs.h"

#include <cfloat>
#include <map>
#include <climits>


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
    if (disp_subdivided_quads.size() != original_quads.size()) {
        for (i = 0; i < disp_subdivided_quads.size(); i++) {
            Face *f = disp_subdivided_quads[i];
            removeFaceEdges(f);
            delete f;
        }
    }
    if (disp_subdivided_tris.size() > 0)
    {
        for (i = 0; i < disp_subdivided_tris.size(); i++)
        {
            Triangle *t = disp_subdivided_tris[i];
            removeTriEdges(t);//removeFaceEdges(f);
            delete t;
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

void Mesh::addFace(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material, enum FACE_TYPE face_type)
{
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
    } else if (face_type == FACE_TYPE_SUBDIVIDED) {
        //assert (face_type == FACE_TYPE_SUBDIVIDED);
        subdivided_quads.push_back(f);
    } else {//(face_type == FACE_TYPE_DISPLACEMENT){
        assert (face_type == FACE_TYPE_DISPLACEMENT);
        disp_subdivided_quads.push_back(f);
    }
    // if it's a light, add it to that list too
    if (glm::length(material->getEmittedColor()) > 0 && face_type == FACE_TYPE_ORIGINAL) {
        original_lights.push_back(f);
    }
}

Triangle* Mesh::addTri(Vertex *a, Vertex *b, Vertex *c, Material *material, enum TRI_TYPE tri_type)
{
    // create the triangle
    Triangle *t = new Triangle(material);// = new Triangle(material);
    // create the edges
    Edge *ea = new Edge(a,b,t);
    Edge *eb = new Edge(b,c,t);
    Edge *ec = new Edge(c,a,t);
    // point the triangle to one of its edges
    t->setEdge(ea);
    // connect the edges to each other
    ea->setNext(eb);
    eb->setNext(ec);
    ec->setNext(ea);
    // verify these edges aren't already in the mesh
    // (which would be a bug, or a non-manifold mesh)
    assert (edges.find(std::make_pair(a,b)) == edges.end());
    assert (edges.find(std::make_pair(b,c)) == edges.end());
    assert (edges.find(std::make_pair(c,a)) == edges.end());
    // add the edges to the master list
    edges[std::make_pair(a,b)] = ea;
    edges[std::make_pair(b,c)] = eb;
    edges[std::make_pair(c,a)] = ec;
    // connect up with opposite edges (if they exist)
    edgeshashtype::iterator ea_op = edges.find(std::make_pair(b,a));
    edgeshashtype::iterator eb_op = edges.find(std::make_pair(c,b));
    edgeshashtype::iterator ec_op = edges.find(std::make_pair(a,c));
    if (ea_op != edges.end()) { ea_op->second->setOpposite(ea); }
    if (eb_op != edges.end()) { eb_op->second->setOpposite(eb); }
    if (ec_op != edges.end()) { ec_op->second->setOpposite(ec); }
    
    //printf("EA: %p EB: %p EC: %p\n", ea->getOpposite(), eb->getOpposite(), ec->getOpposite());
    
    // add the triangle to the appropriate master list
    if (tri_type == TRI_TYPE_ORIGINAL)
    {
        original_tris.push_back(t);
    }
    else if(tri_type == TRI_TYPE_SUBDIVIDED)
    {
        t->setSubdivIndex(subdivided_tris.size());
        subdivided_tris.push_back(t);
    }
    else if(tri_type == TRI_TYPE_DISPLACEMENT)
    {
        disp_subdivided_tris.push_back(t);
    }
    
    assert (triangles.find(t->getID()) == triangles.end());
    triangles[t->getID()] = t;
    
    assert(t->getEdge() != NULL);
    assert((*t)[0] != NULL);
    assert((*t)[1] != NULL);
    assert((*t)[2] != NULL);
    
    
    //printf("\tEA: %p EB: %p EC: %p\n", ea->getOpposite(), eb->getOpposite(), ec->getOpposite());
    
    return t;
    
    //return t->getID();
    
    // if it's a light, add it to that list too //we don't support light triangles currently
    /*if (glm::length(material->getEmittedColor()) > 0)
    {
        original_lights.push_back(f);
    }*/
}

void Mesh::removeFaceEdges(Face *f) {
    // helper function for face deletion
    const Edge *ea = f->getEdge();
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

void Mesh::removeTriEdges(Triangle *t)
{
    
    // helper function for face deletion
    Edge *ea = t->getEdge();
    Edge *eb = ea->getNext();
    Edge *ec = eb->getNext();
    assert (ec->getNext() == ea);
    Vertex *a = ea->getStartVertex();
    Vertex *b = eb->getStartVertex();
    Vertex *c = ec->getStartVertex();
    // remove elements from master lists
    edges.erase(std::make_pair(a,b));
    edges.erase(std::make_pair(b,c));
    edges.erase(std::make_pair(c,a));
    // clean up memory
    delete ea;
    delete eb;
    delete ec;
}

void Mesh::removeSubdivideTriangle(Triangle *t)
{
    subdivided_tris.erase(subdivided_tris.begin() + t->getSubdivIndex());
    for (int i = t->getSubdivIndex(); i<subdivided_tris.size(); i++)
    {
        subdivided_tris[i]->setSubdivIndex(i);
    }
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
            
            ////////printf("\n\nToken: %s\n", token.c_str());
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
            
            ////////printf("DisplaceFile: %s\n", displace_file.c_str());
            
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
    
    if (camera == NULL) {
        // if not initialized, position a perspective camera and scale it so it fits in the window
        assert (bbox != NULL);
        glm::vec3 point_of_interest; bbox->getCenter(point_of_interest);
        float max_dim = bbox->maxDim();
        glm::vec3 camera_position = point_of_interest + glm::vec3(0,0,4*max_dim);
        glm::vec3 up = glm::vec3(0,1,0);
        camera = new PerspectiveCamera(camera_position, point_of_interest, up, 20 * M_PI/180.0);
    }
    ////printf("Finished camera initialization\n");
    //for subdivision using tris
    ConvertOriginalQuadsToTris();
    ////printf("Finished converting mesh to tris\n");
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
    //mesh_tri_verts_VBO  ************************************PUT THE VERT IN HERE**************
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
        
        assert (getEdge(ab,mid) != NULL);
        assert (getEdge(mid,ab) != NULL);
        
        assert (getEdge(bc,mid) != NULL);
        assert (getEdge(mid,bc) != NULL);
        
        assert (getEdge(cd,mid) != NULL);
        assert (getEdge(mid,cd) != NULL);
        
        assert (getEdge(da,mid) != NULL);
        assert (getEdge(mid,da) != NULL);
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
        
        //////printf("\t(tempDisp[%f] - currDisp[%f]) = %f\n", tempDisp, currDisp, (tempDisp - currDisp));
        
        if (std::fabs(tempDisp - currDisp) > maxHeight)
        {
            maxHeight = std::fabs(tempDisp - currDisp);
            pos = tempPos;
            s = tempS;
            t = tempT;
        }
        
    }
    //////printf("MaxHeight: %f\n", maxHeight);
    
    v = addVertex(pos);
    v->setTextureCoordinates(s,t);
    setParentsChild(a,b,v);
    
    returnValue.newVert = v;
    
    if (std::fabs(maxHeight) < epsilon)
    {
        ////printf("EdgeVert: Negligible height difference, no vert should be created.\n");
        returnValue.significant = false;
    }
    else
    {
        returnValue.significant = false;
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
        
        //////printf("\t(tempDisp[%f] - currDisp[%f]) = %f\n", tempDisp, currDisp, (tempDisp - currDisp));
        
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
        
        //////printf("\t(tempDisp[%f] - currDisp[%f]) = %f\n", tempDisp, currDisp, (tempDisp - currDisp));
        
        if (std::fabs(tempDisp - currDisp) > maxHeight)
        {
            maxHeight = std::fabs(tempDisp - currDisp);
            pos = tempPos;
            s = tempS;
            t = tempT;
        }
        
    }
    
    //////printf("MaxHeight: %f\n", maxHeight);
    
    Vertex *v = addVertex(pos);
    v->setTextureCoordinates(s,t);
    
    returnValue.newVert = v;
    
    if (std::fabs(maxHeight) < epsilon)
    {
        ////printf("EdgeVert: Negligible height difference, no vert should be created.\n");
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
        //////printf("Face: %d\n", i);
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
        //////printf("A_ST: (%f %f) B_ST: (%f %f) C_ST: (%f %f) D_ST: (%f %f)", a->get_s(), a->get_t(), b->get_s(), b->get_t(), c->get_s(), c->get_t(), d->get_s(), d->get_t());
        
        // add new vertices on the edges
        //Vertex *ab = AddEdgeVertexDisp(a,b, f->getMaterial());
        //Vertex *bc = AddEdgeVertexDisp(b,c, f->getMaterial());
        //Vertex *cd = AddEdgeVertexDisp(c,d, f->getMaterial());
        //Vertex *da = AddEdgeVertexDisp(d,a, f->getMaterial());
        // add new point in the middle of the patch
        //Vertex *mid = AddMidVertexDisp(ab,bc,cd,da, f->getMaterial());
        
        Mesh::newVertReturn abTemp = AddEdgeVertexDisp(a,b, f->getMaterial());
        Mesh::newVertReturn bcTemp = AddEdgeVertexDisp(b,c, f->getMaterial());
        Mesh::newVertReturn cdTemp = AddEdgeVertexDisp(c,d, f->getMaterial());
        Mesh::newVertReturn daTemp = AddEdgeVertexDisp(d,a, f->getMaterial());
        Vertex *ab = abTemp.newVert;// = new Vertex();
        Vertex *bc = bcTemp.newVert;// = new Vertex();
        Vertex *cd = cdTemp.newVert;//= new Vertex();
        Vertex *da = daTemp.newVert;// = new Vertex();
        
        Mesh::newVertReturn midTemp = AddMidVertexDisp(ab,bc,cd,da, f->getMaterial());
        // add new point in the middle of the patch
        Vertex *mid = midTemp.newVert;// = new Vertex();
        
        //////printf("AB_ST: (%f %f) BC_ST: (%f %f) CD_ST: (%f %f) DA_ST: (%f %f)", ab->get_s(), ab->get_t(), bc->get_s(), bc->get_t(), cd->get_s(), cd->get_t(), da->get_s(), da->get_t());
        
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
        //////printf("Made it here1\n");
        
        if(!abTemp.significant &&
           !bcTemp.significant &&
           !cdTemp.significant &&
           !daTemp.significant &&
           !midTemp.significant)
        {
            ////printf("This subdivided section is nearly useless");
            //addSubdividedQuad(a, b, c, d, material);
            //continue;
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
    }//*/
    
    ////printf("\n");
}

//this is called when 'D' is pressed
void Mesh::DisplacementSubdivision()
{
    bool first_subdivision = false;
    if (original_quads.size() == subdivided_quads.size()) {
        first_subdivision = true;
    }
    
    //subdivide the mesh based on current subdivided mesh
    SubdivisionAltered();
    //make a mesh with height values adjusted for the displacement from this subdivided mesh
    
    //first clear out old mesh
    disp_subdivided_quads.clear();
    //now copy over the subdivided mesh
    disp_subdivided_quads = subdivided_quads;
    
    //for each face, push it in the normal direction equal to the displacement value
    for (int i = 0; i<subdivided_quads.size(); i++)
    {
        //make a temp face
        Face *f = subdivided_quads[i];
        glm::vec3 norm = glm::normalize(f->computeNormal());
        Material* m = f->getMaterial();
        
        //make some temp verts for that face
        Vertex *a = addVertex((*f)[0]->get());
        a->setTextureCoordinates((*f)[0]->get_s(), (*f)[0]->get_t());
        float aDisp = glm::length(m->getDisplacementValue(a->get_s(), a->get_t()));
        
        Vertex *b = addVertex((*f)[1]->get());
        b->setTextureCoordinates((*f)[1]->get_s(), (*f)[1]->get_t());
        float bDisp = glm::length(m->getDisplacementValue(b->get_s(), b->get_t()));
        
        Vertex *c = addVertex((*f)[2]->get());
        c->setTextureCoordinates((*f)[2]->get_s(), (*f)[2]->get_t());
        float cDisp = glm::length(m->getDisplacementValue(c->get_s(), c->get_t()));
        
        Vertex *d = addVertex((*f)[3]->get());
        d->setTextureCoordinates((*f)[3]->get_s(), (*f)[3]->get_t());
        float dDisp = glm::length(m->getDisplacementValue(d->get_s(), d->get_t()));
        
        
        a->set(a->get() + norm*aDisp);
        b->set(b->get() + norm*bDisp);
        c->set(c->get() + norm*cDisp);
        d->set(d->get() + norm*dDisp);
        
        addDisplacementQuad(a, b, c, d, m);
    }//*/
}

// =================================================================
// DISPLACEMENT SUBDIVISION FOR TRIANGLES
// =================================================================
Triangle* Mesh::getOppositeTriangle(Triangle* currentTri, Vertex* a, Vertex* b)
{
    if (getEdge(a, b) == NULL)
    {
        printf("\tEdge %d-%d is NULL, switch it\n", a->getIndex(), b->getIndex());
        Vertex * temp = a;
        a = b;
        b = temp;
        if (getEdge(a, b) == NULL)
        {
            printf("\tCan't find edge between these two points... that's a problem.\n");
            return NULL;
        }
    }
    Edge* currEdge = getEdge(a, b);
    Edge* oppositeEdge = getEdge(b, a);
    
    //I was having some issues with edges not establishing opposites
    if(oppositeEdge != NULL && currEdge != NULL)
    {
        if (currEdge->getOpposite() == NULL)
        {
            currEdge->setOpposite(oppositeEdge);
        }
        if (oppositeEdge->getOpposite() == NULL)
        {
            oppositeEdge->setOpposite(currEdge);
        }
    }

    if (currEdge->getTriangle() == currentTri) //if the a-b edge is on this side, set the mate to the opposite triangle (the b-a side)
    {
        printf("1: Edge %d-%d is the current triangle\n", a->getIndex(), b->getIndex());
        
        //set this new triangles mate
        if (currEdge->getOpposite() == NULL)
        {
            printf("\tOpposite Edge is NULL\n");
            return NULL;
        }
        else
        {
            //assert(oppositeEdge == NULL);
            printf("\t%d-%d's triangle is %d\n", a->getIndex(), b->getIndex(), getEdge(a,b)->getOpposite()->getTriangle()->getID());
            //currentTri->setMateEdge(a, b);
            //currEdge->getOpposite()->getTriangle()->setMateEdge(b, a);
            return currEdge->getOpposite()->getTriangle();
            
        }
        
    }
    else //otherwise edge a-b is on the other side of this triangle, so just set mate to that triangle
    {
        printf("2: Edge %d-%d is the oppposite triangle\n", b->getIndex(), a->getIndex());
        if (currEdge == NULL)
        {
            printf("\tCurr Edge is NULL\n");
            return NULL;
        }
        else
        {
            printf("\t%d-%d's triangle is %d\n", b->getIndex(), a->getIndex(), getEdge(a,b)->getTriangle()->getID());
            //currentTri->setMateEdge(currEdge->getOpposite());
            //currEdge->getTriangle()->setMateEdge(a, b);
            return currEdge->getTriangle();
        }
    }
}

Triangle* Mesh::makeSubTriangle(Vertex *center, Vertex *a, Vertex* b, Triangle* mate, Triangle* tri)
{
    //Triangle* Mesh::addTri(Vertex *a, Vertex *b, Vertex *c, Material *material, enum TRI_TYPE tri_type)
    //make the triangle
    Triangle* subTri = addTri(center, a, b, tri->getMaterial(), TRI_TYPE_SUBDIVIDED);
    //set the generation of the new triangle
    subTri->setGeneration(tri->getGeneration()+1);
    subTri->setMate(mate);
    //subTri->setMateEdge(mate->getMateVert1(), mate->getMateVert2());
    //set the opposite triangles mate to this new triangle if the old one was the parent triangle
    if (mate != NULL && mate->getMate() == tri)
    {
        mate->setMate(subTri);
    }
    return subTri;
}

Edge* Mesh::findMateEdge(Triangle* t1, Triangle* t2)
{
    std::vector<Vertex*> mateVerts;
    for (int i = 0; i<3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if ((*t1)[i] == (*t2)[j])
            {
                mateVerts.push_back((*t1)[i]);
                j = 3;//essentially just break;
            }
        }
    }
    assert(mateVerts.size() == 2);
    return getEdge(mateVerts[0], mateVerts[1]);
}

bool Mesh::splitR3(Triangle* tri)
{
    //root3 split
    /*
     split(T)
        if (T.index is even) then
            compute midpoint P
            split T(A,B,C) into T[1](P,A,B),T[2](P,B,C),T[3](P,C,A)
            for i = 1,2,3 do
                T[i].index = T.index + 1
                if (T[i].mate[1].index == T[i].index) then
                swap(T[i],T[i].mate[1])
        else
            if (T.mate[1].index == T.index - 2)
                split(T.mate[1])

            split(T.mate[1])//  ... triggers edge swap
     */

    if (tri->getGeneration()%2 == 0)
    {
        ////printf("Tri has Even generation\n");
        //get the original verts we will be using
        Vertex* a = (*tri)[0];
        Vertex* b = (*tri)[1];
        Vertex* c = (*tri)[2];

        //get the center of the face
        Vertex* center = addVertex((a->get() + b->get() + c->get())/3.0f);

        //set the new S and T coords for this new face
        center->setTextureCoordinates((a->get_s() + b->get_s() + c->get_s())/3.0f,
                                     (a->get_t() + b->get_t() + c->get_t())/3.0f);

        ////printf("Verts Created\n");

        //first store some information on tri before deleting it
        //Material* mat = tri->getMaterial();

        Triangle* s1Mate = getOppositeTriangle(tri, a, b); //gets the mate triangle for s1
        Triangle* s2Mate = getOppositeTriangle(tri, b, c); //gets the mate triangle for s2
        Triangle* s3Mate = getOppositeTriangle(tri, c, a); //gets the mate triangle for s3
        
        printf("s1Mate: %p, s2Mate: %p, s3Mate: %p\n", s1Mate, s2Mate, s3Mate);
        
        ////printf("Got the mates\n");

        //need to delete the old triangle before making new ones, so the old edges don't
        //conflict with the new triangle edges
        removeTriEdges(tri);
        ////printf("Deleted Old Triangle edges\n");


        //make the new 3 triangles, and set their material, generation, and mate
        Triangle* s1 = makeSubTriangle(center, a, b, s1Mate, tri);
        Triangle* s2 = makeSubTriangle(center, b, c, s2Mate, tri);
        Triangle* s3 = makeSubTriangle(center, c, a, s3Mate, tri);
        
        
        
        ////printf("Made the small triangles\n");

        //change the mates' mates before deleting tri and swap if the mates are the same generation
        if(s1Mate != NULL && s1Mate->getMate() == s1)
        {
            s1Mate->setMate(s1);
            //printf("s1Mate generation = %d == s1 generation = %d\n", s1Mate->getGeneration(), s1->getGeneration());
            if (s1Mate->getGeneration() == s1->getGeneration())
            {
                swap(findMateEdge(s1, s1Mate));
            }
            
            //check if s1Mate leads to a border, then so does s1
            if (s1Mate->leadsToBorderTri())
            {
                s1->setLeadsToBorder(true);
            }
        }
        else
        {
            s1->setLeadsToBorder(true);
        }
        if(s2Mate != NULL && s2Mate->getMate() == s2)
        {
            //printf("s2Mate generation = %d == s2 generation = %d\n", s2Mate->getGeneration(), s2->getGeneration());
            s2Mate->setMate(s2);
            if (s2Mate->getGeneration() == s2->getGeneration())
            {
                swap(findMateEdge(s2, s2Mate));
            }
            
            
            //check if s2Mate leads to a border, then so does s2
            if (s2Mate->leadsToBorderTri())
            {
                s2->setLeadsToBorder(true);
            }
        }
        else
        {
            s2->setLeadsToBorder(true);
        }
        if(s3Mate != NULL && s3Mate->getMate() == s3)
        {
            //printf("s3Mate generation = %d == s3 generation = %d\n", s3Mate->getGeneration(), s3->getGeneration());
            s3Mate->setMate(s3);
            if (s3Mate->getGeneration() == s3->getGeneration())
            {
                swap(findMateEdge(s3, s3Mate));
            }
            
            
            //check if s3Mate leads to a border, then so does s3
            if (s3Mate->leadsToBorderTri())
            {
                s3->setLeadsToBorder(true);
            }
        }
        else
        {
            s2->setLeadsToBorder(true);
        }
        ////printf("Reassigned Mates\n");

        //now we can delete tri (if it isn't from the original tris)
        //if(tri->getGeneration() > 0)
        //{
        removeSubdivideTriangle(tri);
        ////printf("Deleted Old Triangle\n");
        //}

        /*
        for i = 1,2,3 do
            T[i].index = T.index + 1
            if (T[i].mate[1].index == T[i].index) then
                swap(T[i],T[i].mate[1])
         */

        return true;

    }
    else
    {
        if (tri->getMate() == NULL)
        {
            ////printf("Invalid Triangle, try again.\n");
            return false;
        }
        //split the mate if two generations behind
        if (tri->getMate()->getGeneration() == tri->getGeneration() - 2)
        {
            ////printf("Starting Recursion1\n");
            splitR3(tri->getMate());
        }
        assert(tri->getMate() != NULL);
        //split the mate unconditionally (always resulting in an odd generation result)
        ////printf("Starting Recursion2\n");
        splitR3(tri->getMate());
    }
    ////printf("Finished a split\n");
    //initialize VBOs (so that the screen updates)
    return true;
}


//root3 swap edges
/*
 swap(T1,T2)
    change T1(A,B,C), T2(B,A,D) into T1(C,A,D), T2(D,B,C)
    T1.index++
    T2.index++
 */
bool Mesh::swap(Edge* sharedEdge)//, Vertex *a, Vertex *b, Vertex *c, Vertex *d) // turn ABC and ABD into ACD and BCD
{
    //printf("Trying to swap\n");
    Vertex *a = sharedEdge->getEndVertex();
    Vertex *b = sharedEdge->getStartVertex();
    Vertex *c = sharedEdge->getOpposite()->getNext()->getEndVertex();
    Vertex *d = sharedEdge->getNext()->getEndVertex();
    //printf("Made the verts\n");
    Material *material = sharedEdge->getTriangle()->getMaterial();//this would not work if this edge was a boundary between 2 materials, that boundry would move
    //printf("Got the material\n");
    //store triangles for deletion
    Triangle *t1 = sharedEdge->getOpposite()->getTriangle();
    Triangle *t2 = sharedEdge->getTriangle();
    //printf("Reconstructed triangles from edge\n");
    //collect mates
    //Triangle *t1Mate = t1->getMate();
    //Triangle *t2Mate = t2->getMate();
    //printf("got the mates\n");
    
    int newGen = t1->getGeneration() + 1;
    
    //remove old edges
    removeTriEdges(t1);
    removeTriEdges(t2);
    //printf("Removed old edges\n");
    //delete t1;
    //delete t2;
    
    //delete triangles
    removeSubdivideTriangle(t1);
    removeSubdivideTriangle(t2);
    //printf("Removed old Triangles\n");
    t1 = NULL;
    t2 = NULL;
    
    //add triangles
    t1 = addTri(c, a, d, material, TRI_TYPE_SUBDIVIDED);
    t2 = addTri(d, b, c, material, TRI_TYPE_SUBDIVIDED);
    
    //printf("Added the new triangles\n");

    t1->setGeneration(newGen);
    t2->setGeneration(newGen);
    
    //printf("Set the new generations\n");
    
    return true;
    //set the mates of the new edges to what they used to be usable in future iterations

    /*if (t1Mate != NULL) //either has a or b as a vert shared with it's new destination
    {
        //printf("starting t1mate triangle mate assignment\n");
        Edge* t1Edge = t1Mate->getEdge();
        for (int i = 0; i < 3; i++)
        {
            if (t1Edge->getStartVertex() == a) //if a is in the mate's list of verts then it's part of t1
            {
                //printf("t1mate is part of the new t1\n");
                t1Mate->setMate(t1);
                break;
            }
            else if (t1Edge->getStartVertex() == b) //if b is in the mate's list of verts then it's part of t2
            {
                //printf("t1mate is part of the new t2\n");
                t1Mate->setMate(t2);
                break;
            }
            t1Edge = t1Edge->getNext();
        }
    }
    if (t2Mate != NULL) //either has a or b as a vert shared with it's new desination
    {
        //printf("starting t2mate triangle mate assignment\n");
        Edge* t2Edge = t2Mate->getEdge();
        //printf("Assigned t2Edge\n");
        assert(t2Edge->getStartVertex() != NULL);
        assert(t2Edge->getEndVertex() != NULL);
        for (int i = 0; i < 3; i++)          // *****************There is something wrong here****************
        {
            //printf("CurrVert: %d, a: %d, b: %d\n", t2Edge->getStartVertex()->getIndex(), a->getIndex(), b->getIndex());
            if (t2Edge->getStartVertex() == a) //if a is in the mate's list of verts then it's part of t1
            {
                //printf("t2mate is part of the new t1\n");
                t2Mate->setMate(t1);
                break;
            }
            else if (t2Edge->getStartVertex() == b) //if b is in the mate's list of verts then it's part of t2
            {
                //printf("t2mate is part of the new t2\n");
                t2Mate->setMate(t2);
                break;
            }
            t2Edge = t2Edge->getNext();
        }
        //printf("finished t2mate triangle mate assignment\n");
    }*/
}

void Mesh::DisplacementSubdivisionTriangles() //this method can be changed to only update disp_subdivided_tris
{
#if 1
    bool first_subdivision = false;
    if (original_tris.size() == subdivided_tris.size())
    {
        first_subdivision = true;
    }
    
    //subdivide the mesh based on current subdivided mesh
    Triangle *tri = subdivided_tris[args->getRand(0, subdivided_tris.size()-1)];
    
    //for now subdivide based on least subdivided area.
    int smallestGen = INT_MAX;
    for (int i = 0; i < subdivided_tris.size(); i++)
    {
        if (subdivided_tris[i]->getGeneration() < smallestGen && !subdivided_tris[i]->isBorderEdge() && !subdivided_tris[i]->leadsToBorderTri())
        {
            smallestGen = subdivided_tris[i]->getGeneration();
            tri = subdivided_tris[i];
        }
    }
    
    
    
    ////printf("Chose triangle to split\n");
    splitR3(tri);
    ////printf("Split triangle\n");
    //make a mesh with height values adjusted for the displacement from this subdivided mesh
    
    //first clear out old displacement mesh
    clearTriangles(&disp_subdivided_tris);

    //now copy over the subdivided mesh for the displacement mesh
    //disp_subdivided_tris = subdivided_tris;

    ////printf("Copied subdiv to disp vector\n");

    std::map<int, bool> vertMoved;

    assert(subdivided_tris.size() > 0);

    //for each vertex, push it in the normal direction equal to the displacement value
    for (int i = 0; i<subdivided_tris.size(); i++)
    {
        //make a temp face
        Triangle *tri = subdivided_tris[i];
        ////printf("set tri\n");
        assert(tri != NULL);
        ////printf("Normal Value: (%f %f %f)\n", tri->computeNormal().x, tri->computeNormal().y, tri->computeNormal().z);
        glm::vec3 norm = glm::normalize(tri->computeNormal());
        ////printf("set norm\n");
        Material* m = tri->getMaterial();
        ////printf("Got some info on the current triangle\n");
        
        //make some new verts for a new triangle
        Vertex *a = addVertex((*tri)[0]->get());
        a->setTextureCoordinates((*tri)[0]->get_s(), (*tri)[0]->get_t());
        float aDisp = glm::length(m->getDisplacementValue(a->get_s(), a->get_t()));
        
        Vertex *b = addVertex((*tri)[1]->get());
        b->setTextureCoordinates((*tri)[1]->get_s(), (*tri)[1]->get_t());
        float bDisp = glm::length(m->getDisplacementValue(b->get_s(), b->get_t()));
        
        Vertex *c = addVertex((*tri)[2]->get());
        c->setTextureCoordinates((*tri)[2]->get_s(), (*tri)[2]->get_t());
        float cDisp = glm::length(m->getDisplacementValue(c->get_s(), c->get_t()));

        ////printf("Made 3 new verts\n");

        if (!vertMoved[a->getIndex()])
        {
            vertMoved[a->getIndex()] = true;
            a->set(a->get() + norm*aDisp); //*********************OPPORTUNITY TO CHANGE HOW MUCH DISPLACEMENT********
        }
        if (!vertMoved[b->getIndex()])
        {
            vertMoved[b->getIndex()] = true;
            b->set(b->get() + norm*bDisp); //*********************OPPORTUNITY TO CHANGE HOW MUCH DISPLACEMENT********
        }
        if (!vertMoved[c->getIndex()])
        {
            vertMoved[c->getIndex()] = true;
            c->set(c->get() + norm*cDisp); //*********************OPPORTUNITY TO CHANGE HOW MUCH DISPLACEMENT********
        }
        ////printf("Moved those verts\n");

        addTri(a, b, c, m, TRI_TYPE_DISPLACEMENT);
        ////printf("Added the triangle to the displacement vector\n");
    }//*/

    ////printf("Finished the subdivision\n");
    
#endif
}


void Mesh::ConvertOriginalQuadsToTris()
{
    clearNeighbors();
    clearTriangles(&subdivided_tris);
    clearTriangles(&disp_subdivided_tris);
    clearTriangles(&original_tris);
    
    for (int i = 0; i < original_quads.size(); i++)
    {
        //if the face is not a light source
        if (glm::length(original_quads[i]->getMaterial()->getEmittedColor()) == 0)
        {
            //printf("here1\n");
            
            //store the face verts
            Vertex *quadVertA = (*original_quads[i])[0];
            Vertex *quadVertB = (*original_quads[i])[1];
            Vertex *quadVertC = (*original_quads[i])[2];
            Vertex *quadVertD = (*original_quads[i])[3];
            
            //printf("here2\n");
            
            //make some verts for the new triangles
            Vertex *a = NULL;// = addVertex(quadVertA->get());
            Vertex *b = NULL;// = addVertex(quadVertB->get());
            Vertex *c = NULL;// = addVertex(quadVertC->get());
            Vertex *d = NULL;// = addVertex(quadVertD->get());
            
            //printf("here3\n");
            
            //if the a-b edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles a-b vertices are created already
            if (getEdge(quadVertA, quadVertB)->getOpposite() != NULL && getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdge() != NULL)
            {
                printf("A-B Existed already\n");
                if (a == NULL)
                {
                    a = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (b == NULL)
                {
                    b = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            //if the b-c edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles b-c vertices are created already
            if (getEdge(quadVertB, quadVertC)->getOpposite() != NULL && getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdge() != NULL)
            {
                printf("B-C Existed already\n");
                if (b == NULL)
                {
                    b = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (c == NULL)
                {
                    c = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            //if the c-d edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles c-d vertices are created already
            if (getEdge(quadVertC, quadVertD)->getOpposite() != NULL && getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdge() != NULL)
            {
                printf("C-D Existed already\n");
                if (c == NULL)
                {
                    c = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (d == NULL)
                {
                    d = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            //if the d-a edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles d-a vertices are created already
            if (getEdge(quadVertD, quadVertA)->getOpposite() != NULL && getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdge() != NULL)
            {
                printf("D-A Existed already\n");
                if (d == NULL)
                {
                    d = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (a == NULL)
                {
                    a = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            
            //printf("here4\n");
            //if a wasn't found above, make a new vert for it
            if (a == NULL)
            {
                printf("A Hasn't Been Created Yet\n");
                a = addVertex(quadVertA->get());
            }
            //if b wasn't found above, make a new vert for it
            if (b == NULL)
            {
                printf("B Hasn't Been Created Yet\n");
                b = addVertex(quadVertB->get());
            }
            //if c wasn't found above, make a new vert for it
            if (c == NULL)
            {
                printf("C Hasn't Been Created Yet\n");
                c = addVertex(quadVertC->get());
            }
            //if d wasn't found above, make a new vert for it
            if (d == NULL)
            {
                printf("D Hasn't Been Created Yet\n");
                d = addVertex(quadVertD->get());
            }

            //printf("here5\n");
            //make triangle 1 of the quad (abc)
            addTri(a, b, c, original_quads[i]->getMaterial(), TRI_TYPE_ORIGINAL);
            
            //make triangle 2 of the quad (acd)
            addTri(a, c, d, original_quads[i]->getMaterial(), TRI_TYPE_ORIGINAL);
            
            //printf("here6\n");
            //set the triangleEdges of the quad's edges
            getEdge(quadVertA, quadVertB)->setTriangleEdge(getEdge(a, b));
            getEdge(quadVertB, quadVertC)->setTriangleEdge(getEdge(b, c));
            getEdge(quadVertC, quadVertD)->setTriangleEdge(getEdge(c, d));
            getEdge(quadVertD, quadVertA)->setTriangleEdge(getEdge(d, a));
                
            

            //printf("here7\n");
            //Also add to subdiv_tri
            Vertex *aS = NULL;// = addVertex(quadVertA->get());
            Vertex *bS = NULL;// = addVertex(quadVertB->get());
            Vertex *cS = NULL;// = addVertex(quadVertC->get());
            Vertex *dS = NULL;// = addVertex(quadVertD->get());
            
            //printf("here8\n");
            if (getEdge(quadVertA, quadVertB)->getOpposite() != NULL && getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                aS = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                bS = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
            }
            //printf("here8.1\n");
            if (getEdge(quadVertB, quadVertC)->getOpposite() != NULL && getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                //printf("here8.11 %p\n", getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv());
                if (bS == NULL)
                {
                    bS = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                }
                //printf("here8.12\n");
                cS = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
            }
            //printf("here8.2\n");
            if (getEdge(quadVertC, quadVertD)->getOpposite() != NULL && getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                if (cS == NULL)
                {
                    cS = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                }
                dS = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
            }
            //printf("here8.3\n");
            if (getEdge(quadVertD, quadVertA)->getOpposite() != NULL && getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                if (dS == NULL)
                {
                    dS = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                }
                if (aS == NULL)
                {
                    aS = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
                }
            }
            
            //printf("here9\n");
            if (aS == NULL)
            {
                aS = addVertex(quadVertA->get());
            }
            if (bS == NULL)
            {
                bS = addVertex(quadVertB->get());
            }
            if (cS == NULL)
            {
                cS = addVertex(quadVertC->get());
            }
            if (dS == NULL)
            {
                dS = addVertex(quadVertD->get());
            }
            
            //printf("here10\n");
            //make triangle 1 of the quad (abc)
            addTri(aS, bS, cS, subdivided_quads[i]->getMaterial(), TRI_TYPE_SUBDIVIDED);
            
            //make triangle 2 of the quad (acd)
            addTri(aS, cS, dS, subdivided_quads[i]->getMaterial(), TRI_TYPE_SUBDIVIDED);
            
            //printf("here11\n");
            //set the triangleEdges of the quad's edges
            getEdge(quadVertA, quadVertB)->setTriangleEdgeSubdiv(getEdge(aS, bS));
            getEdge(quadVertB, quadVertC)->setTriangleEdgeSubdiv(getEdge(bS, cS));
            getEdge(quadVertC, quadVertD)->setTriangleEdgeSubdiv(getEdge(cS, dS));
            getEdge(quadVertD, quadVertA)->setTriangleEdgeSubdiv(getEdge(dS, aS));
        }
    }
    
    addAllNeighbors();
    printNeighbors();

    //initialize->setup->helper->
}

void Mesh::ConvertSubdividedQuadsToTris()// ********* Maybe make a recursive version of this that makes tris out of a quad, then tries to recurse on the quads around it
{
    clearTriangles(&subdivided_tris);
    clearTriangles(&disp_subdivided_tris);
    clearTriangles(&original_tris);
    //subdivided_tris.clear();
    //disp_subdivided_tris.clear();
    
    clearNeighbors();
    for (int i = 0; i < subdivided_quads.size(); i++)
    {
        //if the face is not a light source
        if (glm::length(subdivided_quads[i]->getMaterial()->getEmittedColor()) == 0)
        {
            //printf("_________________________\n");
            //printf("SHere1\n");
            Vertex *quadVertA = (*subdivided_quads[i])[0];
            Vertex *quadVertB = (*subdivided_quads[i])[1];
            Vertex *quadVertC = (*subdivided_quads[i])[2];
            Vertex *quadVertD = (*subdivided_quads[i])[3];
            
            //printf("SHere2\n");
            Vertex *a = NULL;// = addVertex(quadVertA->get());
            Vertex *b = NULL;// = addVertex(quadVertB->get());
            Vertex *c = NULL;// = addVertex(quadVertC->get());
            Vertex *d = NULL;// = addVertex(quadVertD->get());
            
            //if the a-b edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles a-b vertices are created already
            if (getEdge(quadVertA, quadVertB)->getOpposite() != NULL && getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdge() != NULL)
            {
                //printf("A-B Existed already\n");
                if (a == NULL)
                {
                    a = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (b == NULL)
                {
                    b = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            //if the b-c edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles b-c vertices are created already
            if (getEdge(quadVertB, quadVertC)->getOpposite() != NULL && getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdge() != NULL)
            {
                //printf("B-C Existed already\n");
                if (b == NULL)
                {
                    b = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (c == NULL)
                {
                    c = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            //if the c-d edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles c-d vertices are created already
            if (getEdge(quadVertC, quadVertD)->getOpposite() != NULL && getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdge() != NULL)
            {
                //printf("C-D Existed already\n");
                if (c == NULL)
                {
                    c = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (d == NULL)
                {
                    d = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            //if the d-a edge of the face has an opposite, check to see if a triangle edge already exists, if so the triangles d-a vertices are created already
            if (getEdge(quadVertD, quadVertA)->getOpposite() != NULL && getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdge() != NULL)
            {
                //printf("D-A Existed already\n");
                if (d == NULL)
                {
                    d = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdge()->getEndVertex();
                }
                if (a == NULL)
                {
                    a = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdge()->getStartVertex();
                }
            }
            
            //printf("here4\n");
            //if a wasn't found above, make a new vert for it
            if (a == NULL)
            {
                //printf("A Hasn't Been Created Yet\n");
                a = addVertex(quadVertA->get());
            }
            //if b wasn't found above, make a new vert for it
            if (b == NULL)
            {
                //printf("B Hasn't Been Created Yet\n");
                b = addVertex(quadVertB->get());
            }
            //if c wasn't found above, make a new vert for it
            if (c == NULL)
            {
                //printf("C Hasn't Been Created Yet\n");
                c = addVertex(quadVertC->get());
            }
            //if d wasn't found above, make a new vert for it
            if (d == NULL)
            {
                //printf("D Hasn't Been Created Yet\n");
                d = addVertex(quadVertD->get());
            }
            
            //printf("SHere5\n");
            //make triangle 1 of the quad (abc)
            addTri(a, b, c, subdivided_quads[i]->getMaterial(), TRI_TYPE_ORIGINAL);
            
            //make triangle 2 of the quad (acd)
            addTri(a, c, d, subdivided_quads[i]->getMaterial(), TRI_TYPE_ORIGINAL);
            
            //printf("SHere6\n");
            //set the triangleEdges of the quad's edges
            getEdge(quadVertA, quadVertB)->setTriangleEdge(getEdge(a, b));
            getEdge(quadVertB, quadVertC)->setTriangleEdge(getEdge(b, c));
            getEdge(quadVertC, quadVertD)->setTriangleEdge(getEdge(c, d));
            getEdge(quadVertD, quadVertA)->setTriangleEdge(getEdge(d, a));
            
            //printf("SHere7\n");
            //Also add to subdiv_tri
            Vertex *aS = nullptr;// = NULL;// = addVertex(quadVertA->get());
            Vertex *bS = nullptr;// = NULL;// = addVertex(quadVertB->get());
            Vertex *cS = nullptr;// = NULL;// = addVertex(quadVertC->get());
            Vertex *dS = nullptr;// = NULL;// = addVertex(quadVertD->get());
            
            //Edge * fe = subdivided_quads[i]->getEdge();
            
            //printf("SHere8\n");
            if (getEdge(quadVertA, quadVertB)->getOpposite() != NULL && getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                //printf("As-Bs Existed already\n");
                aS = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                bS = getEdge(quadVertA, quadVertB)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
            }
            //printf("SHere8.1\n");
            if (getEdge(quadVertB, quadVertC)->getOpposite() != NULL && getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                //printf("Bs-Cs Existed already\n");
                //printf("SHere8.11 %p\n", getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv());
                if (bS == NULL)
                {
                    bS = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                }
                //printf("SHere8.12\n");
                if (cS == NULL)
                {
                    cS = getEdge(quadVertB, quadVertC)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
                }
            }
            //printf("SHere8.2\n");
            if (getEdge(quadVertC, quadVertD)->getOpposite() != NULL && getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                //printf("Cs-Ds Existed already\n");
                if (cS == NULL)
                {
                    cS = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                }
                if (dS == NULL)
                {
                    dS = getEdge(quadVertC, quadVertD)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
                }
            }
            //printf("SHere8.3\n");
            if (getEdge(quadVertD, quadVertA)->getOpposite() != NULL && getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdgeSubdiv() != NULL)
            {
                //printf("Ds-As Existed already\n");
                if (dS == NULL)
                {
                    dS = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdgeSubdiv()->getEndVertex();
                }
                if (aS == NULL)
                {
                    aS = getEdge(quadVertD, quadVertA)->getOpposite()->getTriangleEdgeSubdiv()->getStartVertex();
                }
            }
            
            //printf("SHere9\n");
            if (aS == NULL)
            {
                //printf("As Hasn't Been Created Yet\n");
                aS = addVertex(quadVertA->get());
            }
            if (bS == NULL)
            {
                //printf("Bs Hasn't Been Created Yet\n");
                bS = addVertex(quadVertB->get());
            }
            if (cS == NULL)
            {
                //printf("Cs Hasn't Been Created Yet\n");
                cS = addVertex(quadVertC->get());
            }
            if (dS == NULL)
            {
                //printf("Ds Hasn't Been Created Yet\n");
                dS = addVertex(quadVertD->get());
            }
            
            assert(aS != NULL);
            assert(bS != NULL);
            assert(cS != NULL);
            assert(dS != NULL);
            
            //printf("SHere10\n");
            //make triangle 1 of the quad (abc)
            addTri(aS, bS, cS, subdivided_quads[i]->getMaterial(), TRI_TYPE_SUBDIVIDED);
            //printf("SHere10.5\n");
            //make triangle 2 of the quad (acd)
            addTri(aS, cS, dS, subdivided_quads[i]->getMaterial(), TRI_TYPE_SUBDIVIDED);
            
            //printf("SHere11\n");
            //set the triangleEdges of the quad's edges
            getEdge(quadVertA, quadVertB)->setTriangleEdgeSubdiv(getEdge(aS, bS));
            //printf("SHere11.1\n");
            getEdge(quadVertB, quadVertC)->setTriangleEdgeSubdiv(getEdge(bS, cS));
            //printf("SHere11.2\n");
            getEdge(quadVertC, quadVertD)->setTriangleEdgeSubdiv(getEdge(cS, dS));
            //printf("SHere11.3\n");
            getEdge(quadVertD, quadVertA)->setTriangleEdgeSubdiv(getEdge(dS, aS));
            //printf("SHere11.4\n");
            //printf("Finished converting quad %d/%d\n", i+1, (int)subdivided_quads.size());
        }
        else
        {
            //printf("Skipped quad %d/%d, it was a light.\n", i+1, (int)subdivided_quads.size());
        }
    }
    
    //printf("%d subdiv tris currently\n", (int)subdivided_tris.size());
    //lets check that everything is in working order
    for (int i = 0; i < subdivided_tris.size(); i++)
    {
        //printf("Triangle %d\n", subdivided_tris[i]->getSubdivIndex());
        assert(subdivided_tris[i]->getEdge() != NULL);
        Edge * e = subdivided_tris[i]->getEdge();
        for (int j = 0; j<3; j++)
        {
            //printf("Checking edge %d of triangle %d\n",j+1, i+1);
            //check edges
            if (e->getOpposite() == NULL) //if the opposite edge is null, make sure an opposite edge doesn't actually exist
            {
                //printf("\tOpposite of Edge %d is null\n", j+1);
                assert(getEdge(e->getEndVertex(), e->getStartVertex()) == NULL);
            }
            else //otherwise make sure the opposite edge has the same verts and such.
            {
                //printf("\tOpposite of Edge %d is not null\n", j+1);
                //make sure this is the right edge for oppositeEdge()
                assert(getEdge(e->getEndVertex(), e->getStartVertex()) == e->getOpposite());
                //make sure it has a link back to the current edge
                assert(e->getOpposite()->getOpposite() == e);
                //make sure the verts are the same
                assert(e->getOpposite()->getStartVertex() == e->getEndVertex());
                assert(e->getOpposite()->getEndVertex() == e->getStartVertex());
                
            }
            
            //check next edge
            assert(e->getNext() != NULL);
            e = e->getNext();
            
            //check verts
            assert((*subdivided_tris[i])[j] != NULL);
        }
    }
    
    
    //printf("SHere11.45\n");
    //addAllNeighbors();
    //printf("SHere11.5\n");
    //printNeighbors();
    //printf("SHere12\n");
    //initialize->setup->helper->
}

// boundary edges are red, crease edges are yellow
glm::vec3 EdgeColor(Edge *e) {
    if (e->getOpposite() == NULL) {
        return glm::vec3(1,0,0);
    }
    else
    {
        return glm::vec3(0,0,0.0);
    }
}

glm::vec3 ComputeNormal(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
    glm::vec3 v12 = p2;
    v12 -= p1;
    glm::vec3 v23 = p3;
    v23 -= p2;
    glm::vec3 normal = glm::cross(v12,v23);
    //std::cout << "\t\t\tFaceNormal: " << normal.x << ", " << normal.y << ", " << normal.z << "\n";
    normal = glm::normalize(normal);
    //std::cout << "\t\t\tFaceNormal: " << normal.x << ", " << normal.y << ", " << normal.z << "\n";
    return normal;
}

void Mesh::TriVBOHelper( std::vector<glm::vec3> &indexed_verts,
                        std::vector<unsigned int> &mesh_tri_indices,
                        const glm::vec3 &pos_a,
                        const glm::vec3 &pos_b,
                        const glm::vec3 &pos_c,
                        const glm::vec3 &normal_a,
                        const glm::vec3 &normal_b,
                        const glm::vec3 &normal_c,
                        const glm::vec3 &color_ab,
                        const glm::vec3 &color_bc,
                        const glm::vec3 &color_ca) {

    /*
     // To create a wireframe rendering...
     // Each mesh triangle is actually rendered as 3 small triangles
     //           b
     //          /|\
     //         / | \
     //        /  |  \
     //       /   |   \
     //      /    |    \
     //     /    .'.    \
     //    /  .'     '.  \
     //   /.'           '.\
     //  a-----------------c
     //
     */

    // the center is white, the colors of the two vertices depend on
    // whether the edge is a boundary edge (red) or crease edge (yellow)
    glm::vec3 center_color(1,1,1);
    // use simple averaging to find centroid & average normal
    glm::vec3 centroid = 1.0f / 3.0f * (pos_a + pos_b + pos_c);
    glm::vec3 normal = normal_a + normal_b + normal_c;
    glm::normalize(normal);

    int i = indexed_verts.size()/3;

    if (args->wireframe) {
        // WIREFRAME

        // make the 3 small triangles
        indexed_verts.push_back(pos_a);
        indexed_verts.push_back(normal_a);
        indexed_verts.push_back(color_ab);
        indexed_verts.push_back(pos_b);
        indexed_verts.push_back(normal_b);
        indexed_verts.push_back(color_ab);
        indexed_verts.push_back(centroid);
        indexed_verts.push_back(normal);
        indexed_verts.push_back(center_color);

        indexed_verts.push_back(pos_b);
        indexed_verts.push_back(normal_b);
        indexed_verts.push_back(color_bc);
        indexed_verts.push_back(pos_c);
        indexed_verts.push_back(normal_c);
        indexed_verts.push_back(color_bc);
        indexed_verts.push_back(centroid);
        indexed_verts.push_back(normal);
        indexed_verts.push_back(center_color);

        indexed_verts.push_back(pos_c);
        indexed_verts.push_back(normal_c);
        indexed_verts.push_back(color_ca);
        indexed_verts.push_back(pos_a);
        indexed_verts.push_back(normal_a);
        indexed_verts.push_back(color_ca);
        indexed_verts.push_back(centroid);
        indexed_verts.push_back(normal);
        indexed_verts.push_back(center_color);

        // add all of the triangle vertices to the indices list
        for (int j = 0; j < 9; j++) {
            mesh_tri_indices.push_back(i+j);
        }
    } else {
        // NON WIREFRAME
        // Note: gouraud shading with the mini triangles looks bad... :(

        // make the 3 small triangles
        indexed_verts.push_back(pos_a);
        indexed_verts.push_back(normal_a);
        indexed_verts.push_back(center_color);
        indexed_verts.push_back(pos_b);
        indexed_verts.push_back(normal_b);
        indexed_verts.push_back(center_color);
        indexed_verts.push_back(pos_c);
        indexed_verts.push_back(normal_c);
        indexed_verts.push_back(center_color);

        // add all of the triangle vertices to the indices list
        for (int j = 0; j < 3; j++) {
            mesh_tri_indices.push_back(i+j);
        }
    }
    
}

void Mesh::initializeVBOs() {
    // create a pointer for the vertex & index VBOs
    glGenBuffers(1, &mesh_tri_verts_VBO);
    glGenBuffers(1, &mesh_tri_indices_VBO);
    glGenBuffers(1, &mesh_textured_tri_indices_VBO);
}

void Mesh::setupVBOs() //Radiosity setupVBOs copy pasta
{
    HandleGLError("enter radiosity setupVBOs()");
    mesh_tri_verts.clear();
    mesh_tri_indices.clear();
    mesh_textured_tri_indices.clear();
    ////////printf("here1\n");

    std::vector<Triangle*> trianglesToShow;
    if (disp_subdivided_tris.size() > 0)
    {
        trianglesToShow = disp_subdivided_tris;
    }
    else if(subdivided_tris.size() > 0)
    {
        trianglesToShow = subdivided_tris;
    }
    else
    {
        trianglesToShow = original_tris;
    }
    
    printf("Got the list of things to show\n"); // ****************************THE ISSUE IS IN HERE SOMEWHERE********************
    // initialize the data in each vector
    assert (trianglesToShow.size() > 0);
    for (int i = 0; i < trianglesToShow.size(); i++)
    {
        printf("Me?.1\n");
        Triangle *triangle = trianglesToShow[i];
        printf("Me?.25\n");
        Edge *e = triangle->getEdge();
        printf("Me?.5\n");
        glm::vec3 normal = triangle->computeNormal();
        printf("Me?1\n");
        double avg_s = 0;
        double avg_t = 0;
        glm::vec3 avg_color(0,0,0);
        
        int start = mesh_tri_verts.size();
        
        // wireframe is normally black, except when it's the special
        // patch, then the wireframe is red
        glm::vec4 wireframe_color(0,0,0,0.5);
        //if (args->render_mode == RENDER_FORM_FACTORS && i == max_undistributed_patch) {
        //    wireframe_color = glm::vec4(1,0,0,1);
        //}
        ////////printf("here1.5\n");
        // add the 3 corner vertices
        for (int j = 0; j < 3; j++)
        {
            glm::vec3 pos = ((*triangle)[j])->get();
            double s = (*triangle)[j]->get_s();
            double t = (*triangle)[j]->get_t();
            printf("Me?2\n");
            ////////printf("here1.6\n");
            //glm::vec3 color = setupHelperForColor(t,i,j);
            //color = glm::vec3(linear_to_srgb(color.r),
            //                  linear_to_srgb(color.g),
            //                 linear_to_srgb(color.b));
            
            glm::vec3 color;
            if (triangle->getMaterial() != NULL)
            {
                color = triangle->getMaterial()->getDiffuseColor();
            }
            else
            {
                color = glm::vec3(.8, .8, .8);
            }
            ////////printf("here1.65\n");
            avg_color += (1.0f/3.0f) * color;
            ////////printf("here1.7\n");
            mesh_tri_verts.push_back(VBOPosNormalColor(pos,normal,
                                                       glm::vec4(color.r,color.g,color.b,1.0),
                                                       wireframe_color,
                                                       s,t));
            printf("Me?3\n");
            ////////printf("here1.75\n");
            avg_s += (1.0f/3.0f) * s;
            avg_t += (1.0f/3.0f) * t;
            e = e->getNext();
        }
        ////////printf("here2\n");
        
        // the centroid (for wireframe rendering)
        glm::vec3 centroid = triangle->computeCentroid();
        mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,
                                                   glm::vec4(avg_color.r,avg_color.g,avg_color.b,1),
                                                   glm::vec4(1,1,1,1),
                                                   avg_s,avg_t));
        printf("Me?4\n");
        ////////printf("here3\n");
        
        if (triangle->getMaterial()->hasTextureMap()) {
            mesh_textured_tri_indices.push_back(VBOIndexedTri(start+0,start+1,start+3));
            mesh_textured_tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+3));
            mesh_textured_tri_indices.push_back(VBOIndexedTri(start+2,start+0,start+3));
        } else {
            mesh_tri_indices.push_back(VBOIndexedTri(start+0,start+1,start+3));
            mesh_tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+3));
            mesh_tri_indices.push_back(VBOIndexedTri(start+2,start+0,start+3));
        }
        printf("Me?5\n");
        ////////printf("here4\n");
    }
    //assert ((int)mesh_tri_verts.size() == num_tris*5); //should this be a 4?
    //assert ((int)mesh_tri_indices.size() + (int)mesh_textured_tri_indices.size() == num_tris*4); //should this be a 3?
    
    printf("Me?5.1\n");
    // copy the data to each VBO
    glBindBuffer(GL_ARRAY_BUFFER,mesh_tri_verts_VBO);
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(VBOPosNormalColor) * trianglesToShow.size() * 5, //should this be a 4?
                 &mesh_tri_verts[0],
                 GL_STATIC_DRAW);
    printf("Me?5.25\n");
    ////////printf("here5\n");
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_tri_indices_VBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(VBOIndexedTri) * mesh_tri_indices.size(),
                 &mesh_tri_indices[0], GL_STATIC_DRAW);
    printf("Me?5.75\n");
    ////////printf("here6\n");
    if (mesh_textured_tri_indices.size() > 0) {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_textured_tri_indices_VBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     sizeof(VBOIndexedTri) * mesh_textured_tri_indices.size(),
                     &mesh_textured_tri_indices[0], GL_STATIC_DRAW);
        printf("Me?6\n");
        ////////printf("here7\n");
        
    }
    
    HandleGLError("radiosity setupVBOs() just before texture");
    
    // WARNING: this naive VBO implementation only allows a single texture
    // FIXME: something still buggy about textures
    int num_textured_materials = 0;
    for (unsigned int mat = 0; mat < materials.size(); mat++) {
        Material *m = materials[mat];
        if (m->hasTextureMap()) {
            // FIXME: old gl...
            //glBindTexture(GL_TEXTURE_2D,m->getTextureID());
            num_textured_materials++;
        }
    }
    printf("Me?7\n");
    assert (num_textured_materials <= 1);
    
    printf("Finished Setting up VBOs\n");
    
    HandleGLError("leave radiosity setupVBOs()");
}

void Mesh::drawVBOs() //Radiosity drawVBO
{
    
    // =====================
    // DRAW ALL THE POLYGONS
    
    //assert ((int)mesh_tri_indices.size() + (int)mesh_textured_tri_indices.size() == num_tris*3);
    
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

void Mesh::cleanupVBOs() //Radiosity Cleanup
{
    glDeleteBuffers(1, &mesh_tri_verts_VBO);
    glDeleteBuffers(1, &mesh_tri_indices_VBO);
    glDeleteBuffers(1, &mesh_textured_tri_indices_VBO);
    
    glDeleteTextures(1, &GLCanvas::textureID);
}

void CollectTrisWithVertex(Vertex *have, Triangle *tri, std::vector<Triangle*> &tris)
{
    for (unsigned int i = 0; i < tris.size(); i++)
    {
        if (tris[i] == tri) return;
    }
    if (have != (*tri)[0] && have != (*tri)[1] && have != (*tri)[2] && have != (*tri)[3]) return;
    tris.push_back(tri);
    for (int i = 0; i < 3; i++)
    {
        Edge *ea = tri->getEdge()->getOpposite();
        Edge *eb = tri->getEdge()->getNext()->getOpposite();
        Edge *ec = tri->getEdge()->getNext()->getNext()->getOpposite();
        if (ea != NULL) CollectTrisWithVertex(have,ea->getTriangle(),tris);
        if (eb != NULL) CollectTrisWithVertex(have,eb->getTriangle(),tris);
        if (ec != NULL) CollectTrisWithVertex(have,ec->getTriangle(),tris);
    }
}






#if 0
void Mesh::drawVBOs(const glm::mat4 &ProjectionMatrix,const glm::mat4 &ViewMatrix,const glm::mat4 &ModelMatrix) {
    HandleGLError("enter drawVBOs");
    
    // prepare data to send to the shaders
    glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;
    
    glm::vec3 lightPos = glm::vec3(4,4,4);
    HandleGLError("Error1");
    glUniform3f(GLCanvas::LightID, lightPos.x, lightPos.y, lightPos.z); //This one is possible, has the same error
    HandleGLError("Error2");
    glUniformMatrix4fv(GLCanvas::MatrixID, 1, GL_FALSE, &MVP[0][0]);
    HandleGLError("Error3");
    glUniformMatrix4fv(GLCanvas::ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
    HandleGLError("Error4");
    glUniformMatrix4fv(GLCanvas::ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);
    HandleGLError("Error5");
    glUniform1i(GLCanvas::wireframeID, args->wireframe);
    HandleGLError("Error6");
    
    /* all the above gl commands
     GL_INVALID_OPERATION is generated if there is no current program object.
     
     GL_INVALID_OPERATION is generated if the size of the uniform variable declared in the shader does not
     match the size indicated by the glUniform command.
     
     GL_INVALID_OPERATION is generated if one of the signed or unsigned integer variants of this function
     is used to load a uniform variable of type float, vec2, vec3, vec4, or an array of these, or if one
     of the floating-point variants of this function is used to load a uniform variable of type int,
     ivec2, ivec3, ivec4, unsigned int, uvec2, uvec3, uvec4, or an array of these.
     
     GL_INVALID_OPERATION is generated if one of the signed integer variants of this function is used to
     load a uniform variable of type unsigned int, uvec2, uvec3, uvec4, or an array of these.
     
     GL_INVALID_OPERATION is generated if one of the unsigned integer variants of this function is used to
     load a uniform variable of type int, ivec2, ivec3, ivec4, or an array of these.
     
     GL_INVALID_OPERATION is generated if location is an invalid uniform location for the current program
     object and location is not equal to -1.
     
     GL_INVALID_OPERATION is generated if count is greater than 1 and the indicated uniform variable is
     not an array variable.
     
     GL_INVALID_OPERATION is generated if a sampler is loaded using a command other than glUniform1i and
     glUniform1iv.
     */
    
    // triangle vertex positions
    glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);//******************THIS IS AN ISSUE*********
    HandleGLError("Error7");
    glEnableVertexAttribArray(0);
    HandleGLError("Error8");
    glVertexAttribPointer(0,                  // attribute
                          3,                  // size
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          3*sizeof(glm::vec3),// stride
                          (void*)0            // array buffer offset
                          );
    HandleGLError("Error9");
    // triangle vertex normals
    glEnableVertexAttribArray(1);
    HandleGLError("Error10");
    glVertexAttribPointer(1,                      // attribute
                          3,                      // size
                          GL_FLOAT,               // type
                          GL_FALSE,               // normalized?
                          3*sizeof(glm::vec3),    // stride
                          (void*)sizeof(glm::vec3)// array buffer offset
                          );
    HandleGLError("Error11");
    // triangle vertex colors
    glEnableVertexAttribArray(2);
    HandleGLError("Error12");
    glVertexAttribPointer(2,                          // attribute
                          3,                          // size
                          GL_FLOAT,                   // type
                          GL_FALSE,                   // normalized?
                          3*sizeof(glm::vec3),        // stride
                          (void*)(sizeof(glm::vec3)*2)// array buffer offset
                          );
    HandleGLError("Error13");
    /*glEnableVertexAttribArray
     GL_INVALID_OPERATION is generated if size is GL_BGRA and type is not GL_UNSIGNED_BYTE, GL_INT_2_10_10_10_REV or GL_UNSIGNED_INT_2_10_10_10_REV.
     
     GL_INVALID_OPERATION is generated if type is GL_INT_2_10_10_10_REV or GL_UNSIGNED_INT_2_10_10_10_REV and size is not 4 or GL_BGRA.
     
     GL_INVALID_OPERATION is generated if type is GL_UNSIGNED_INT_10F_11F_11F_REV and size is not 3.
     
     GL_INVALID_OPERATION is generated by glVertexAttribPointer if size is GL_BGRA and noramlized is GL_FALSE.
     
     GL_INVALID_OPERATION is generated if zero is bound to the GL_ARRAY_BUFFER buffer object binding point and the pointer argument is not NULL.
     */
    
    
    // triangle indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);//************THIS IS AN ISSUE*********
    HandleGLError("Error14");
    glDrawElements(GL_TRIANGLES,         // mode
                   num_mini_triangles*3, // count
                   GL_UNSIGNED_INT,      // type
                   (void*)0              // element array buffer offset
                   );
    HandleGLError("Error15");
    /*
     GL_INVALID_OPERATION is generated if a non-zero buffer object name is bound to an enabled array or the element array and the buffer object's data store is currently mapped.
     
     GL_INVALID_OPERATION is generated if glDrawElements is executed between the execution of glBegin and the corresponding glEnd.
     */
    glDisableVertexAttribArray(0);
    HandleGLError("Error16");
    glDisableVertexAttribArray(1);
    HandleGLError("Error17");
    glDisableVertexAttribArray(2);
    HandleGLError("Error18");
    
    // =================================
    // draw the different types of edges
    //  if (args->wireframe) {
    
    HandleGLError("leaving drawVBOs");
}

void Mesh::cleanupVBOs() {
    glDeleteBuffers(1, &mesh_VAO);
    glDeleteBuffers(1, &mesh_tri_verts_VBO);
    glDeleteBuffers(1, &mesh_tri_indices_VBO);
}

void Mesh::setupVBOs()
{
 HandleGLError("enter setupVBOs");
 
 std::vector<glm::vec3> indexed_verts;
 std::vector<unsigned int> mesh_tri_indices;
 
 //for keeping track of verts we have already calculated normals for
 //std::map<int, bool> vertsNormalized;
 
 //int i = 0;
 // write the vertex & triangle data
 ////printf("NumTriangles: %d\n", (int)triangles.size());
 for (triangleshashtype::iterator iter = triangles.begin();
 iter != triangles.end(); iter++)
 {
 
 Triangle *tri = iter->second;
 Edge *e = tri->getEdge();
 // grab the vertex positions
 //glm::vec3 a = (*tri)[0]->get();
 //glm::vec3 b = (*tri)[1]->get();
 //glm::vec3 c = (*tri)[2]->get();
 
 // determine edge colors (when wireframe is enabled)
 //glm::vec3 edgecolor_ab = EdgeColor(tri->getEdge());
 //glm::vec3 edgecolor_bc = EdgeColor(tri->getEdge()->getNext());
 //glm::vec3 edgecolor_ca = EdgeColor(tri->getEdge()->getNext()->getNext());
 
 
 //STUFF FROM RADIOSITY**********
 double avg_s = 0;
 double avg_t = 0;
 glm::vec3 avg_color(0,0,0);
 
 int start = mesh_tri_verts.size(); // ***********THIS NEEDS TO BE INITIALIZED******
 
 // wireframe is normally black, except when it's the special
 // patch, then the wireframe is red
 glm::vec4 wireframe_color(0,0,0,0.5);
 //if (args->render_mode == RENDER_FORM_FACTORS && i == max_undistributed_patch) {
 //    wireframe_color = glm::vec4(1,0,0,1);
 //}
 
 glm::vec3 normal = tri->computeNormal();
 
 // add the 3 corner vertices
 for (int j = 0; j < 3; j++)
 {
 glm::vec3 pos = ((*tri)[j])->get();
 double s = (*tri)[j]->get_s();
 double t = (*tri)[j]->get_t();
 //glm::vec3 color = setupHelperForColor(tri,i,j);
 //color = glm::vec3(linear_to_srgb(color.r),
 //                  linear_to_srgb(color.g),
 //                  linear_to_srgb(color.b));
 glm::vec3 color = iter->second->getMaterial()->getDiffuseColor();
 //avg_color += 0.25f * color;
 mesh_tri_verts.push_back(VBOPosNormalColor(pos,normal,
 glm::vec4(color.r,color.g,color.b,1.0),
 wireframe_color,
 s,t));
 avg_s += 0.33333 * s;
 avg_t += 0.33333 * t;
 e = e->getNext();
 }
 
 glm::vec3 centroid = tri->computeCentroid();
 mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,
 glm::vec4(avg_color.r,avg_color.g,avg_color.b,1),
 glm::vec4(1,1,1,1),
 avg_s,avg_t));
 
 if (tri->getMaterial()->hasTextureMap())
 {
 mesh_textured_tri_indices.push_back(VBOIndexedTri(start+0,start+1,start+4));
 mesh_textured_tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+4));
 mesh_textured_tri_indices.push_back(VBOIndexedTri(start+2,start+3,start+4));
 mesh_textured_tri_indices.push_back(VBOIndexedTri(start+3,start+0,start+4));
 }
 else
 {
 mesh_tri_indices_vec.push_back(VBOIndexedTri(start+0,start+1,start+4));
 mesh_tri_indices_vec.push_back(VBOIndexedTri(start+1,start+2,start+4));
 mesh_tri_indices_vec.push_back(VBOIndexedTri(start+2,start+3,start+4));
 mesh_tri_indices_vec.push_back(VBOIndexedTri(start+3,start+0,start+4));
 }
 if (args->gouraud)
 {
 // =====================================
 // ASSIGNMENT: complete this functionality
 // =====================================
 //Runs at O(F) where F is the number of faces.
 
 //I want to get rid of the redundant computations
 std::map<int, bool> vertsProccessed;
 
 //go through each vertex of the current triangle and calculate vert norms from them using the adjacent faces of each
 for (int i = 0; i<3; i++)
 {
 //to keep track of where we are
 Edge::Edge* ptr = t->getEdge();
 
 //(*t)[i] target vertex, lets ensure we are there.
 while (ptr->getStartVertex() != (*t)[i])
 {
 ptr = ptr->getNext();
 }
 
 //averages all the face normals adjacent to the current edges startVertex and stores it in the vertex
 (*t)[i]->setNorm(calculateVertNormal(ptr));
 
 TriVBOHelper(indexed_verts,mesh_tri_indices,
 a,b,c,
 (*t)[0]->getNorm(),(*t)[1]->getNorm(),(*t)[2]->getNorm(),
 //(*t)[0]->getNorm(),(*t)[1]->getNorm(),(*t)[0]->getNorm());
 //a, b, c);
 edgecolor_ab,edgecolor_bc,edgecolor_ca);
 
 }
 }
 else
 {
 // for flat shading, use the triangle normal at each vertex
 // use the normal of the triangl
 
 
 //glm::vec3 normal = ComputeNormal(a,b,c);
 //TriVBOHelper(indexed_verts,mesh_tri_indices,
 //                 a,b,c,
 //                 normal,normal,normal,
 //                 edgecolor_ab,edgecolor_bc,edgecolor_ca);
 
 //}
 //i++;
 }
 
 // the vertex data
 glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);
 glBufferData(GL_ARRAY_BUFFER, indexed_verts.size() * sizeof(glm::vec3), &indexed_verts[0], GL_STATIC_DRAW);
 // the index data (refers to vertex data)
 glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);
 glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices.size() * sizeof(unsigned int), &mesh_tri_indices[0] , GL_STATIC_DRAW);
 
 num_mini_triangles = mesh_tri_indices.size();
 
 HandleGLError("leaving setupVBOs");
 }
}
void Mesh::initializeVBOs() {
        HandleGLError("enter initialize VBOs");
        
        // create a pointer for the vertex & index VBOs
        glGenVertexArrays(1, &mesh_VAO);
        glBindVertexArray(mesh_VAO);
        glGenBuffers(1, &mesh_tri_verts_VBO);
        glGenBuffers(1, &mesh_tri_indices_VBO);
        // and the data to pass to the shaders
        GLCanvas::MatrixID = glGetUniformLocation(GLCanvas::programID, "MVP");
        GLCanvas::LightID = glGetUniformLocation(GLCanvas::programID, "LightPosition_worldspace");
        GLCanvas::ViewMatrixID = glGetUniformLocation(GLCanvas::programID, "V");
        GLCanvas::ModelMatrixID = glGetUniformLocation(GLCanvas::programID, "M");
        
        GLCanvas::wireframeID = glGetUniformLocation(GLCanvas::programID, "wireframe");
        
        // call this the first time...
        setupVBOs();
        HandleGLError("leaving initializeVBOs");
    }

#endif

