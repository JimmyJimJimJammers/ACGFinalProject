#ifndef MESH_H
#define MESH_H

#include "glCanvas.h"

#include <vector>
#include "hash.h"
#include "material.h"
#include "triangle.h"

class Vertex;
class Edge;
class BoundingBox;
class Face;
class Primitive;
class ArgParser;
class Ray;
class Hit;
class Camera;

enum FACE_TYPE {FACE_TYPE_ORIGINAL, FACE_TYPE_RASTERIZED, FACE_TYPE_SUBDIVIDED, FACE_TYPE_DISPLACEMENT};
enum TRI_TYPE {TRI_TYPE_ORIGINAL, TRI_TYPE_SUBDIVIDED, TRI_TYPE_DISPLACEMENT};

// ======================================================================
// ======================================================================
// A class to store all objects in the scene.  The quad faces of the
// mesh can be subdivided to improve the resolution of the radiosity
// solution.  The original mesh is maintained for efficient occlusion
// testing.

class Mesh {

public:

  // ===============================
  // CONSTRUCTOR & DESTRUCTOR & LOAD
  Mesh() { bbox = NULL; }
  virtual ~Mesh();
  void Load(ArgParser *_args);
    
  // ========
  // VERTICES
  int numVertices() const { return vertices.size(); }
  Vertex* addVertex(const glm::vec3 &pos);
  // look up vertex by index from original .obj file
  Vertex* getVertex(int i) const {
    assert (i >= 0 && i < numVertices());
    return vertices[i]; }
  // this creates a relationship between 3 vertices (2 parents, 1 child)
  void setParentsChild(Vertex *p1, Vertex *p2, Vertex *child);
  // this accessor will find a child vertex (if it exists) when given
  // two parent vertices
  Vertex* getChildVertex(Vertex *p1, Vertex *p2) const;

  // =====
  // EDGES
  int numEdges() const { return edges.size(); }
  // this efficiently looks for an edge with the given vertices, using a hash table
  Edge* getEdge(Vertex *a, Vertex *b) const;
  const edgeshashtype& getEdges() const { return edges; }

  // =================
  // ACCESS THE LIGHTS
  std::vector<Face*>& getLights() { return original_lights; }

  // ==================================
  // ACCESS THE QUADS (for ray tracing)
  int numOriginalQuads() const { return original_quads.size(); }
  Face* getOriginalQuad(int i) const {
    assert (i < numOriginalQuads());
    return original_quads[i]; }
    
    int numDispQuads() const { return disp_subdivided_quads.size(); }
    
    Face* getDispQuad(int i) const {
        assert (i < numDispQuads());
        return disp_subdivided_quads[i]; }

    
  // =======================================
  // ACCESS THE PRIMITIVES (for ray tracing)
  int numPrimitives() const { return primitives.size(); }
  Primitive* getPrimitive(int i) const {
    assert (i >= 0 && i < numPrimitives()); 
    return primitives[i]; }
  // ACCESS THE PRIMITIVES (for radiosity)
  int numRasterizedPrimitiveFaces() const { return rasterized_primitive_faces.size(); }
  Face* getRasterizedPrimitiveFace(int i) const {
    assert (i >= 0 && i < numRasterizedPrimitiveFaces());
    return rasterized_primitive_faces[i]; }

  // ==============================================================
  // ACCESS THE SUBDIVIDED QUADS + RASTERIZED FACES (for radiosity)
  int numFaces() const { return subdivided_quads.size() + rasterized_primitive_faces.size(); }//possibly need disp_sub_quads

    int numTris() const { return triangles.size(); }//possibly need disp_sub_quads
    
    
  Face* getFace(int i) const {
    int num_faces = numFaces();
    assert (i >= 0 && i < num_faces);
    if (i < (int)subdivided_quads.size()) return subdivided_quads[i];
    else return getRasterizedPrimitiveFace(i-subdivided_quads.size()); }
    
    Triangle* getTri(int i)
    {
        return triangles[i];
    }

  // ============================
  // CREATE OR SUBDIVIDE GEOMETRY
  void addRasterizedPrimitiveFace(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material) {
    addFace(a,b,c,d,material,FACE_TYPE_RASTERIZED); }
  void addOriginalQuad(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material) {
    addFace(a,b,c,d,material,FACE_TYPE_ORIGINAL); }
  void addSubdividedQuad(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material) {
    addFace(a,b,c,d,material,FACE_TYPE_SUBDIVIDED); }
    //DISPLACEMENT ADD
    void addDisplacementQuad(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material)
    {
        addFace(a,b,c,d,material,FACE_TYPE_DISPLACEMENT);
    }

  // ===============
  // OTHER ACCESSORS
  BoundingBox* getBoundingBox() const { return bbox; }

  // ===============
  // OTHER FUNCTIONS
  void Subdivision();
    
    //displacement subdivision stuff
    void DisplacementSubdivision();
    void setMaxFaceDisps();
    float findLargestDispBetween2(Vertex *a, Vertex *b, Material *m);
    struct newVertReturn
    {
        Vertex* newVert;
        bool significant;
        float dispSize;
    };

    //TRIANGLE DISPLACEMENT STUFF
    void TriVBOHelper( std::vector<glm::vec3> &indexed_verts,
                            std::vector<unsigned int> &mesh_tri_indices,
                            const glm::vec3 &pos_a,
                            const glm::vec3 &pos_b,
                            const glm::vec3 &pos_c,
                            const glm::vec3 &normal_a,
                            const glm::vec3 &normal_b,
                            const glm::vec3 &normal_c,
                            const glm::vec3 &color_ab,
                            const glm::vec3 &color_bc,
                            const glm::vec3 &color_ca);
    void initializeVBOs();
    void setupVBOs();
    void drawVBOs(const glm::mat4 &ProjectionMatrix,const glm::mat4 &ViewMatrix,const glm::mat4 &ModelMatrix);
    void drawVBOs();
    void cleanupVBOs();
    void removeSubdivideTriangle(Triangle *t);
    Triangle* getOppositeTriangle(Triangle* currentTri, Vertex* a, Vertex* b);
    void DisplacementSubdivisionTriangles();
    Triangle* makeSubTriangle(Vertex *center, Vertex *a, Vertex* b, Triangle* mate, Triangle* tri);
    void splitR3(Triangle* tri);
    void swap(Edge* sharedEdge);
    int numDispTriangles()
    {
        //if (disp_subdivided_tris.size() > 0)
        //{
        //    return disp_subdivided_tris.size();
        //}
        //else
        //{
            return subdivided_tris.size();
        //}
    }

private:

  // ==================================================
  // HELPER FUNCTIONS FOR CREATING/SUBDIVIDING GEOMETRY
  Vertex* AddEdgeVertex(Vertex *a, Vertex *b);
  Vertex* AddMidVertex(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
  void addFace(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material, enum FACE_TYPE face_type);
  void removeFaceEdges(Face *f);
  void addPrimitive(Primitive *p);
    
    Triangle* addTri(Vertex *a, Vertex *b, Vertex *c, Material *material, enum TRI_TYPE tri_type);
    void removeTriEdges(Triangle *t);
    void ConvertOriginalQuadsToTris();
    
    //helper functions for displacement
    newVertReturn AddMidVertexDisp(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *m);
    newVertReturn AddEdgeVertexDisp(Vertex *a, Vertex *b, Material *m);
    //Vertex* AddEdgeVertexDisp(Vertex *a, Vertex *b, Material *m);
    //Vertex* AddMidVertexDisp(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *m);
    void SubdivisionAltered();
    void TriangleSimplification(int target_tri_count);
    void SubdivisionAlteredTriangles();

  // ==============
  // REPRESENTATION
  ArgParser *args;
 public:
  std::vector<Material*> materials;
  glm::vec3 background_color;
  Camera *camera;
 private:

  // the bounding box of all rasterized faces in the scene
  BoundingBox *bbox; 

  // the vertices & edges used by all quads (including rasterized primitives)
  std::vector<Vertex*> vertices;  
  edgeshashtype edges;
  vphashtype vertex_parents;

  // the quads from the .obj file (before subdivision)
  std::vector<Face*> original_quads;
  // the quads from the .obj file that have non-zero emission value
  std::vector<Face*> original_lights; 
  // all primitives (spheres, etc.)
  std::vector<Primitive*> primitives;
  // the primitives converted to quads
  std::vector<Face*> rasterized_primitive_faces;
  // the quads from the .obj file after subdivision
  std::vector<Face*> subdivided_quads;
    std::vector<Face*> disp_subdivided_quads;
    
    std::vector<Triangle*> original_tris;
    std::vector<Triangle*> subdivided_tris; // for tris
    std::vector<Triangle*> disp_subdivided_tris; // for tris

    //Stuff added from hw 1 (triangle related stuff)
    
    
    triangleshashtype triangles;
    
    // VBOs (Radiosity stuff)
    GLuint mesh_tri_verts_VBO;
    GLuint mesh_tri_indices_VBO;
    GLuint mesh_textured_tri_indices_VBO;
    
    std::vector<VBOPosNormalColor> mesh_tri_verts;
    std::vector<VBOIndexedTri> mesh_tri_indices;
    std::vector<VBOIndexedTri> mesh_textured_tri_indices;
    
    
    
    
    /*
    GLuint mesh_VAO;
    GLuint mesh_tri_verts_VBO;
    GLuint mesh_tri_indices_VBO;
    std::vector<VBOPosNormalColor> mesh_tri_verts;
    std::vector<VBOIndexedTri> mesh_tri_indices;
    std::vector<VBOIndexedTri> mesh_textured_tri_indices;
    int num_mini_triangles;*/
};

// ======================================================================
// ======================================================================

#endif




