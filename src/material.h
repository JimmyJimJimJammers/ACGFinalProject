#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "glCanvas.h"

#include <cassert>
#include <string>

#include "image.h"

class ArgParser;
class Ray;
class Hit;

// ====================================================================
// ====================================================================
// A simple Phong-like material 

class Material {

public:

  Material(const std::string &texture_file, const glm::vec3 &d_color,
	   const glm::vec3 &r_color, const glm::vec3 &e_color, float roughness_) {
    textureFile = texture_file;
    if (textureFile != "") {
      image = new Image(textureFile);
      ComputeAverageTextureColor();
    } else {
      diffuseColor = d_color;
      image = NULL;
    }
    reflectiveColor = r_color;
    emittedColor = e_color;
    roughness = roughness_;
    // need to initialize texture_id after glut has started
    texture_id = 0;
  }
    
    
    //                          texture_file, displace_file, diffuse, displace, reflective,emitted,roughness
    Material(const std::string &texture_file, const std::string &disp_file, const glm::vec3 &d_color, const glm::vec3 &disp_color,
             const glm::vec3 &r_color, const glm::vec3 &e_color, float roughness_) {
        textureFile = texture_file;
        if (textureFile != "") {
            image = new Image(textureFile);
            ComputeAverageTextureColor();
        } else {
            diffuseColor = d_color;
            image = NULL;
        }
        //printf("Disp_File: %s\n", disp_file.c_str());
        if (disp_file != "")
        {
            displacementImage = new Image(disp_file);
            displace_file = disp_file;
        }
        else
        {
            displacementColor = disp_color;
            displacementImage = NULL;
        }
        
        reflectiveColor = r_color;
        emittedColor = e_color;
        roughness = roughness_;
        // need to initialize texture_id after glut has started
        texture_id = 0;
    }
  
  ~Material();

  // ACCESSORS
    //ADDED FOR DISPLACEMENT
    const glm::vec3 getDisplacementValue(float s, float t) const;
    const glm::vec3 &getDisplacementValue() const { return displacementColor;}
    
  const glm::vec3& getDiffuseColor() const { return diffuseColor; }
  const glm::vec3 getDiffuseColor(float s, float t) const;
  const glm::vec3& getReflectiveColor() const { return reflectiveColor; }
  const glm::vec3& getEmittedColor() const { return emittedColor; }  
  float getRoughness() const { return roughness; } 
  bool hasTextureMap() const { return (textureFile != ""); }
    
    //ADDED FOR DISPLACEMENT
    bool hasDisplacementMap() const
    {
        //printf("DisplaceFile: %s\n", displace_file.c_str());
        return (displace_file != "");
    }
    
  GLuint getTextureID();

  // SHADE
  // compute the contribution to local illumination at this point for
  // a particular light source
  glm::vec3 Shade(const Ray &ray, const Hit &hit, const glm::vec3 &dirToLight, const glm::vec3 &lightColor, ArgParser *args) const;
  
protected:

  Material() { exit(0); }
  Material(const Material&) { exit(0); }
  const Material& operator=(const Material&) { exit(0); }

  void ComputeAverageTextureColor();

  // REPRESENTATION
    //ADDED FOR DISPLACEMENT
    glm::vec3 displacementColor;
    
  glm::vec3 diffuseColor;
  glm::vec3 reflectiveColor;
  glm::vec3 emittedColor;
  float roughness;
    //ADDED FOR DISPLACEMENT
    std::string displace_file;
  std::string textureFile;
  GLuint texture_id;
  Image *image;
    
    //ADDED FOR DISPLACEMENT (MIGHT NOT BE NECESSARY)
    Image *displacementImage;
    
    
    //To Do:
    //Look into where material is set
    //make displacement map includable from command line
    //see what Image is
    //look into texturefile and how it's implemented
    //rewrite const glm::vec3 getDiDisplacementValue(float s, float t) const;
    //  so that it actually works correctly
};

// ====================================================================
// ====================================================================

#endif

