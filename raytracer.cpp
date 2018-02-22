
#include "raytracer.h"
#include "scene_types.h"
#include "ray.h"
#include "image.h"
#include "kdtree.h"
#include <stdio.h>


/// acne_eps is a small constant used to prevent acne when computing intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *plane) {
  vec3 n = plane->geom.plane.normal;
  vec3 d = ray->dir;
  float dist = plane->geom.plane.dist;
  float dotDN = dot(d, n);
  
  if (abs(dotDN) > acne_eps) {
    point3 o = ray->orig;
    float t;

    t = -(dot(o,n) + dist) / dotDN;
    if (t >= ray->tmin && t <= ray->tmax) {
      ray->tmax = t;
      // normal intersection = normal au plan
      // pt intersect : o + t* direction
      intersection->normal = n ;
      intersection->position = o + t * d;
      return true;
    }
  }
  return false;

}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {
  float a = 1;
  float b = 2 * dot(ray->dir, (ray->orig - obj->geom.sphere.center));

  vec3 oMinusC = (ray->orig - obj->geom.sphere.center);
  float c = dot(oMinusC,oMinusC) 
              - obj->geom.sphere.radius * obj->geom.sphere.radius;

  float delta = b * b - 4.f * a * c;

  if (delta == 0.f) {     // une solution
    float t = (- b) / 2.f * a;

    if (t >= ray->tmin && t <= ray->tmax) {
      ray->tmax = t;
      intersection->position = ray->orig + t * ray->dir;
      intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
      return true;
    }

  } else if (delta > 0.f) {    // deux solutions 
      float racineDelta = sqrt(delta);
      float t1 = (- b - racineDelta) / 2.f * a;
      float t2 = (- b + racineDelta) / 2.f * a;
      bool t1Valide = t1 >= ray->tmin && t1 <= ray->tmax;
      bool t2Valide = t2 >= ray->tmin && t2 <= ray->tmax;
      
      if (t1Valide && t2Valide) {
        if (t1 < t2) {
          intersection->position = ray->orig + t1 * ray->dir;
          ray->tmax = t1;
        } else {
          intersection->position = ray->orig + t2 * ray->dir;
          ray->tmax = t2;
        }
        intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
        return true;

      } else if (t1Valide) {    // seul t1 valide
        intersection->position = ray->orig + t1 * ray->dir;
        intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
        ray->tmax = t1;
        return true;

      } else if (t2Valide) {    // seul t2 valide
        intersection->position = ray->orig + t2 * ray->dir;
        intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
        ray->tmax = t2;
        return true;
      }
  }
  
  // pas de pts d'intersection valide
  return false;

}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  size_t objectCount = scene->objects.size();
  bool hasIntersection = false;

  //!\todo loop on each object of the scene to compute intersection
  for (size_t i = 0 ; i < objectCount ; i++) {
    if (scene->objects[i]->geom.type == SPHERE) {
      if (intersectSphere(ray, intersection, scene->objects[i])) {
        hasIntersection = true;
      }
    } else if (scene->objects[i]->geom.type == PLANE){
      if (intersectPlane(ray, intersection, scene->objects[i])) {
        hasIntersection = true;
      }
    }
  }

  return hasIntersection;
}

/* -------------------- ------------------------------------------------------- */
/*
 *	The following functions are coded from Cook-Torrance bsdf model description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba renderer)
 */

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha) {


  //! \todo compute Beckmann normal distribution
  return 0.5f;

}

// Fresnel term computation. Implantation of the exact computation. we can use the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {

  //! \todo compute Fresnel term
  return 0.5f;

}


// Shadowing and masking function. Linked with the NDF. Here, Smith function, suitable for Beckmann NDF
float RDM_chiplus(float c) {
  return (c > 0.f) ? 1.f : 0.f;
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {

  //!\todo compute G1 term of the Smith fonction
  return 0.5f;

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN, float alpha) {

  //!\todo the Smith fonction
  return 0.5f;


}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {

  //!\todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G = RDM_Smith
  return color3(.5f);

  
}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {

  //!\todo compute diffuse component of the bsdf
  return color3(.5f);

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {

  //! \todo compute bsdf diffuse and specular term
  return color3(0.f);

}




/* --------------------------------------------------------------------------- */

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat ){
  color3 ret = color3(0.f);

  //! \todo compute bsdf, return the shaded color taking into account the
  //! lightcolor
  

  return ret;
	    
}

//! if tree is not null, use intersectKdTree to compute the intersection instead of intersect scene
color3 trace_ray(Scene * scene, Ray *ray, KdTree *tree) {  
  Intersection intersection;

  if (intersectScene(scene, ray, &intersection)) {
    return color3(0.5f * intersection.normal + 0.5f);
  } else {
    return scene->skyColor;
  }
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing and kdtree initializaion
  float aspect = 1.f/scene->cam.aspect;
    
  KdTree *tree =  NULL;


  //! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f); //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step 
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) * aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x = (0.5f - img->width * 0.5f) / (img->width * 0.5f) *scene->cam.xdir;
  
    
  for(size_t j=0; j<img->height; j++) {
    if(j!=0) printf("\033[A\r");
    float progress = (float)j/img->height*100.f;
    printf("progress\t[");
    int cpt = 0;
    for(cpt = 0; cpt<progress; cpt+=5) printf(".");
    for(       ; cpt<100; cpt+=5) printf(" ");
    printf("]\n");
#pragma omp parallel for
    for(size_t i=0; i<img->width; i++) {
      color3 *ptr = getPixelPtr(img, i,j);
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y + float(i)*dx + float(j)*dy;

      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      *ptr = trace_ray(scene, &rx, tree);

    }
  }
}
