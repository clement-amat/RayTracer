
#include "raytracer.h"
#include "scene_types.h"
#include "ray.h"
#include "image.h"
#include "kdtree.h"
#include <stdio.h>
#include <cmath>


/// acne_eps is a small constant used to prevent acne when computing intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;
const float  INV_PI_F= 0.318309f;

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
      intersection->mat = &plane->mat;
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
      intersection->mat = &obj->mat;
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
        intersection->mat = &obj->mat;
        intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
        return true;

      } else if (t1Valide) {    // seul t1 valide
        intersection->position = ray->orig + t1 * ray->dir;
        intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
        ray->tmax = t1;
        intersection->mat = &obj->mat;
        return true;

      } else if (t2Valide) {    // seul t2 valide
        intersection->position = ray->orig + t2 * ray->dir;
        intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
        ray->tmax = t2;
        intersection->mat = &obj->mat;
        return true;
      }
      
  }
  
  // pas de pts d'intersection valide
  return false;

}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  size_t objectCount = scene->objects.size();
  bool hasIntersection = false;

  // loop on each object of the scene to compute intersection
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
  float NdotH2 = NdotH * NdotH;
  float alpha2 = alpha * alpha;
  float tan2NDotH = (1.f - NdotH2) / NdotH2;
  return exp(-tan2NDotH / alpha2) / (M_PI * alpha2 * NdotH2 * NdotH2);
}

// Fresnel term computation. Implantation of the exact computation. we can use the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {
  float rs, rp;
  float cosTetaI = LdotH;
  float n1 = extIOR;
  float n2 = intIOR;

  float sin2TetaT = pow(n1/n2, 2) * (1 - cosTetaI * cosTetaI);
  if (sin2TetaT > 1.f) {
    return 1;
  }

  float cosTetaT = sqrt(1 - sin2TetaT);
  rs = pow(n1 * cosTetaI - n2 * cosTetaT, 2) / pow( n1 * cosTetaI + n2 * cosTetaT , 2);
  rp = pow( n1 * cosTetaT - n2 * cosTetaI , 2) / pow( n1 * cosTetaT + n2 * cosTetaI, 2);

  return 0.5 * (rs + rp);
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

/**
 * n : normale à la surface
 * v : direction de l'observation
 * l : direction de l'éclairage
 * lc : couleur de la lumière
 * mat : materiau
 */ 
color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat ){
  color3 ret = color3(0);

  float scalaire = dot(l, n);
  if (scalaire > 0.f) {
    ret = color3(mat->diffuseColor * INV_PI_F * scalaire * lc);
  } 

  return ret;
	    
}

//! if tree is not null, use intersectKdTree to compute the intersection instead of intersect scene
color3 trace_ray(Scene * scene, Ray *ray, KdTree *tree) {  
  Intersection intersection;
  Intersection intersectionOmbre;

  if (intersectScene(scene, ray, &intersection)) {   // intersection entre le rayon est un objet de la scene
    color3 color = color3(0);
    Ray rayonOmbre;

    // pour chaque source de lumière additionner les contributions :
    for (size_t i = 0 ; i < scene->lights.size() ; i++) {
      vec3 l = normalize(scene->lights.at(i)->position-intersection.position);
      rayInit(&rayonOmbre, intersection.position + acne_eps * l, l);
      
      // si l'objet n'est pas caché par une ombre
      if (!intersectScene(scene, &rayonOmbre, &intersectionOmbre)){
        color += shade(-intersection.normal, - ray->dir, -l,
                            scene->lights.at(i)->color, intersection.mat);
      }
    }
    
    return color;
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
