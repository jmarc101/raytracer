/*
JeanMarc Prudhomme, MAT 20137035
Jean-Daniel Toupin


*/

#include "object.hpp"

#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>  

const double EPSILON = 0.00001;



bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // Assure une valeur correcte pour la coordonnée W de l'origine et de la direction
	// Vous pouvez commentez ces lignes si vous faites très attention à la façon de construire vos rayons.
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    Ray local_ray(i_transform * ray.origin, i_transform * ray.direction);
	//!!! NOTE UTILE : pour calculer la profondeur dans localIntersect(), si l'intersection se passe à
	// ray.origin + ray.direction * t, alors t est la profondeur
	//!!! NOTE UTILE : ici, la direction peut êytre mise à l'échelle, alors vous devez la renormaliser
	// dans localIntersect(), ou vous aurez une profondeur dans le système de coordonnées local, qui
	// ne pourra pas être comparée aux intersection avec les autres objets.
    if (localIntersect(local_ray, hit)) 
	{
        // Assure la valeur correcte de W.
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
		// Transforme les coordonnées de l'intersection dans le repère global.
        hit.position = transform * hit.position;
        hit.normal = (n_transform * hit.normal).normalized();
        
		return true;
    }

    return false;
}


bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const 
{
    // @@@@@@ VOTRE CODE ICI
	// Vous pourriez aussi utiliser des relations géométriques pures plutôt que les
	// outils analytiques présentés dans les slides.
	// Ici, dans le système de coordonées local, la sphère est centrée en (0, 0, 0)
	// avec un radius de 1.0.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.

	//Algo inspire par http://viclw17.github.io/2018/07/16/raytracing-ray-sphere-intersection/


	/*
	Sphere: dot((P−C),(P−C))=r2
	Ray: p(t)=A+tB
	Combined: dot((A+tB−C),(A+tB−C))=r2
	
	a =  || rayDirection ||
	b = 2 * cross(B,A−C)
	c = cross(A−C,A−C)−r2
	*/

	
	double r = this->radius;

	double a = ray.direction.length2();
	double b = 2 * ray.origin.dot(ray.direction);
	double c = ray.origin.length2() - pow(r, 2);


	

	//equation quadratic pour calculer notre t
	double discriminant = (pow(b, 2)) - (4 * a * c);
	double t;
	if (discriminant < 0) { // no solution
		return false;
	}
	else if (discriminant > 0) { //two solutions
		double t1 = (-b + sqrt(discriminant)) / (2 * a);
		double t2 = (-b - sqrt(discriminant)) / (2 * a);
		t = fmax(0, fmin(t1, t2));
		if (t <= 0){
			return false;
		}
	} else if (discriminant == 0) {
		double t1 = (-b) / (2 * a);
		if (t1 < 0){
			return false;
		}
		t = t1;
	}
	
	if (hit.depth > t) {

		hit.depth = t;
		Vector rayPosition = ray.origin + ray.direction  * t;
		hit.position = rayPosition;
		hit.normal = rayPosition.normalized();
		return true;
	}
	
	return false;
	
	
	
}



bool Plane::localIntersect(Ray const &ray, Intersection &hit) const
{
	
	// Nous avons que l'equiation de la plane est z = 0 
	// Donc z0 + tzd = 0 -> -z0/zd = t
	double t = -ray.origin[2] / ray.direction[2];
	
	if (t < 0 || hit.depth <= t){
		return false;
	}


	hit.depth = t;
	Vector positionRay = ray.origin + ray.direction * t;
	hit.position = positionRay;
	hit.normal = Vector(0,0,1);

	return true;

}

// Intersections !
bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
	// Test de la boite englobante
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++) {
		if (ray.direction[i] == 0.0) {
			if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) {
				// Rayon parallèle à un plan de la boite englobante et en dehors de la boite
				return false;
			}
			// Rayon parallèle à un plan de la boite et dans la boite: on continue
		}
		else {
			double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Assure t1 <= t2

			if (t1 > tNear) tNear = t1; // On veut le plus lointain tNear.
			if (t2 < tFar) tFar = t2; // On veut le plus proche tFar.

			if (tNear > tFar) return false; // Le rayon rate la boite englobante.
			if (tFar < 0) return false; // La boite englobante est derrière le rayon.
		}
	}
	// Si on arrive jusqu'ici, c'est que le rayon a intersecté la boite englobante.

	// Le rayon interesecte la boite englobante, donc on teste chaque triangle.
	bool isHit = false;
	for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
		Triangle const &tri = triangles[tri_i];

		if (intersectTriangle(ray, tri, hit)) {
			isHit = true;
		}
	}
	return isHit;
}

double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}

bool Mesh::intersectTriangle(Ray const &ray,
	Triangle const &tri,
	Intersection &hit) const
{
	// Extrait chaque position de sommet des données du maillage.
	Vector const &p0 = positions[tri[0].pi];
	Vector const &p1 = positions[tri[1].pi];
	Vector const &p2 = positions[tri[2].pi];

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Vous pourriez trouver utile d'utiliser la routine implicitLineEquation()
	// pour calculer le résultat de l'équation de ligne implicite en 2D.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.
	//!!! NOTE UTILE : pour le point d'intersection, sa normale doit satisfaire hit.normal.dot(ray.direction) < 0


	// Inpirer de https://courses.cs.washington.edu/courses/csep557/10au/lectures/triangle_intersection.pdf
	// ++ slides

	Vector edge0 = p1 - p0;
	Vector edge1 = p2 - p1;
	Vector edge2 = p0 - p2;
	Vector normal = edge1.cross(edge2);
	normal.normalize();


		// n * x(p0) = 0
	float d = normal.dot(p0);

	// t= (d - n*p)/n*d 
	// Si parallele ne se croise pas /
	float nDotD = normal.dot(ray.direction);
	float t = (d - normal.dot(ray.origin) ) / nDotD;


	if (fabs(nDotD) < EPSILON) {return false;} 

	if (t < 0) {return false;}

	Vector pointQ = ray.origin + t * ray.direction;

	double bMinusA = edge0.cross(pointQ - p0).dot(normal);
	double cMinusB = edge1.cross(pointQ - p1).dot(normal);
	double aMinusC = edge2.cross(pointQ - p2).dot(normal);

	if (bMinusA < 0 || cMinusB < 0 || aMinusC < 0){
		return false;
	}

	if (hit.depth > t){
		hit.depth = t;
		hit.position = pointQ;
		hit.normal =  normal;
		return true;
	}
	return false;
}


bool Plane::is_Plane() const {return true;}
bool Sphere::is_Plane() const { return false;}
bool Mesh::is_Plane() const { return false;}

