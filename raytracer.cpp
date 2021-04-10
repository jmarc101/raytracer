/*
JeanMarc Prudhomme, MAT 20137035
Jean-Daniel Toupin

Ray tracer
super sampling ligne 94
shadowsampling ligne 326 ne fonctionne pas parfaitement, fait des pico donc mis le nombre a 1 pour l'instant
texture mapping  = 420
refraction + transparence + recursion a la fin

**/

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <random>

#include "bitmap_image.hpp"
#include "raytracer.hpp"
#include "image.hpp"

const double EPSILON = 0.00001;

void Raytracer::render(const char *filename, const char *depth_filename,
                       Scene const &scene)
{
    // Alloue les deux images qui seront sauvegardées à la fin du programme.
    Image colorImage(scene.resolution[0], scene.resolution[1]);
    Image depthImage(scene.resolution[0], scene.resolution[1]);
    
    // Crée le zBuffer.
    double *zBuffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        zBuffer[i] = DBL_MAX;
    }

	// @@@@@@ VOTRE CODE ICI
	// Calculez les paramètres de la caméra pour les rayons. Référez-vous aux slides pour les détails.
	//!!! NOTE UTILE : Le plan de vue peut être n'importe où, mais il sera implémenté différement.
	// Vous trouverez des références dans le cours.
	
	// Info de base
	// position de loeil
	Vector posInit = scene.camera.position;
	// Center point focus
	Vector focus = scene.camera.center;
	
	// Angle 
	double theta = scene.camera.fovy;
	double aspectRatio = scene.camera.aspect;
	// Near  / far
	double zNear = scene.camera.zNear;
	double zFar = scene.camera.zFar;


	// Vecteur resultant
	Vector upVect = scene.camera.up.normalized();;
	Vector forwardVect = (focus - posInit).normalized();
	Vector rightVect = forwardVect.cross(upVect);

	
	//!!! NOTE UTILE : tan() prend des radians plutot que des degrés.
	// Utilisez deg2rad() pour la conversion.
	double angle = tan(deg2rad(scene.camera.fovy / 2));
	double height = zNear * angle;
	double width = height * aspectRatio;


	// Calculer les longueurs entre le point central
	Vector centerVect = posInit + (zNear * forwardVect);
	Vector zeroVect = centerVect + (-width * rightVect) + (-height * upVect);
	double WidthVariation = (2* width) / scene.resolution[0];
	double HeigthVariation = (2* height) / scene.resolution[1];



	

	

    // Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {

        for(int x = 0; x < scene.resolution[0]; x++) {


			/*
			* ************************
			* SUPERSAMPLING RAYTRACER
			* ************************
			*/
			Vector color;
            double depth = scene.camera.zFar;

			// Chercher le parametre de super sampling, s'il y en a aucune alors on mets 1 par default
			int sampling = scene.supersampling <= 0 ? 1 : scene.supersampling ;
			



			for (int i = 0; i < sampling; i++)
			{

			// random for x,y si sampleNum = 1 alors on fait le raytracer de default sinon aleatoire
			double random = sampling == 1 ? 0.5 : (float)rand() / RAND_MAX;

            // Génère le rayon approprié pour ce pixel.
			Ray ray;

			if (scene.objects.empty())
			{
				// Pas d'objet dans la scène --> on rend la scène par défaut.
				// Pour celle-ci, le plan de vue est à z = 640 avec une largeur et une hauteur toute deux à 640 pixels.
				ray = Ray(scene.camera.position, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - scene.camera.position).normalized());
			}
			else
			{
				
				// Calculation des points pour creer le ray --> Random car nous faisons du super sampling
				Vector xFactor = (x + 1 * random) * WidthVariation * rightVect;
				Vector yFactor = (y + 1 * random) * HeigthVariation * upVect;
				Vector planePosition = zeroVect + xFactor + yFactor;
				Vector directionRay = planePosition - posInit;
				directionRay.normalize();

				//le ray resultant
				ray = Ray(posInit, directionRay);

			}

            // Initialise la profondeur de récursivité du rayon.
            int rayDepth = 0;
           
            // Notre lancer de rayons récursif calculera la couleur et la z-profondeur.
            Vector colorSampler;
            double depthSampler = scene.camera.zFar;

            // Ceci devrait être la profondeur maximum, correspondant à l'arrière plan.
            // NOTE : Ceci supVectpose que la direction du rayon est de longueur unitaire (normalisée)
			//        et que l'origine du rayon est à la position de la caméra.

            // Calcule la valeur du pixel en lançant le rayon dans la scène.
            trace(ray, rayDepth, scene, colorSampler, depthSampler);
			depth = min(depth, depthSampler);

			//additionne les couleurs ensemble de chaque ray
			color += colorSampler;
			}

			//divise par le nombre de ray cree
			color = color / sampling;

            // Test de profondeur
            if(depth >= scene.camera.zNear && depth <= scene.camera.zFar && 
                depth < zBuffer[x + y*scene.resolution[0]]) {
                zBuffer[x + y*scene.resolution[0]] = depth;

                // Met à jour la couleur de l'image (et sa profondeur)
                colorImage.setPixel(x, y, color);
                depthImage.setPixel(x, y, (depth-scene.camera.zNear) / 
                                        (scene.camera.zFar-scene.camera.zNear));
            }
        }

		// Affiche les informations de l'étape
		if (y % 100 == 0)
		{
			printf("Row %d pixels finished.\n", y);
		}
    }

	// Sauvegarde l'image
    colorImage.writeBMP(filename);
    depthImage.writeBMP(depth_filename);

	printf("Ray tracing finished with images saved.\n");

    delete[] zBuffer;
}


bool Raytracer::trace(Ray const &ray, 
                 int &rayDepth,
                 Scene const &scene,
                 Vector &outColor, double &depth)
{
    // Incrémente la profondeur du rayon.
    rayDepth++;

    // - itérer sur tous les objets en appelant calling Object::intersect.
    // - ne pas accepter les intersections plus lointaines que la profondeur donnée.
    // - appeler Raytracer::shade avec l'intersection la plus proche.
    // - renvoyer true ssi le rayon intersecte un objet.
	if (scene.objects.empty())
	{
		// Pas d'objet dans la scène --> on rend la scène par défaut :
		// Par défaut, un cube est centré en (0, 0, 1280 + 160) avec une longueur de côté de 320, juste en face de la caméra.
		// Test d'intersection :
		double x = 1280 / ray.direction[2] * ray.direction[0] + ray.origin[0];
		double y = 1280 / ray.direction[2] * ray.direction[1] + ray.origin[1];
		if ((x <= 160) && (x >= -160) && (y <= 160) && (y >= -160))
		{
			// S'il y a intersection :
			Material m; m.emission = Vector(16.0, 0, 0); m.reflect = 0; // seulement pour le matériau par défaut ; vous devrez utiliser le matériau de l'objet intersecté
			Intersection intersection;	// seulement par défaut ; vous devrez passer l'intersection trouvée par l'appel à Object::intersect()
			outColor = shade(ray, rayDepth, intersection, m, scene);
			depth = 1280;	// la profondeur devrait être mise à jour dans la méthode Object::intersect()
		}
	}
	else
	{
		// @@@@@@ VOTRE CODE ICI
		// Notez que pour Object::intersect(), le paramètre hit correspond à celui courant.
		// Votre intersect() devrait être implémenté pour exclure toute intersection plus lointaine que hit.depth
		

		Intersection hit;
        hit.depth = depth;
        
        for (int i = 0 ; i < scene.objects.size(); i++) {
            if (scene.objects[i]->intersect(ray, hit)){
				
					depth = fmin(hit.depth, depth);
					outColor = shade(ray, rayDepth, hit, scene.objects[i]->material, scene);
					
				

            }
        }
		return true;
	}

    // Décrémente la profondeur du rayon.
    rayDepth--;

    return false; 
}


Vector Raytracer::shade(Ray const &ray,
                 int &rayDepth,
                 Intersection const &intersection,
                 Material const &material,
                 Scene const &scene)
{
    // - itérer sur toutes les sources de lumières, calculant les contributions ambiant/diffuse/speculaire
    // - utiliser les rayons d'ombre pour déterminer les ombres
    // - intégrer la contribution de chaque lumière
    // - inclure l'émission du matériau de la surface, s'il y a lieu
    // - appeler Raytracer::trace pour les couleurs de reflection/refraction
    // Ne pas réfléchir/réfracter si la profondeur de récursion maximum du rayon a été atteinte !
	//!!! NOTE UTILE : facteur d'atténuation = 1.0 / (a0 + a1 * d + a2 * d * d)..., la lumière ambiante ne s'atténue pas, ni n'est affectée par les ombres
	//!!! NOTE UTILE : pour chaque type de rayon, i.e. rayon d'ombre, rayon reflechi, et rayon primaire, les profondeurs maximales sont différentes
	Vector diffuse(0);
	Vector ambient(0);
	Vector specular(0);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(1, 1000);



	// Vector important qui dependent pas du loop
	Vector rayDirection = ray.direction.normalized();
	Vector normal = intersection.normal.normalized();
	
	// Slide 31 reflection = 2((vdot N)N- V)
	Vector reflectionRay = (2 * -rayDirection.dot(normal)) * normal + rayDirection;
	reflectionRay.normalize();


	for (auto lightIter = scene.lights.begin(); lightIter != scene.lights.end(); lightIter++)
	{
		/*
		Lumiere ambiant L = l_a*k_a*s(lambda) --> ** calculer 1 seul fois **
		l_a = intensite ambiante lumiere = lightIter->ambient
		k_a = material propriete reflexion = materiel.ambient
		S(lambda) = couleur

		Reflexion lambertienne  ->  L = k_d*L_l*(N*L(l-lumiere)) * s(lambda)
		L_l = intensite de la lumiere = lightIter->diffuse
		k_d = propriete diffuse du material = material.diffuse

		Reflexion Speculaire = (k_s*L_l (N*H)^n)/ r^2  
		k_s = material.speculaire
		L_l = intensite speculaire = lightIter->specular
		n = index de shininess du material = material.shininess
		H = distance entre camera et objet
		r^2 = distance du point de lumiere 

		diffuse + speculaire == ***SOMMES***
		*/
		
	
	
		// Vecteur lumiere Point source - intersection
		Vector lightRay = (lightIter->position - intersection.position).normalized();
		Vector reflectionLight = (2 * lightRay.dot(normal)) * normal - lightRay;
		reflectionLight.normalize();

		
		// donnnees importante
		double oneOverR2 = 1 / (lightIter->position - intersection.position).length2();
		double nDotL = max(normal.dot(lightRay), 0.0);


		Vector ambientIndex = material.ambient * lightIter->ambient;
		Vector diffuseIndex = material.diffuse * lightIter->diffuse *
				nDotL * oneOverR2;
		Vector specularIndex = material.specular * lightIter-> specular *
			 powf(max(0.0, reflectionRay.dot(lightRay)), material.shininess)  * oneOverR2;
		
		
		
		// shadowRay = intersection -> source lumiere
		// Distance ||source - hit||
		/*
		*******************
		SHADOW SUPERSAMPLING
		*********************
		* ne fonctionne pas parfaitement, mais voici ma solution pour des points partiels
		* changer le sampleNum pour avoir plus de super sampling
		*/
		Intersection hitShadow;	
		double sampleNum = 1;
		double radius = lightIter->radius;
		if (radius < 1) {
			sampleNum = 1;
		}
			for (int i = 0; i < sampleNum; i++)
			{
				Ray shadowRay;
				

				// petit pico noir bizarre si on ajoute pas un epsilon
				if (sampleNum == 1) {
					shadowRay = Ray(intersection.position + (lightRay * EPSILON), lightRay);


				}
				else {

					double x = 1000;
					double y = 1000;
					double z = 1000;
					while (true) {

						x = (double)dis(gen) / 1000 * 2 - 1;
						y = (double)dis(gen) / 1000 * 2 - 1;
						z = (double)dis(gen) / 1000 * 2 - 1;
						//cout << x << endl;
						double normaliseXYZ = 1 / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
						//cout << normaliseXYZ << endl;
						double distro = powf((double)dis(gen) / 1000, 1 / 3);
						//cout << distro << endl;

						if (sqrt(x * x + y * y + z * z) < radius) { break; }
					}

					x = x + lightIter->position[0];
					y = y + lightIter->position[1];
					z = z + lightIter->position[2];


					Vector bob = Vector(x, y, z);
					Vector bobRay = (bob - intersection.position).normalized();

					shadowRay = Ray(intersection.position + bobRay * EPSILON, bobRay);

				}

				// les elements ne commence pas dans lombre par default
				bool inShade = false;

				//calculer la distance max entre objet et source de lumiere
				//!!! NOTE UTILE : n'acceptez pas les intersection des rayons
				// d'ombre qui sont plus loin que la position de la lumière

				double maxShadow = (intersection.position - lightIter->position).length() - radius;

				for (int j = 0; j < scene.objects.size(); ++j) {
					if (scene.objects[j]->intersect(shadowRay, hitShadow)) {

						// distance entre intersection et point qui donne l'ombrage
						double objectDistance = (intersection.position - hitShadow.position).length();
						if (objectDistance >= maxShadow) { break; }

						inShade = true;
					}
				}

				// On additionne pas ambiant, pareil partout
				ambient = ambientIndex;



				if (!inShade) {

					diffuse += diffuseIndex / sampleNum;
					specular += specularIndex / sampleNum;
				}
				

			}
		

	}


	/*
	***************
	Texture MAPPING 
	***************
	*/

	if (material.texture.height() != 0) {
		int height = material.texture.height();
		int width = material.texture.width();
		int pointX = round(intersection.position[0] * width/12 -.5 );
		int pointY = round(intersection.position[2] * height/12-.5);
		int x = pointX % width >= 0 ? pointX % width : abs(pointX % width) ;
		int y = pointY % height >= 0 ? pointY % height : abs(pointY % height) ;
		rgb_t foo = material.texture.get_pixel(x, y);
		int red = foo.red;
		int green = foo.green;
		int blue = foo.blue;

		ambient = Vector(red, green, blue, 1) * .0025;
	}



	/*
	* ******************************************************
	* *** REFRACTION +++ RECURSION ++ TRANSPARENCY ****
	* ******************************************************
	*/
	
	// les indices de refraction, reflection, transparence
	Vector reflectedLight(0);
	Vector refractedLight(0);
	Vector transparentLight(0);

	if (rayDepth < 10)
	{
		//1280 par default car c'est ce quils utilisent plus haut
		double depth = 1280;
		
		// REFRACTION
		if (material.refract == 0) {
			int x = 0;
		}
		else if (material.refract != 1 && rayDepth<10) {

			double ref = 1 / material.refract;
			double nDotV = normal.dot(-ray.direction) < 0 ? -normal.dot(ray.direction) : normal.dot(ray.direction);
			
			Vector refract = normal * (ref * nDotV - sqrtf(1 - powf(ref, 2) * (1 - powf(nDotV, 2)))) - ref * -rayDirection;

			Ray refractRay(intersection.position + refract * EPSILON, refract);


			trace(refractRay, rayDepth, scene, refractedLight, depth);

		}

		//RECURSION
		if (!(ABS_FLOAT(material.reflect) < EPSILON) && material.refract == 0 ) {
			//Reflexion calcule page 30-31 des notes de cours aucune refraction ici je crois
			Ray recursiveRay(intersection.position + reflectionRay * EPSILON, reflectionRay);
			

			//trace avec les parametre donnee 
			trace(recursiveRay, rayDepth, scene, reflectedLight, depth);
		}

		//TRANSPARENCE
		if (material.transparency != 0 && material.refract == 0) {
			double refIndex = 1 ;
			double nDotV = normal.dot(-rayDirection) < 0 ? -normal.dot(-rayDirection) : normal.dot(-rayDirection);
			Vector refract = normal * (refIndex * nDotV - sqrtf(1 - powf(refIndex, 2) * (1 - powf(nDotV, 2)))) - refIndex * -rayDirection;
			Ray refractRay(intersection.position + refract * EPSILON, refract);


			trace(refractRay, rayDepth, scene, refractedLight, depth);
		}



	}

		// ce que nous retournons pour faire la couleur --> ajouter refraction + transparence
		return material.emission + ambient + diffuse + specular + (material.reflect * reflectedLight) +
			refractedLight + material.transparency * transparentLight;


}