#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gd.h>
#include <stdbool.h>
#include <time.h>

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define MAX_RAY_DEPTH 5
#define M_PI 3.141592653589793

gdImagePtr image;
int noOfEllipse;

typedef struct
{
    double x;
    double y;
	double z;
} vect;

typedef struct 
{
	vect* origin;
	vect* direction;
} ray;

typedef struct
{
	vect* center;
	double radius;
	double transparency;
	double reflection;
	vect* surfaceColor;
	vect* emission;
} ellipse;

double mix(double a, double b, double c) 
{ 
    return b * c + a * (1 - c); 
} 

double length(vect* vector)
{
	double length = (vector->x * vector->x) + (vector->y * vector->y) + (vector->z * vector->z);
	return sqrt(length);
}

void normalize(vect* vector)
{
	if(length(vector) != 0)
	{
		double inverseLength = 1/length(vector);
		vector->x *= inverseLength;
		vector->y *= inverseLength;
		vector->z *= inverseLength;
	}
}

double dot(vect* vect1, vect* vect2)
{
	return (vect1->x * vect2->x) + (vect1->y * vect2->y) + (vect1->z * vect2->z);
}

ellipse* createEllipse(double centerX, double centerY, double centerZ, double radius, double transp, double ref, vect* col, vect* emiss)
{
	ellipse* ell = malloc(sizeof(vect*) * 4 + sizeof(double) * 3);
	vect* vector = malloc(sizeof(vect*));
	
	vector->x = centerX; vector->y = centerY; vector->z = centerZ;

	ell->center = vector;
	ell->radius = radius;
	ell->transparency = transp;
	ell->reflection = ref;
	ell->surfaceColor = malloc(sizeof(vect*));
	ell->surfaceColor = col;
	ell->emission = malloc(sizeof(vect*));
	ell->emission = emiss;

	return ell;
}

bool intersect(ellipse* ell, vect* rayorig, vect* raydir, double* t0, double* t1)
{
	vect* vector = malloc(sizeof(vect*));
	vector->x = ell->center->x - rayorig->x;
	vector->y = ell->center->y - rayorig->y;
	vector->z = ell->center->z - rayorig->z;

	double tca = dot(vector, raydir);
	
	if(tca < 0)
	{
		return false;
	}

	double x = dot(vector, vector) - (tca*tca);

	if(x > (ell->radius * ell->radius))
	{
		return false;
	}

	double thc = sqrt((ell->radius * ell->radius) - x); 
    *t0 = tca - thc; 
    *t1 = tca + thc;

	return true;
}

vect* trace(vect* rayorig, vect* raydir, ellipse* ellipses[], int depth)
{ 
    double tnear = INFINITY; 
	ellipse* ell = NULL;
	vect* surfaceColor = malloc(sizeof(vect*)); 

	for (int i = 0; i < noOfEllipse; i++) 
	{ 
    	double t0 = INFINITY, t1 = INFINITY; 
    	if (intersect(ellipses[i], rayorig, raydir, &t0, &t1)) { 
        	if (t0 < 0)
			{ 
				t0 = t1;
			}

        	if (t0 < tnear) 
			{ 
            	tnear = t0; 
            	ell = ellipses[i]; 
        	} 
    	} 
	} 
    
    if (!ell)
	{ 
		surfaceColor->x = 2;
		surfaceColor->y = 2;
		surfaceColor->z = 2;
		return surfaceColor;
	} 
    surfaceColor->x = 0; surfaceColor->y = 0; surfaceColor->z = 0; 
	
	vect* phit = malloc(sizeof(vect*));
	phit->x = rayorig->x + raydir->x * tnear;
	phit->y = rayorig->y + raydir->y * tnear;
	phit->z = rayorig->z + raydir->z * tnear;

	vect* nhit = malloc(sizeof(vect*));
	nhit->x = phit->x - ell->center->x;
	nhit->y = phit->y - ell->center->y;
	nhit->z = phit->z - ell->center->z;  
	normalize(nhit);  
	
	double bias = 1e-4; 
	bool inside = false; 
    if (dot(raydir, nhit) > 0.0)
	{ 
		nhit->x = -nhit->x; nhit->y = -nhit->y; nhit->z = -nhit->z; 
		inside = true; 
	} 
    if ((ell->transparency > 0 || ell->reflection > 0) && depth < MAX_RAY_DEPTH) 
	{ 
    	double facingratio = -dot(raydir, nhit); 
        
        double fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1); 
        
        vect* refldir = malloc(sizeof(vect*)); 
		refldir->x = raydir->x - nhit->x * 2 * dot(raydir, nhit); 
		refldir->y = raydir->y - nhit->y * 2 * dot(raydir, nhit);
		refldir->z = raydir->z - nhit->z * 2 * dot(raydir, nhit);
        normalize(refldir); 

		vect* reflection = malloc(sizeof(vect*));
		reflection->x = phit->x + nhit->x * bias;
		reflection->y = phit->y + nhit->y * bias;
		reflection->z = phit->z + nhit->z * bias;
        reflection = trace(reflection, refldir, ellipses, depth + 1); 
        
		vect* refraction = malloc(sizeof(vect*));
		refraction->x = 0; refraction->y = 0; refraction->z = 0;  

        if (ell->transparency) 
		{ 
            double ior = 1.1, eta = (inside) ? ior : 1 / ior; 
            double cosi = -dot(nhit, raydir); 
            double k = 1 - eta * eta * (1 - cosi * cosi); 
            vect* refrdir = malloc(sizeof(vect*)); 
			refrdir->x = raydir->x * eta + nhit->x * (eta *  cosi - sqrt(k));
			refrdir->y = raydir->y * eta + nhit->y * (eta *  cosi - sqrt(k));
			refrdir->z = raydir->z * eta + nhit->z * (eta *  cosi - sqrt(k)); 
            normalize(refrdir);

			refraction->x = phit->x + nhit->x * bias;
			refraction->y = phit->y + nhit->y * bias;
			refraction->z = phit->z + nhit->z * bias; 
            refraction = trace(refraction, refrdir, ellipses, depth + 1); 
        }
        surfaceColor->x = ( 
            reflection->x * fresneleffect + 
            refraction->x * (1 - fresneleffect) * ell->transparency) * ell->surfaceColor->x; 
		surfaceColor->y = ( 
            reflection->y * fresneleffect + 
            refraction->y * (1 - fresneleffect) * ell->transparency) * ell->surfaceColor->y;
		surfaceColor->z = ( 
            reflection->z * fresneleffect + 
            refraction->z * (1 - fresneleffect) * ell->transparency) * ell->surfaceColor->z;
    } 
    else 
	{ 
        for (int i = 0; i < noOfEllipse; i++) 
		{ 
            if (ellipses[i]->emission->x > 0) 
			{ 
                vect* transmission = malloc(sizeof(vect*));
				transmission->x = 1; transmission->y = 1; transmission->z = 1;

                vect* lightDirection = malloc(sizeof(vect*));
				lightDirection->x = ellipses[i]->center->x - phit->x;
				lightDirection->y = ellipses[i]->center->y - phit->y;
				lightDirection->z = ellipses[i]->center->z - phit->z; 
                normalize(lightDirection);

                for (int j = 0; j < noOfEllipse; j++) 
				{ 
                    if (i != j) 
					{ 
                        double t0 = 0.0, t1 = 0.0;
						vect* temp = malloc(sizeof(vect*));
						temp->x = phit->x + nhit->x * bias;
						temp->y = phit->y + nhit->y * bias;
						temp->z = phit->z + nhit->z * bias;

                        if (intersect(ellipses[j], temp, lightDirection, &t0, &t1)) 
						{ 
                            transmission->x = 0; transmission->y = 0; transmission->z = 0; 
                            break; 
                        } 
                    } 
                } 
                surfaceColor->x = ell->surfaceColor->x * transmission->x * 
                max(0.0, dot(nhit, lightDirection)) * ellipses[i]->emission->x; 

				surfaceColor->y = ell->surfaceColor->y * transmission->y * 
                max(0.0, dot(nhit, lightDirection)) * ellipses[i]->emission->y;

				surfaceColor->z = ell->surfaceColor->z * transmission->z * 
                max(0.0, dot(nhit, lightDirection)) * ellipses[i]->emission->z;
            } 
        } 
    } 
	vect* temp2 = malloc(sizeof(vect*));
	temp2->x = surfaceColor->x + ell->emission->x;
	temp2->y = surfaceColor->y + ell->emission->y;
	temp2->z = surfaceColor->z + ell->emission->z;
    return temp2; 
} 

void render(ellipse* ellipses[]) 
{ 
    int width = 1000, height = 1000; 
    vect* pixels[width * height]; 
    double invWidth = 1 / (double)width, invHeight = 1 / (double)height; 
    double fov = 30, aspectratio = width / (double)height; 
    double angle = tan(M_PI * 0.5 * fov / 180.0); 
    // Trace rays
	int index = 0;
    for (int y = 0; y < height; y++) { 
        for (int x = 0; x < width; x++) { 
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio; 
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle; 
            vect* raydir = malloc(sizeof(vect));
			raydir->x = xx; raydir->y = yy; raydir->z = -1; 
            normalize(raydir);
			vect* temp = malloc(sizeof(vect*));
			temp->x = 0;
			temp->y = 0;
			temp->z = 0; 
            pixels[index] = trace(temp, raydir, ellipses, 0);
			index++; 
        } 
    } 

	index = 0;
	for (int i = 0; i < height; i++) 
	{
		for(int j = 0; j < width; j++)
		{ 
			int red = (int)(min(1.0, pixels[index]->x) * 255);
			int green =  (int)(min(1.0, pixels[index]->y) * 255);
			int blue = (int)(min(1.0, pixels[index]->z) * 255);
        	int color = red;
			color = (color<<8)+green;
			color = (color<<8)+blue;
			gdImageSetPixel(image, j, i, color);
			index++;
		}
    }
}

void createEllipses(ellipse* ellipses[])
{
	vect* color = malloc(sizeof(vect*));			
	vect* emiss = malloc(sizeof(vect*));

	color->x = 0.90; color->y = 0.76; color->z = 0.46;
	emiss->x = 0; emiss->y = 0; emiss->z = 0;

	ellipse* ell = createEllipse(0.0, -1, -15, 2, 1.0, 0, color, emiss);
	ellipses[0] = ell;

	vect* color2 = malloc(sizeof(vect*));			

	color2->x = 0; color2->y = 1.0; color2->z = 0;

	ellipse* ell2 = createEllipse(1.0, 2.0, -20, 1, 0.5, 0, color2, emiss);
	ellipses[1] = ell2;

	vect* color3 = malloc(sizeof(vect*));			

	color3->x = 1.0; color3->y = 0.0; color3->z = 1.0;

	ellipse* ell3 = createEllipse(-1.0, -1.0, -10, 1, 0.5, 0.5, color3, emiss);
	ellipses[2] = ell3;

	vect* color4 = malloc(sizeof(vect*));			

	color4->x = 0.0; color4->y = 1.0; color4->z = 0.0;

	ellipse* ell4 = createEllipse(2.0, 2.0, -10, 1, 0.5, 0.5, color4, emiss);
	ellipses[3] = ell4;

	vect* color5 = malloc(sizeof(vect*));			

	color5->x = 0.5; color5->y = 1.0; color5->z = 0.25;

	ellipse* ell5 = createEllipse(-2.0, -2.0, -10, 1, 0.5, 1.0, color4, emiss);
	ellipses[4] = ell5;

	vect* color6 = malloc(sizeof(vect*));			

	color6->x = 0.5; color6->y = 0.0; color6->z = 0.25;

	ellipse* ell6 = createEllipse(0.0, 3.0, -10, 1, 0.5, 1.0, color5, emiss);
	ellipses[5] = ell6;

	vect* color7 = malloc(sizeof(vect*));			

	color7->x = 0.0; color7->y = 0.5; color7->z = 0.25;

	ellipse* ell7 = createEllipse(3.0, 1.0, -10, 2, 0.5, 1.0, color7, emiss);
	ellipses[6] = ell7;

	vect* color8 = malloc(sizeof(vect*));			

	color8->x = 0.0; color8->y = 0.5; color8->z = 0.25;

	ellipse* ell8 = createEllipse(-3.0, 2.0, -10, 1, 0.5, 1.0, color8, emiss);
	ellipses[7] = ell8;

	vect* color9 = malloc(sizeof(vect*));			

	color9->x = 0.0; color9->y = 0.5; color9->z = 0.25;

	ellipse* ell9 = createEllipse(3.0, -2.0, -10, 1, 0.5, 1.0, color9, emiss);
	ellipses[8] = ell9;

	vect* light = malloc(sizeof(vect*));
	light->x = 0; light->y = 0; light->z = 0;
	vect* emiss2 = malloc(sizeof(vect*));
	emiss2->x = 3; emiss2->y = 3; emiss->z = 3;

	ellipse* lightBulb = createEllipse(0, 30, -10, 2, 0, 0, light, emiss2);
	ellipses[9] = lightBulb;
}

int main()
{
	FILE *out;

	image = gdImageCreateTrueColor(1000, 1000);

	ellipse* ellipses[10];
	createEllipses(ellipses);

	noOfEllipse = 10;
	render(ellipses);

	out = fopen("spheres.jpg", "w");
	gdImageJpeg(image, out, 100);

	fclose(out);
	gdImageDestroy(image);
}
