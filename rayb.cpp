#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
using namespace std;


/*
The majority code is taken from www.scratchapixel.com
It is licenced under the GNU General Public Licence.
I do not own any rights on the code except the modified parts.

Code is modified for CS405 Ray Tracing homework.
Added configuration info and modified functions in order to understand the basics
of the Ray Tracing. 
Most of the parts will be commented based on their sources on various websites.
What I did was, re-writing the whole code with looking the geometric functions
and tried to implement the functions with different parameters.

In order to run the program you need read README and fill the input based on the instructions.
-Width
-Height
-Field of View
-Object informations
should be entered in order.

Berkant Deniz Aktaş
October 2018
*/


#if defined __linux__ || defined __APPLE__
#else
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

template<typename T>
class Vec
{
public:
	T x,y,z;
	Vec(){ x=T(0),y=T(0),z=T(0);} 							// Creating 0,0,0
	Vec(T p){x=p,y=p,z=p;}									// Creating p,p,p
	Vec(T xVal, T yVal, T zVal){x=xVal, y=yVal, z=zVal;}	// Creating x,y,z

	Vec &normalize()										//Normalizing the vector
	{
		T normalSquare = x*x + y*y + z*z;					//x^2 + y^2 + z^2
		if(normalSquare>0)
		{
			T inverseNormal = 1 / sqrt(normalSquare);
			x = x * inverseNormal;
			y = y * inverseNormal;
			z = z * inverseNormal;
		}
		return * this;
	}

	Vec<T> operator * (const T &value) const 				//Multiply by a value
	{return Vec<T>(x*value, y*value, z*value);}
	Vec<T> operator * (const Vec<T> &vec) const				//Multiply with a vector 
	{return Vec<T>(x*vec.x, y*vec.y, z*vec.z);}
	T dotProduct(const Vec<T> &vec) const					//Scalar product of vectors
	{return x*vec.x + y*vec.y + z*vec.z;}
	Vec<T> operator - (const Vec<T> &vec) const				//Substract from a vector
	{return Vec<T>(x-vec.x, y-vec.y, z-vec.z);}
	Vec<T> operator + (const Vec<T> &vec) const				//Add to a vector
	{return Vec<T>(x+vec.x, y+vec.y, z+vec.z);}
	Vec<T> &operator += (const Vec<T> &vec) 				//Add to itself
	{ Vec(x+=vec.x, y+=vec.y, z+=vec.z); return*this;}
	Vec<T> &operator *= (const Vec<T> &vec) 
	{ Vec(x*=vec.x, y*=vec.y, z*=vec.z); return*this;}		//Multiply by itself
	Vec<T> operator - () const { return Vec<T>(-x,-y,-z);}	//Inverse direction of a vector
	friend ostream & operator << (ostream &os, const Vec<T> &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};
typedef Vec<float> Vecf;
class Sphere 												//Sphere object defined with
{															//Color,Emission(light),Trans,Reflection
public:
	Vecf center, surfaceColor, emissionColor;
	float radius, radiusSquare, transparency, reflection;
	Sphere(const Vecf &c, const float &r, const Vecf &sc, const float &refl=0, // Constructor of a sphere
		const float &transp=0, const Vecf &ec =0)
	{
		center = c,
		surfaceColor = sc,
		emissionColor = ec,
		radius = r,
		radiusSquare = r * r,
		transparency = transp,
		reflection = refl;
	}

	bool intersect(const Vecf &rayOrigin, const Vecf &rayDirection, float &t0, float &t1) const // Intersection of a sphere and ray
	{	
		// two solutions exist - geometric and analytic this one is geometric
		Vecf l = center - rayOrigin;   				// formula information : 
		float tca = l.dotProduct(rayDirection);		// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
		if(tca<0)									//tca and other parameters are visualised in website
		{
			return false;
		}
		float d2 = l.dotProduct(l) - tca * tca;
		if(d2>radiusSquare)
		{
			return false;
		}
		float thc = sqrt(radiusSquare-d2);
		t0 = tca - thc; 
        t1 = tca + thc;
        
        return true;
	}
};

#define MAX_RAY_DEPTH 5 // In order to limit the rendering time reflection and transparency recursion
float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
} // Problematic equation that I could not find explaination 
	//TODO: Find and implement it by hand

Vecf trace(const Vecf &rayOrigin,	const Vecf &rayDirection, const vector<Sphere> &sphereVector, const int &depth)
{
	float tnear = INFINITY;
	const Sphere *sphere = NULL;

	for(unsigned i=0; i<sphereVector.size(); ++i)
	{
		float t0 = INFINITY; // A ray can intersect at two points, take the closest one 
		float t1 = INFINITY;
		if(sphereVector[i].intersect(rayOrigin, rayDirection, t0, t1)) // Check for intersection for each object 
		{
			if(t0<0)
			{
				t0=t1;
			}
			if(t0<tnear)
			{
				tnear = t0;
				sphere = &sphereVector[i];
			}
		}
	}
	if(!sphere) // if it does not intersect return black or bgc
	{
		return Vecf(2);
	}
	Vecf surfaceColor = 0; 
	Vecf phit = rayOrigin + rayDirection*tnear; // intersection point 
	Vecf nhit = phit - sphere->center;			// normal of the intersection point 
	nhit.normalize();

	float bias=1e-4; //Adding a bias 
					//TODO: Find, why we are adding a bias 
	bool inside = false;
	if(rayDirection.dotProduct(nhit)>0)
	{
		nhit = -nhit;
		inside = true;
	}
	if((sphere->transparency>0 || sphere->reflection>0) 
		&&
		depth < MAX_RAY_DEPTH)
	{
		// This is advanced part could be deleted if transparency and reflection is not wanted
		// detailed information about reflection, refraction and fresnel effect taken from : 
		// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
		float facingRatio = -rayDirection.dotProduct(nhit);
		float fresnelEffect = mix(pow(1- facingRatio,3),1,0.1); //Adding a fresnelEffect,
		Vecf reflectionDirection = rayDirection - nhit * 2 * rayDirection.dotProduct(nhit);
		reflectionDirection.normalize();
		Vecf reflection = trace(phit + nhit * bias, reflectionDirection, sphereVector, depth+1);
		Vecf refraction = 0;
		if(sphere->transparency) // Good example slides for this equations
			// http://web.cs.wpi.edu/~emmanuel/courses/cs563/S10/talks/wk11_p2_wadii_simple_transparency.pdf
			// taken from Worcester Polytechnic Institute (WPI) CS 563 Advanced Topics in Computer Graphics
			// Credit : Wadii Bellamine
			// idea is tracing again and again, transmission   - refraction ray
		{
			float ior = 1.1;
			float eta = (inside) ? ior : 1 / ior;
			float cosi = -nhit.dotProduct(rayDirection);
			float k = 1 - eta * eta * (1-cosi*cosi);
			Vecf refractionDirection = rayDirection * eta + nhit * (eta * cosi - sqrt(k)); 
			refractionDirection.normalize();
			refraction = trace(phit - nhit - bias, refractionDirection, sphereVector, depth+1);
		}
		surfaceColor = (reflection * fresnelEffect + refraction * (1-fresnelEffect) * sphere-> transparency)
		* sphere->surfaceColor; // Could be printed with all details with this crazy geometric equations

		//surfaceColor =  sphere->surfaceColor; // delete reflection refraction transparency and just set shadows 
	}
	else
	{
		// stop tracing if diffuse object
		for(unsigned i = 0; i<sphereVector.size();++i)
		{
			if(sphereVector[i].emissionColor.x>0) // light !
			{
				Vecf transmission = 1;
				Vecf lightDirection = sphereVector[i].center - phit;
				lightDirection.normalize();
				for(unsigned j = 0; j<sphereVector.size(); ++j)
				{
					if(i!=j)
					{
						float t0,t1;
						if(sphereVector[j].intersect(phit+nhit*bias, lightDirection, t0, t1))
						{
							transmission = 0;
							break;
						}
					}
				}
				surfaceColor += sphere->surfaceColor * transmission 
				* max(float(0), nhit.dotProduct(lightDirection)) * sphereVector[i].emissionColor;
			}
		}
	}
	return surfaceColor + sphere->emissionColor;
}

void render(const vector<Sphere> &sphereVector, int W, int H, float fieldOfView )
// main tracing function setted a window for camera 
// camera functions are taken from : 
// https://www.scratchapixel.com/code.php?id=7&origin=/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays
// modified main code to set FoV, if it is greater than 180, it inverses the camera !
{
	unsigned windowWidth = (unsigned) W;
	unsigned windowHeight = (unsigned) H;
	Vecf *image = new Vecf[windowWidth*windowHeight];
	Vecf *pixel = image;
	float inverseWidth = 1/(float)windowWidth;
	float inverseHeight = 1/(float)windowHeight;
	//float fov = 30;
	float aspectRatio = windowWidth / float(windowHeight);
	float angle = tan(3.141592653589793 * 0.5 * fieldOfView / 180.);

	for(unsigned y = 0; y < windowHeight; ++y)
	{
		for(unsigned x = 0; x < windowWidth; ++x, ++pixel)
		{
			float xx = (2*((x+0.5) * inverseWidth)-1) * angle * aspectRatio;
			float yy = (1-2 * ((y+0.5) * inverseHeight)) * angle;
			Vecf rayDirection(xx,yy,-1);
			rayDirection.normalize();
			*pixel = trace(Vecf(0),rayDirection,sphereVector,0);
		}
	}
	// it saves title as ppm
	// TODO: TIFF format 
	std::ofstream ofs("./untitled.TIFF", std::ios::out | std::ios::binary);
	ofs << "P6\n" << windowWidth << " " << windowHeight << "\n255\n";
	 for (unsigned i = 0; i < windowWidth * windowWidth; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
               (unsigned char)(std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
    delete [] image;

}



int main(int argc, char **argv)
{
	// Main is mostly written by me, it has basic I/O operations
	ifstream getInput;
	int windowHeight,windowWidth;
	float lightX, lightY, lightZ, lightRadius,fieldOfView; 
	vector<Sphere> sphereVector;

	string filename = "input.txt";
	// TODO: open it from command line
	getInput.open(filename);
	getInput >> windowWidth;
	getInput >> windowHeight;
	getInput >> fieldOfView;
	getInput >> lightX >> lightY >> lightZ >> lightRadius;

	sphereVector.push_back(Sphere(Vecf(lightX,lightY,lightZ),lightRadius,Vecf(0,0,0),0,0.0,Vecf(3)));

	float sphereX, sphereY, sphereZ, sphereRadius, sphereReflect, sphereTrans;
	float sphereRedValue, sphereGreenValue, sphereBlueValue; 
	while(getInput >> sphereX >>sphereY >> sphereZ >> sphereRadius 
		>> sphereRedValue >> sphereGreenValue >>sphereBlueValue 
		>> sphereReflect >> sphereTrans)
	{
		sphereRedValue = sphereRedValue/255;		// changed RGB - 1/0 values for easy input entering
		sphereGreenValue = sphereGreenValue/255;
		sphereBlueValue = sphereBlueValue/255;

		sphereVector.push_back(Sphere(Vecf(sphereX,sphereY,sphereZ),sphereRadius,Vecf(sphereRedValue,sphereGreenValue,sphereBlueValue),sphereReflect,sphereTrans));

		
	}
	render(sphereVector,windowWidth,windowHeight,fieldOfView); // renders the image
	
	return 0;
}