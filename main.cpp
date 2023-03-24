#include<cmath>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<vector>
#include<utility>
#include<cstdlib>
#include"CImg/CImg.h"
#include<cstring>
#include<string>
using namespace cimg_library;
using namespace std;

// defines

#define RES_WIDTH 1920									//
#define RES_HEIGHT 1080									// the resolution of the final picture
#define C 1												// speed of light
#define V ((double)C * 0.8)								// velocity of the plane
#define Beta ((double)V  / (double)C)					// coefficient
#define DEFAULT_COLOR {60, 60, 60}						// color used when the sampling ray fails to hit the target
#define INITIAL_TIME -50								// the 0 is when the plane moves to the origin
#define SCREEN_DISTANCE 1.0								// the distance of the rendering screen to the observer
#define SCREEN_WIDTH  SCREEN_DISTANCE * 39.6 / 180.0	// the actual width of the rendering screen, here I used the 35mm Movie Film Standard.
#define SCREEN_HEIGHT SCREEN_DISTANCE * 27.0 / 180.0	// the actual height of the rendering screen.
#define TILING_FACTOR 0.5								// how long is a color square in the world coordinates

// structures

struct RGB
{
	unsigned char R;
	unsigned char G;
	unsigned char B;
};
#
struct spectrum_1step
{
	double radiance[400] = { 0 }; // 380 to 779, unit is nanometer
	double scaling_factor = 0;

	double operator()(int const& wavelength)const
	{
		return radiance[wavelength - 380] * scaling_factor;
	}
};
struct Vector3
{
	double x, y, z;
	Vector3 operator * (double const& a)const
	{
		return { x * a, y * a, z * a };
	}
	double operator * (Vector3 const& a)const
	{
		return x * a.x + y * a.y + z * a.z;
	}
	bool is_normalized()
	{
		return ((x * x + y * y + z * z - 1) <= 0.0005);
	}
	void normalize()
	{
		double mod = x * x + y * y + z * z;
		mod = sqrt(mod);
		x /= mod;
		y /= mod;
		z /= mod;
	}
	double mod()
	{
		return sqrt(x * x + y * y + z * z);
	}
};

double gamma;												// coefficient 2
double global_time;											// the time when the image is captured
RGB res_matrix[1920][1080];									// the image matrix, where (0, 0) is the top left corner of the rendering screen (from the observer side).
bool sample_matrix[1920][1080];								// the unprocessed matrix for rendering
double color_table[400][3], buf[400][3];
spectrum_1step spec1, spec2;

void read_spectrum(int i)
{
	string s = "spectrum_";
	s += (char)(i + '0');
	s += ".txt";
	ifstream is(s);
	if (is.bad())
	{
		cerr << "failed reading spectrum preset file " << i << "!";
		exit(1);
	}
	spectrum_1step* pt = nullptr;
	if (i == 1)
	{
		pt = &spec1;
	}
	else
	{
		pt = &spec2;
	}
	for (int i = 0; i < 400; i++)
	{
		cin >> pt->radiance[i];
	}
	cin >> pt->scaling_factor;
	is.close();
}
void read_color_table()
{
	ifstream is("color_match_table.txt");
	if (is.bad())
	{
		cerr << "failed reading color matching table file!";
		exit(1);
	}
	double temp;
	for (int i = 0; i < 400; i++)
	{
		std::cin >> temp;
		std::cin >> color_table[i][0];
		std::cin >> color_table[i][1];
		std::cin >> color_table[i][2];
	}
	is.close();
}

spectrum_1step convert_spectrum(Vector3 pos, spectrum_1step& const sp)
{
	spectrum_1step res;
	memset(res.radiance, 0, sizeof(res.radiance));
	res.scaling_factor = sp.scaling_factor;
	double headlight_coeff = (Beta * C * global_time) * (Beta * C * global_time) / (gamma * Beta * C * global_time - pos.x / gamma) / (gamma * Beta * C * global_time - pos.x / gamma) / gamma;
	for (int i = 0; i < 400; i++)
	{
		double freq_orig = (((double)i + 380.0) / 1e9) / C;
		double freq_trans = freq_orig * (Beta * C * global_time - pos.x) / (gamma * Beta * C * global_time - pos.x / gamma);
		double wl = freq_trans * C * 1e9 - 380;
		int wl_int = wl;
		double wl_remain = wl - wl_int;
		if (wl_int >= 400)break;


	}
}

inline double length(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}

Vector3 get_hitpoint(Vector3& direction) // get the location of the point where the ray intersects with the transformed object
{
	if (!direction.is_normalized())direction.normalize();
	double lambda = (C * global_time) / (1.0 + (direction.x / Beta));// coefficient for the direction vector. formula from paper
	return direction * lambda;
}
Vector3 get_original_point(Vector3& const hit_pos)
{
	Vector3 res = { 0,0,0 };
	res.x = -gamma * Beta * (C * global_time - hit_pos.mod()) + gamma * hit_pos.x;
	res.y = hit_pos.y;
	res.z = hit_pos.z;
	return res;
}

void sample()
{

	// the sampling process, where right is positive y and up is positive z.

	Vector3 ray_direction = { 0,0,0 }, target_pos = { 0,0,0 };

	int start_y = -(RES_WIDTH >> 1);
	int start_z = (RES_HEIGHT >> 1);
	int end_y = RES_WIDTH + start_y;
	int end_z = -RES_HEIGHT + start_z;

	for (int i = start_y; i <= end_y; i++) // from -0.5 to 0.5
	{
		for (int j = start_z; j >= end_z; j--) // from 0.5 to -0.5
		{
			double posy = SCREEN_WIDTH * ((double)i / (double)RES_WIDTH);
			double posz = SCREEN_HEIGHT * ((double)j / (double)RES_HEIGHT);
			double posx = SCREEN_DISTANCE;
			double len = length(posx, posy, posz);
			ray_direction = { posx / len, posy / len, posz / len };
			Vector3 hit = get_hitpoint(ray_direction);
			Vector3 original_pos = get_original_point(hit);
			int yp = original_pos.y / TILING_FACTOR;
			int zp = original_pos.z / TILING_FACTOR;
			if (!((yp + zp) & 1))
			{
				sample_matrix[i + RES_WIDTH >> 1][RES_HEIGHT - j] = 1;
			}
		}
	}

}

int main()
{
	read_color_table();
	read_spectrum(1);
	read_spectrum(2);
	gamma = pow((1.0 - Beta * Beta), -0.5);
	CImg<unsigned char> res_img(RES_WIDTH, RES_HEIGHT, 0, 3, 0);


}