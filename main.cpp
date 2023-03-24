#include<cmath>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<vector>
#include<utility>
#include<cstdlib>
#include"CImg/CImg.h"
using namespace cimg_library;

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
#define TILING_FACTOR 0.05								// how long is a pixel in the world coordinates

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
	double radiance[400]; // 380 to 779, unit is nanometer
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
};


double gamma;												// coefficient 2
double global_time;											// the time when the image is captured
RGB res_matrix[1920][1080];									// the image matrix, where (0, 0) is the top left corner of the rendering screen (from the observer side).
RGB sample_matrix[1920][1080];								// the unprocessed color matrix for rendering
double color_table[400][3], buf[400][3];
spectrum_1step spec1;

void read_color_table()
{
	std::ifstream is("color_match_table.txt");
	if (is.bad())
	{
		std::cout << "failed reading color matching table file!";
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
}

inline double length(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}

Vector3 get_original_pos_in_plane(Vector3& direction)
{
	double lambda = (C * global_time) / (1.0 + (direction.x / Beta));// coefficient for the direction vector. formula from paper
	return direction * lambda;
}

void sample(CImg<unsigned char>& texture)
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
			Vector3 original_pos = get_original_pos_in_plane(ray_direction);
			while (original_pos.y < 0)original_pos.y += texture._width * TILING_FACTOR;
			while (original_pos.z < 0)original_pos.z += texture._height * TILING_FACTOR;
			int y = original_pos.y / TILING_FACTOR;
			int z = original_pos.z / TILING_FACTOR;
			sample_matrix[y][z].R = texture(y, z, 0, 0);
			sample_matrix[y][z].G = texture(y, z, 0, 1);
			sample_matrix[y][z].B = texture(y, z, 0, 2);
		}
	}
}

int main()
{
	gamma = pow((1.0 - Beta * Beta), -0.5);
	CImg<unsigned char> res_img(RES_WIDTH, RES_HEIGHT, 0, 3, 0);


}