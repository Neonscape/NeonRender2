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
#define IlluminantD65 0.3127, 0.3291					
#define GAMMA_REC709 0

// structures

struct RGB
{
	unsigned char R;
	unsigned char G;
	unsigned char B;
};
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
struct color_system
{
	double xred, yred,
		xgreen, ygreen,
		xblue, yblue,
		xwhite, ywhite,
		gamma;
};


double gamma;												// coefficient 2
double global_time;											// the time when the image is captured
RGB res_matrix[1920][1080];									// the image matrix, where (0, 0) is the top left corner of the rendering screen (from the observer side).
bool sample_matrix[1920][1080];								// the unprocessed matrix for rendering
double color_table[400][3], buf[400][3];
spectrum_1step spec1, spec2;
color_system SMPTE_SYSTEM = {
	0.630, 0.340,
	0.310, 0.595,
	0.155, 0.070,
	IlluminantD65,
	GAMMA_REC709
};

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
		res.radiance[wl_int] = sp.radiance[i] * headlight_coeff * (1 - wl_remain);
		if (wl_int < 399)
		{
			res.radiance[wl_int + 1] = sp.radiance[i] * headlight_coeff * wl_remain;
		}
	}
	return res;
}
RGB spectrum_to_rgb(spectrum_1step& const sp, color_system& const cs)
{
	double	xr, xg, xb,
		yr, yg, yb,
		zr, zg, zb,
		xw, yw, zw,
		rw, gw, bw,
		rvxr, rvxg, rvxb,
		rvyr, rvyg, rvyb,
		rvzr, rvzg, rvzb;

	// the original matrix
	//
	// xr xg xb		R		X
	// yr yg yb  *	G   =	Y
	// zr zg zb		B		Z
	//

	xr = cs.xred;
	xg = cs.xgreen;
	xb = cs.xblue;
	xw = cs.xwhite;
	yr = cs.yred;
	yg = cs.ygreen;
	yb = cs.yblue;
	yw = cs.ywhite;
	zr = 1 - (xr + yr);
	zg = 1 - (xg + yg);
	zb = 1 - (xb + yb);
	zw = 1 - (xw + yw);

	rvxr = (yg * zb - yb * zg);
	rvxg = -(yr * zb - yb * zr);
	rvxb = (yr * zg - yg * zr);
	rvyr = -(xg * zb - zg * xb);
	rvyg = (xr * zb - xb * zr);
	rvyb = -(xr * zg - xg * zr);
	rvzr = (xg * yb - xb * yg);
	rvzg = -(xr * yb - xb * yr);
	rvzb = (xr * yg - xg * yr);

	rw = (rvxr * xw + rvyr * yw + rvzr * zw) / yw;
	gw = (rvxg * xw + rvyg * yw + rvzg * zw) / yw;
	bw = (rvxb * xw + rvyb * yw + rvzb * zw) / yw;

	rvxr /= rw;
	rvyr /= rw;
	rvzr /= rw;
	rvxg /= gw;
	rvyg /= gw;
	rvzg /= gw;
	rvxb /= bw;
	rvyb /= bw;
	rvzb /= bw;

	for (int i = 0; i < 400; i++)
	{
		double X, Y, Z;
		X = sp.radiance[i] * color_table[i][0];
		Y = sp.radiance[i] * color_table[i][1];
		Z = sp.radiance[i] * color_table[i][2];

		Vector3 res = { 0,0,0 }; // here the x/y/z corresponds to R/G/B.
		res.x = rvxr * X + rvyr * Y + rvzr * Z;
		res.y = rvxg * X + rvyg * Y + rvzg * Z;
		res.z = rvxb * X + rvyb * Y + rvzb * Z;

		//compensate the values that cannot be represented in the current color systems (indicated as negative values) with white.
#define Min(a, b) a < b ? a : b
		double white_compensate = Min(0, Min(res.x, Min(res.y, res.z)));
		white_compensate = -1;
		if (white_compensate > 0)
		{
			res.x += white_compensate;
			res.y += white_compensate;
			res.z += white_compensate;
		}

		//normalize RGB values to integer within [0, 255]
#define Max(a, b) a > b ? a : b
		double mx = Max(Max(res.x, res.y), res.z);
		res.x /= mx;
		res.y /= mx;
		res.z /= mx;

		res.x = (int)(255 * res.x);
		res.y = (int)(255 * res.y);
		res.z = (int)(255 * res.z);

		return { (unsigned char)res.x, (unsigned char)res.y, (unsigned char)res.z };

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
			ray_direction = { posx, posy, posz };
			ray_direction.normalize();
			Vector3 hit = get_hitpoint(ray_direction);
			Vector3 original_pos = get_original_point(hit);
			if (((original_pos.y - (int)original_pos.y) <= 0.005) && ((original_pos.z - (int)original_pos.z) <= 0.005))
			{
				res_matrix[i + RES_WIDTH >> 1][RES_HEIGHT >> 1 - j] = { 255, 255, 255 };
				continue;
			}
			int yp = original_pos.y / TILING_FACTOR;
			int zp = original_pos.z / TILING_FACTOR;
			spectrum_1step temp;
			if (!((yp + zp) & 1))
			{
				temp = convert_spectrum(hit, spec1);
			}
			else
			{
				temp = convert_spectrum(hit, spec2);
			}
			res_matrix[i + RES_WIDTH >> 1][RES_HEIGHT >> 1 - j] = spectrum_to_rgb(temp, SMPTE_SYSTEM);
		}
	}
}
void output_image(CImg<unsigned char>& img)
{
	for (int i = 0; i < RES_WIDTH; i++)
	{
		for (int j = 0; j < RES_HEIGHT; j++)
		{
			img(i, j, 0, 0) = res_matrix[i][j].R;
			img(i, j, 0, 1) = res_matrix[i][j].G;
			img(i, j, 0, 2) = res_matrix[i][j].B;
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
	sample();
	output_image(res_img);
	res_img.save("result.jpg");
	return 0;
}