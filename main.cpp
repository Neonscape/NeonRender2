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
#include<ctime>
using namespace cimg_library;
using namespace std;

// defines

int RES_WIDTH = 640;									//
int RES_HEIGHT = 360;									// the resolution of the final picture
#define C 1												// speed of light
double speed_coeff = 0.5;
#define V ((double)C * speed_coeff)						// velocity of the plane
#define Beta ((double)V  / (double)C)					// coefficient
#define DEFAULT_COLOR {60, 60, 60}						// color used when the sampling ray fails to hit the target
#define INITIAL_TIME -50								// the 0 is when the plane moves to the origin
#define SCREEN_DISTANCE 0.5							// the distance of the rendering screen to the observer
#define SCREEN_WIDTH  SCREEN_DISTANCE * 39.6 / 180.0	// the actual width of the rendering screen, here I used the 35mm Movie Film Standard.
#define SCREEN_HEIGHT SCREEN_DISTANCE * 27.0 / 180.0	// the actual height of the rendering screen.
#define TILING_FACTOR 0.1								// how long is a color square in the world coordinates
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
double global_time = -20;									// the time when the image is captured
RGB res_matrix[1920][1080];									// the image matrix, where (0, 0) is the top left corner of the rendering screen (from the observer side).
double color_table[400][3];
int divide = 1;
double grid_size = 100;
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
		is >> pt->radiance[i];
	}
	is >> pt->scaling_factor;
	is.close();
	cout << "spectrum file " << i << " read" << endl;
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
		is >> temp;
		is >> color_table[i][0];
		is >> color_table[i][1];
		is >> color_table[i][2];
	}
	is.close();
	cout << "color table read" << endl;
}
spectrum_1step convert_spectrum(Vector3 pos, spectrum_1step& sp)
{
	spectrum_1step res1;
	memset(res1.radiance, 0, sizeof(res1.radiance));
	res1.scaling_factor = sp.scaling_factor;
	double headlight_coeff = (Beta * C * global_time) * (Beta * C * global_time) / (gamma * Beta * C * global_time - pos.x / gamma) / (gamma * Beta * C * global_time - pos.x / gamma) / gamma;
	headlight_coeff /= divide;
	for (int i = 0; i < 400; i++)
	{
		double freq_orig = C / (((double)i + 380.0) / 1e9);
		double freq_trans = freq_orig * (Beta * C * global_time - pos.x) / (gamma * Beta * C * global_time - pos.x / gamma);
		double wl = (C / freq_trans) * 1e9 - 380;
		int wl_int = wl;
		if (wl_int < 0) continue;
		double wl_remain = wl - wl_int;
		if (wl_int >= 400)break;
		res1.radiance[wl_int] = sp.radiance[i] * headlight_coeff * wl_remain > 0 ? (1 - wl_remain) : 1;
		if (wl_int < 399 && wl_remain > 0)
		{
			res1.radiance[wl_int + 1] = sp.radiance[i] * headlight_coeff * wl_remain;
		}
	}
	return res1;
}
RGB spectrum_to_rgb(spectrum_1step& sp, color_system& cs)
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

	Vector3 res = { 0,0,0 }; // here the x/y/z corresponds to R/G/B.
	for (int i = 0; i < 400; i++)
	{
		double X, Y, Z;
		X = sp.radiance[i] * color_table[i][0];
		Y = sp.radiance[i] * color_table[i][1];
		Z = sp.radiance[i] * color_table[i][2];

		res.x += rvxr * X + rvyr * Y + rvzr * Z;
		res.y += rvxg * X + rvyg * Y + rvzg * Z;
		res.z += rvxb * X + rvyb * Y + rvzb * Z;

	}


	//compensate the values that cannot be represented in the current color systems (indicated as negative values) with white.
#define Min(a, b) a < b ? a : b
	double white_compensate = Min(0, Min(res.x, Min(res.y, res.z)));
	white_compensate *= -1;
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

#undef Min
#undef Max

	return { (unsigned char)res.x, (unsigned char)res.y, (unsigned char)res.z };

}

inline double length(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}

Vector3 get_hitpoint(Vector3& direction) // get the location of the point where the ray intersects with the transformed object
{
	if (!direction.is_normalized())direction.normalize();
	double lambda = (double)((double)C * (double)global_time) / (double)(1.0 + (direction.x / Beta));// coefficient for the direction vector. formula from paper
	return direction * lambda;
}
Vector3 get_original_point(Vector3& hit_pos)
{
	Vector3 res = { 0,0,0 };
	res.x = -gamma * Beta * (C * global_time - hit_pos.mod()) + gamma * hit_pos.x;
	res.y = hit_pos.y;
	res.z = hit_pos.z;
	return res;
}

void sample()
{
	cout << "starting the sampling process..." << endl;
	time_t st = 0, c = 0;
	st = clock();

	// the sampling process, where right is positive y and up is positive z.

	Vector3 ray_direction = { 0,0,0 }, target_pos = { 0,0,0 };

	int start_y = -(RES_WIDTH >> 1);
	int start_z = (RES_HEIGHT >> 1);
	int end_y = RES_WIDTH + start_y;
	int end_z = -RES_HEIGHT + start_z;

	int cnt = 0, cnt_grid = 0;

	for (int i = start_y; i <= end_y; i++) // from -0.5 to 0.5
	{
		for (int j = start_z; j >= end_z; j--) // from 0.5 to -0.5
		{
			cnt++;
			c = clock();
			time_t del = c - st;
			double tp = del;
			tp /= ((double)cnt / ((double)RES_WIDTH * (double)RES_HEIGHT));
			tp -= del;
			del = tp;
			del /= 1000;
			cout << "total progress " << (int)((double)cnt / ((double)RES_WIDTH * (double)RES_HEIGHT) * 100.0) << "%, " << "ETA: " << del / 60 << "min, " << del % 60 << "s, " << "sampling (" << i << ", " << j << ")" << endl;

			double posy = SCREEN_WIDTH * ((double)i / (double)RES_WIDTH);
			double posz = SCREEN_HEIGHT * ((double)j / (double)RES_HEIGHT);
			double posx = SCREEN_DISTANCE;
			ray_direction = { posx, posy, posz };
			ray_direction.normalize();
			Vector3 hit = get_hitpoint(ray_direction);
			Vector3 original_pos = get_original_point(hit);
			Vector3 tmp2 = original_pos;
			if (grid_size > 0)
			{
				while (tmp2.y < 0)tmp2.y += TILING_FACTOR * grid_size;
				while (tmp2.z < 0)tmp2.z += TILING_FACTOR * grid_size;
				while (tmp2.y > TILING_FACTOR * grid_size)tmp2.y -= TILING_FACTOR * grid_size;
				while (tmp2.z > TILING_FACTOR * grid_size)tmp2.z -= TILING_FACTOR * grid_size;
				if (tmp2.y < 0.001 || tmp2.z < 0.001)
				{
					res_matrix[i + (RES_WIDTH >> 1)][(RES_HEIGHT >> 1) - j] = { 255, 255, 255 };
					cnt_grid++;
					continue;
				}
			}
			int yp = original_pos.y / TILING_FACTOR;
			int zp = original_pos.z / TILING_FACTOR;
			spectrum_1step temp;
			if ((abs(yp + zp) & 1)/* || ((yp + zp) % 2 == -1)*/)
			{
				temp = convert_spectrum(hit, spec1);
			}
			else
			{
				temp = convert_spectrum(hit, spec2);
			}
			res_matrix[i + (RES_WIDTH >> 1)][(RES_HEIGHT >> 1) - j] = spectrum_to_rgb(temp, SMPTE_SYSTEM);
		}
	}
	cout << cnt_grid << endl;
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
	int cnt = 0;
	while (1)
	{
		memset(res_matrix, 0, sizeof(res_matrix));
		cout << "please input the resolution X and Y" << endl;
		cin >> RES_WIDTH >> RES_HEIGHT;
		cout << "please input the relative speed of the plane to the speed of light" << endl;
		cin >> speed_coeff;
		cout << "please input the global time" << endl;
		cin >> global_time;
		cout << "please input the divide factor of the headlight effect" << endl;
		cin >> divide;
		cout << "please input the grid size" << endl;
		cin >> grid_size;
		gamma = pow((1.0 - Beta * Beta), -0.5);
		CImg<unsigned char> res_img(RES_WIDTH, RES_HEIGHT, 1, 3, 0);
		sample();
		try
		{
			output_image(res_img);
		}
		catch (CImgException& e)
		{
			cerr << e.what() << endl;
			exit(1);
		}
		CImgDisplay disp(res_img, "result");
		while (!disp.is_closed())disp.wait();
		res_img.save("result.bmp");
	}
	//cout << "saving image to \"result.bmp\"" << endl;
	//res_img.save("result.bmp");
	//cout << "image saved to \"result.bmp\"";
	return 0;
}