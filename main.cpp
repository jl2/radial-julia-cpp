#include <iostream>
#include <string>
#include <algorithm>
#include <complex>

#include <png++/png.hpp>

#include <random>
#include <cmath>

#include <omp.h>


typedef std::complex<double> dcomp;

struct Pixel {
	double red;
	double green;
	double blue;
	double alpha;
	Pixel(double r=0.0, double g=0.0, double b=0.0, double a=1.0) : red(r), green(g), blue(b), alpha(a) {}
};

class Image {
public:
	Image(unsigned int width, unsigned int height) : _w(width), _h(height), _pixels(_h, std::vector<Pixel>(_w)) {
	}

	void setPixel(unsigned int x, unsigned int y, const Pixel& px) {
		_pixels[y][x] = px;
	}
	Pixel getPixel(unsigned int x, unsigned int y) {
		return _pixels[y][x];
	}
	void ToFile(const std::string& fname) {
		png::image< png::rgba_pixel > image(_w, _h);

		for (size_t y = 0; y < _h; ++y) {
			for (size_t x = 0; x < _w; ++x) {
				const Pixel val = getPixel(x, y);
				image[y][x] = png::rgba_pixel(int(val.red*255), int(val.green*255), int(val.blue*255), int(val.alpha*255));
			}
		}
		image.write(fname.c_str());
	}
private:
	unsigned int _w, _h;
	std::vector< std::vector<Pixel> > _pixels;
};

Pixel blackAndWhite(size_t iters, size_t iterations, size_t i, size_t j, size_t width, size_t height)
{
    size_t val = iters & 1;
    return Pixel(val, val, val, 1.0);
}

double mapValue(int x, size_t width, double xmin, double xmax)
{
    return xmin + (xmax - xmin) * ((double(x) / width));
}


void radialJuliaLine(size_t j,
                     dcomp c,
                     Image &img,
                     size_t width = 800,
                     size_t height = 800,
                     size_t max_iterations=120,
                     double rmin=0.0, double rmax=2 * 3.141592654,
                     double tmin=0.0, double tmax=1.25)
{
    double rp = mapValue(j, height, rmin, rmax);
    for (size_t i = 0; i<width; ++i) {

        
        double tp = mapValue(i, width, tmin, tmax);
        dcomp pt = std::polar(rp, tp);

        int iteration = 0;
        while (std::norm(pt) < 16 && iteration<max_iterations) {
            pt = pt * pt + c;
            ++iteration;
        }
        Pixel tmp = blackAndWhite(iteration, max_iterations, i, j, width, height);
        img.setPixel(i, j, tmp);
    }
}


void radialJulia(std::string file_name,
                 dcomp c,
                 size_t width,
                 size_t height,
                 size_t max_iterations=100,
                 double rmin=0.0, double rmax=1.25,
                 double tmin=0.0, double tmax=2*3.141592654,
                 size_t thread_count = 4)
{
    Image img(width, height);
#pragma omp parallel for
    for (size_t i = 0; i < width; ++i)
    {
        radialJuliaLine(i, c, img, width, height, max_iterations, rmin, rmax, tmin, tmax);
    }
    img.ToFile(file_name);
}

dcomp randomComplex(const dcomp min_val, const dcomp max_val)
{
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    std::uniform_real_distribution<> dis(0.0, 1.0);

    dcomp diff = max_val - min_val;

    return dcomp(std::real(diff) * dis(gen), std::imag(diff) * dis(gen)) + min_val;
}

void randomWalkJuliaSet(size_t frame_count, std::string output_dir, double dc)
{
    dcomp startPt = randomComplex(dcomp(0.20, 0.25), dcomp(0.50, 0.5));
    double real_dir = 1.0;
    double imag_dir = -1.0;
    dcomp lower_bound( 0.125, -0.5);
    dcomp upper_bound(0.5, 0.5);

    const dcomp min_val = dcomp(0.0, 0.0);
    const dcomp max_val = dcomp(dc, dc);
    
    for (size_t frame = 0; frame < frame_count; ++frame)
    {
        std::cout << "Rendering " << startPt << "\n";
        char buffer[512]={0};
        snprintf(buffer, 512, "%s/juliaset%06ul.png", output_dir.c_str(), frame);
        radialJulia(buffer,
                    startPt,
                    1200, 1200, 100,
                    0.0, 1.25,
                    0.0, 2 * 3.141592654,
                    4);
        dcomp dc = dcomp(real_dir, imag_dir);
        dcomp randVal = randomComplex(min_val, max_val);
        std::cout << "dcomp(real_dir, imag_dir) *  = " << dc << " " << randVal << " = " << (dc * randVal) << "\n";
        startPt += dc * randVal;
        if (std::real(startPt) > std::real(upper_bound))
        {
            real_dir = -real_dir;
        }
        if (std::imag(startPt) > std::imag(upper_bound))
        {
            imag_dir = -imag_dir;
        }

        if (std::real(startPt) < std::real(lower_bound))
        {
            real_dir = -real_dir;
        }
        if (std::imag(startPt) < std::imag(lower_bound))
        {
            imag_dir = -imag_dir;
        }

    }
}

int main (int argc, char *argv[])
{
    randomWalkJuliaSet(1200, argv[1], 0.001);
    return 0;
}
