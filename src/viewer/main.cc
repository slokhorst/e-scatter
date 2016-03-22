/**
 * @file src/viewer/viewer.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <cstdlib>
#include <iostream>
#include <map>
#include <signal.h>
#include <string>
#include <SDL.h>

class image {
public:
	image(int width, int height) {
		_allocate(width,height);
	}
	image(const image& img) {
		*this = img;
	}
	~image() {
		close_window();
		delete[] _pixel_p;
	}
	image& operator=(const image& img) {
		if(this != &img) {
			_allocate(img._width,img._height);
			for(int y = 0; y < _height; y++)
			for(int x = 0; x < _width; x++)
				(*this)(x,y) = img(x,y);
		}
		return *this;
	}
	int width() const {
		return _width;
	}
	int height() const {
		return _height;
	}
	float min() const {
		float i = _pixel_p[0];
		for(int y = 0; y < _height; y++)
		for(int x = 0; x < _width; x++)
			i = std::min(i, (*this)(x,y));
		return i;
	}
	float max() const {
		float i = _pixel_p[0];
		for(int y = 0; y < _height; y++)
		for(int x = 0; x < _width; x++)
			i = std::max(i, (*this)(x,y));
		return i;
	}
	const float& operator()(int x, int y) const {
		x = std::max(0, std::min(x, _width-1));
		y = std::max(0, std::min(y, _height-1));
		return _pixel_p[x+y*_width];
	}
	float& operator()(int x, int y) {
		x = std::max(0, std::min(x, _width-1));
		y = std::max(0, std::min(y, _height-1));
		return _pixel_p[x+y*_width];
	}
	void display_window() {
		display_window(_width,_height);
	}
	void display_window(int width, int height) {
		if(_display_map.empty()) {
			SDL_Init(SDL_INIT_NOPARACHUTE|SDL_INIT_VIDEO);
			signal(SIGINT,SIG_DFL);
		}
		if(_window_p == nullptr) {
			SDL_CreateWindowAndRenderer(width, height, SDL_WINDOW_SHOWN, &_window_p, &_renderer_p);
			_display_map[this] = _window_p;
		} else
			SDL_SetWindowSize(_window_p, width, height);
		SDL_Surface* surface_p = SDL_CreateRGBSurface(0,_width, _height, 32, 0, 0, 0, 0);
		const float min_i = min();
		const float max_i = max();
		for(int y = 0; y < _height; y++)
		for(int x = 0; x < _width; x++) {
			float i = (*this)(x, y);
			if(i > 0)
				i -= min_i;
			if(max_i != min_i)
				i /= max_i-min_i;
			Uint8 alpha = 255;
			Uint8 gray = 255*i;
			Uint8* pixel_p = reinterpret_cast<Uint8*>(surface_p->pixels);
			pixel_p += x*sizeof(Uint32)+y*surface_p->pitch;
			*reinterpret_cast<Uint32*>(pixel_p) = (alpha << 24) | (gray << 16) | (gray << 8) | gray;
		}
		SDL_Texture* texture_p = SDL_CreateTextureFromSurface(_renderer_p, surface_p);
		SDL_RenderCopy(_renderer_p, texture_p, nullptr, nullptr);
		SDL_RenderPresent(_renderer_p);
		SDL_DestroyTexture(texture_p);
		SDL_FreeSurface(surface_p);
	}
	void close_window() {
		if(_renderer_p != nullptr) {
			SDL_DestroyRenderer(_renderer_p);
			_renderer_p = nullptr;
		}
		if(_window_p != nullptr) {
			_display_map.erase(this);
			SDL_DestroyWindow(_window_p);
			_window_p = nullptr;
			if(_display_map.empty())
				SDL_Quit();
		}
	}
	bool has_window() const {
		if(_window_p != nullptr)
			return true;
		return false;
	}
	static int window_count() {
		return _display_map.size();
	}
	static void window_poll() {
		if(_display_map.empty())
			return;
		SDL_Event event;
		while(SDL_PollEvent(&event))
			if(event.window.event == SDL_WINDOWEVENT_CLOSE)
				for(auto cit = _display_map.cbegin(); cit != _display_map.cend(); cit++)
					if(event.window.windowID == SDL_GetWindowID(cit->second)) {
						cit->first->close_window();
						break;
					}
	}
private:
	void _allocate(int width, int height) {
		delete[] _pixel_p;
		_pixel_p = new float[width*height];
		_width = width;
		_height = height;
		for(int y = 0; y < _height; y++)
		for(int x = 0; x < _width; x++)
			(*this)(x,y) = 0;
	}
	float* _pixel_p = nullptr;
	int _width;
	int _height;
	SDL_Window* _window_p = nullptr;
	SDL_Renderer* _renderer_p = nullptr;
	static std::map<image*,SDL_Window*> _display_map;
};

std::map<image*,SDL_Window*> image::_display_map;

int main(const int argc, char* argv[]) {
	image sem_image(128, 1024);
	sem_image.display_window();
	int i = 0;
	while(true) {
		float rx, ry, rz;
		std::cin.read(reinterpret_cast<char*>(&rx), sizeof(rx));
		std::cin.read(reinterpret_cast<char*>(&ry), sizeof(ry));
		std::cin.read(reinterpret_cast<char*>(&rz), sizeof(rz));
		float dx, dy, dz;
		std::cin.read(reinterpret_cast<char*>(&dx), sizeof(dx));
		std::cin.read(reinterpret_cast<char*>(&dy), sizeof(dy));
		std::cin.read(reinterpret_cast<char*>(&dz), sizeof(dz));
		float K;
		std::cin.read(reinterpret_cast<char*>(&K), sizeof(K));
		int px, py;
		std::cin.read(reinterpret_cast<char*>(&px), sizeof(px));
		std::cin.read(reinterpret_cast<char*>(&py), sizeof(py));
		if(std::cin.eof())
			break;
		if(K < 50)
			sem_image(px, py)++;
		if((i++)%10240 == 0) {
			image::window_poll();
			if(!sem_image.has_window())
				break;
			sem_image.display_window();
		}
	}
	if(sem_image.has_window()) {
		sem_image.display_window();
		while(sem_image.has_window())
			image::window_poll();
	}
	return EXIT_SUCCESS;
}