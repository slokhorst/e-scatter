/**
 * @file src/sample-viewer/control.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SAMPLE_VIEWER__CONTROL__HEADER_INCLUDED
#define SAMPLE_VIEWER__CONTROL__HEADER_INCLUDED

#include <cdsem/point3.hh>


class control
{
private:
	point3 requested_movement_dir;
	bool forward=false, backward=false, left=false, right=false, up=false, down=false;
	bool capture_mouse=false;
	std::pair<int,int> mmove = std::pair<int,int>(0,0);
	int mouse_dz=0;
public:
	bool quit_flag=false, reload_flag=false, wireframe_flag=false;
	control() {}
	void handleEvent(const SDL_Event& e) {
		switch(e.type) {
			case SDL_KEYDOWN: {
				switch(e.key.keysym.sym) {
					case SDLK_w: 
						forward = true; break;
					case SDLK_a:
						left = true; break;
					case SDLK_s:
						backward = true; break;
					case SDLK_d:
						right = true; break;

					case SDLK_q:
						quit_flag = true; break;
					case SDLK_r:
						reload_flag = true; break;
					case SDLK_e:
						wireframe_flag = !wireframe_flag; break;
				}
				break;
			}
			case SDL_KEYUP: {
				switch(e.key.keysym.sym) {
					case SDLK_w:
						forward = false; break;
					case SDLK_a:
						left = false; break;
					case SDLK_s:
						backward = false; break;
					case SDLK_d:
						right = false; break;
				}
				break;
			}
			case SDL_MOUSEMOTION: {
				mmove.first += e.motion.xrel;
				mmove.second += e.motion.yrel;
				break;
			}
			case SDL_MOUSEWHEEL: {
				mouse_dz = e.wheel.y; break;
			}
			case SDL_MOUSEBUTTONDOWN: {
				SDL_SetRelativeMouseMode(SDL_TRUE);
				capture_mouse = true; break;
			}
			case SDL_MOUSEBUTTONUP: {
				SDL_SetRelativeMouseMode(SDL_FALSE);
				capture_mouse = false; break;
			}
			case SDL_QUIT: {
				quit_flag = true; break;
			}
		}
	}
	std::pair<double,double> mouse_movement() {
		auto ret = mmove;
		mmove = std::make_pair(0,0);
		if(capture_mouse)
			return ret;
		else
			return std::make_pair(0,0);
	}
	point3 requested_movement() const {
		return point3(
			(right  ?1:0) - (left    ?1:0),
			(forward?1:0) - (backward?1:0),
			(up     ?1:0) - (down    ?1:0)
		);
	}
};

#endif
