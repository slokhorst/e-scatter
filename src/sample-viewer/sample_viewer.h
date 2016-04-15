/**
 * @file src/sample-viewer/sample_viewer.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SAMPLE_VIEWER__SAMPLE_VIEWER__HEADER_INCLUDED
#define SAMPLE_VIEWER__SAMPLE_VIEWER__HEADER_INCLUDED

#include <vector>
#include <GL/glew.h>
#include <SDL.h>
#include <SDL_opengl.h>

//for GLM<0.9.6 (Ubuntu<15.10)
#define GLM_FORCE_RADIANS

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <sample-viewer/control.h>
#include <sample-viewer/camera.h>

#include <cpl/vector3.h>
#include <cpl/triangle.h>

class material_interface {
public:
	material_interface(cpl::triangle t, int8_t i_m, int8_t o_m)
		: triangle(t), in_mat(i_m), out_mat(o_m) {};
	material_interface()
		: triangle({cpl::vector3(),cpl::vector3(),cpl::vector3()}), in_mat(0), out_mat(1) {};
	cpl::triangle triangle;
	int8_t in_mat, out_mat;
};


class sample_viewer {
public:
	sample_viewer(const std::vector<material_interface>* mi_vec_p);

private:
	const std::vector<material_interface>* m_mi_vec_p;
	std::vector<float> vertices;
	std::vector<float> normals;
	std::vector<float> colors;

	const int SCREEN_WIDTH = 1280;
	const int SCREEN_HEIGHT = 720;

	SDL_Window* gWindow;
	SDL_GLContext gContext;

	GLuint points_vbo = 0, normals_vbo = 0, colors_vbo = 0;
	GLuint vao = 0;
	GLuint shaderProgramID;
	GLuint MmatrixID, VmatrixID, PmatrixID, MVPmatrixID;

	glm::mat4 Model, View, Projection, MVP;

	int n_vertices;

	control m_control;
	camera m_camera;

	void load(const std::vector<material_interface>* mi_vec_p);

	void init();
	void initGL();
	void termGL();
	void close();
	void handleEvent(const SDL_Event& e);
	void render() const;
	void update_voxel_mesh();

	void run();
};

#endif
