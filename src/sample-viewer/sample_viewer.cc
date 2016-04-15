/**
 * @file src/sample-viewer/sample_viewer.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <algorithm>
#include <iostream>
#include <sample-viewer/sample_viewer.h>
#include <sample-viewer/shader.h>

glm::vec3 conv(const cpl::vector3& v) { return glm::vec3(v.x, v.y, v.z); }

cpl::vector3 mati2color(int8_t mati) {
	if(mati == 0) // Si (light grey)
		return cpl::vector3(0.9,0.9,0.9);

	if(mati == -128) // 
		return cpl::vector3(0.0,0.0,0.0);
	if(mati == -127) // TERMINATE (RED)
		return cpl::vector3(1.0,0.0,0.0);
	if(mati == -126) // 
		return cpl::vector3(0.0,0.0,0.0);
	if(mati == -125) // SE DETECTOR (GREEN)
		return cpl::vector3(0.0,1.0,0.0);
	if(mati == -124) // BSE DETECTOR (BLUE)
		return cpl::vector3(0.0,0.0,1.0);
	if(mati == -123) // VACUUM (GREY)
		return cpl::vector3(0.3,0.3,0.3);
	if(mati == -122) // MIRROR (YELLOW)
		return cpl::vector3(0.5,0.5,0.0);

	return cpl::vector3(0.1,0.1,0.1);
};

sample_viewer::sample_viewer(const std::vector<material_interface>* mi_vec_p)
	: m_control(), m_camera(m_control)
{
	load(mi_vec_p);
	run();
}

void sample_viewer::init() {
	if(SDL_Init(SDL_INIT_VIDEO) < 0)
		throw std::runtime_error("SDL could not initialize! SDL Error: "+std::string(SDL_GetError()));

	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);

	gWindow = SDL_CreateWindow("Sample viewer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
	if(gWindow == NULL)
		throw std::runtime_error("Window could not be created! SDL Error: "+std::string(SDL_GetError()));

	gContext = SDL_GL_CreateContext(gWindow);
	if(gContext == NULL)
		throw std::runtime_error("OpenGL context could not be created! SDL Error: "+std::string(SDL_GetError()));

	if(SDL_GL_SetSwapInterval(1) < 0)
		std::clog << "Unable to set VSync! SDL Error: " << std::string(SDL_GetError()) << std::endl;

	glewExperimental = GL_TRUE;
	if(glewInit() != GLEW_OK)
		throw std::runtime_error("Failed to initialize GLEW");

	initGL();
}
void sample_viewer::close() {
	termGL();
	SDL_DestroyWindow( gWindow );
	SDL_Quit();
}

void sample_viewer::initGL() {
	glClearColor(0.1f,0.1f,0.1f,1.0f);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glEnable(GL_CULL_FACE);

	shaderProgramID = loadShaders("src/sample-viewer/shader.vertex", "src/sample-viewer/shader.fragment");
	MVPmatrixID = glGetUniformLocation(shaderProgramID, "MVP");
	MmatrixID = glGetUniformLocation(shaderProgramID, "M");
	VmatrixID = glGetUniformLocation(shaderProgramID, "V");
	PmatrixID = glGetUniformLocation(shaderProgramID, "P");

	Model = glm::mat4(1.0f);
	View = glm::mat4(1.0f);
	Projection = glm::perspective(45.0f, (float)SCREEN_WIDTH/SCREEN_HEIGHT, 0.001f, 1000.0f);

	glGenBuffers(1, &points_vbo);
	glGenBuffers(1, &normals_vbo);
	glGenBuffers(1, &colors_vbo);
	glGenVertexArrays(1, &vao);
}
void sample_viewer::termGL() {
	glDeleteBuffers(1, &points_vbo);
	glDeleteBuffers(1, &normals_vbo);
	glDeleteBuffers(1, &colors_vbo);
	glDeleteVertexArrays(1, &vao);
}

void sample_viewer::load(const std::vector<material_interface>* mi_vec_p) {
	m_mi_vec_p = mi_vec_p;
	std::clog << "Loaded sample with " << m_mi_vec_p->size() << " interfaces" << std::endl;
}
void sample_viewer::update_voxel_mesh() {
	vertices.clear();
	normals.clear();
	colors.clear();

	cpl::vector3 max;
	for(const material_interface mi : *m_mi_vec_p) {
		cpl::vector3 tri_normal = mi.triangle.normal()/mi.triangle.normal().norm();
		for(auto i = mi.triangle.vertices().cbegin(); i != mi.triangle.vertices().cend(); ++i) {
			vertices.push_back(i->x);
			vertices.push_back(i->y);
			vertices.push_back(i->z);
			normals.push_back(tri_normal.x);
			normals.push_back(tri_normal.y);
			normals.push_back(tri_normal.z);
			colors.push_back(mati2color(mi.in_mat).x);
			colors.push_back(mati2color(mi.in_mat).y);
			colors.push_back(mati2color(mi.in_mat).z);

			max.x = std::max({max.x, i->x, -(i->x)});
			max.y = std::max({max.y, i->y, -(i->y)});
			max.z = std::max({max.z, i->z, -(i->z)});
		}
		for(auto i = mi.triangle.vertices().crbegin(); i != mi.triangle.vertices().crend(); ++i) {
			vertices.push_back(i->x);
			vertices.push_back(i->y);
			vertices.push_back(i->z);
			normals.push_back(-tri_normal.x);
			normals.push_back(-tri_normal.y);
			normals.push_back(-tri_normal.z);
			colors.push_back(mati2color(mi.out_mat).x);
			colors.push_back(mati2color(mi.out_mat).y);
			colors.push_back(mati2color(mi.out_mat).z);
		}
	}

	n_vertices = vertices.size()/3;
	double max_dim = std::max({max.x,max.y,max.z});

	std::clog << "Generated mesh with " << (n_vertices/3) << " triangles" << std::endl;
	std::clog << "Dimensions: (" << max.x << "," << max.y << "," << max.z << "). Scaling by 0.5/" << max_dim << "." << std::endl;

	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
	glBufferData(GL_ARRAY_BUFFER, normals.size()*sizeof(float), normals.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
	glBufferData(GL_ARRAY_BUFFER, colors.size()*sizeof(float), colors.data(), GL_STATIC_DRAW);

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, NULL);


	Model = glm::scale(glm::mat4(1.0f),glm::vec3(0.5/max_dim));
}

void sample_viewer::render() const {
	if(m_control.wireframe_flag)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	else
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

	glUniformMatrix4fv(MmatrixID, 1, GL_FALSE, &Model[0][0]);
	glUniformMatrix4fv(VmatrixID, 1, GL_FALSE, &View[0][0]);
	glUniformMatrix4fv(PmatrixID, 1, GL_FALSE, &Projection[0][0]);
	glUniformMatrix4fv(MVPmatrixID, 1, GL_FALSE, &MVP[0][0]);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);

	glUseProgram(shaderProgramID);

	glBindVertexArray(vao);
	glDrawArrays(GL_TRIANGLES, 0, n_vertices);

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
}

void sample_viewer::run() {
	try {
		init();
	} catch (const std::exception& e) {
		close();
		std::clog << "Failed to initialize: " << e.what() << std::endl;
		return;
	}

	update_voxel_mesh();

	double t = SDL_GetTicks()/1000.0f;
	double t_prev = t;
	double dt = 0;
	SDL_Event e;
	while( !m_control.quit_flag )
	{
		t_prev = t;
		t = SDL_GetTicks()/1000.0f;
		dt = t - t_prev;

		while( SDL_PollEvent(&e) != 0 )
			m_control.handleEvent(e);
		m_camera.update(dt);

		if (m_control.reload_flag) {
			load(m_mi_vec_p);
			update_voxel_mesh();
			m_control.reload_flag = false;
		}

		View = glm::lookAt(
			conv(m_camera.position()),
			conv(m_camera.direction()+m_camera.position()),
			conv(cpl::vector3(0,0,1))
		);

		MVP = Projection * View * Model;

		render();
		SDL_GL_SwapWindow(gWindow);
	}
	std::clog << std::endl;

	close();
}
