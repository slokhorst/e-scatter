/**
 * @file src/sample-viewer/sample_viewer.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <algorithm>
#include <iostream>
#include <sample-viewer/sample_viewer.h>
#include <sample-viewer/shader.h>

glm::vec3 conv(const point3& v) { return glm::vec3(v.x, v.y, v.z); }

point3 mati2color(int8_t mati) {
	if(mati == 0) // Si (light grey)
		return point3(0.9,0.9,0.9);

	if(mati == -128) // 
		return point3(0.0,0.0,0.0);
	if(mati == -127) // TERMINATE (RED)
		return point3(1.0,0.0,0.0);
	if(mati == -126) // 
		return point3(0.0,0.0,0.0);
	if(mati == -125) // SE DETECTOR (GREEN)
		return point3(0.0,1.0,0.0);
	if(mati == -124) // BSE DETECTOR (BLUE)
		return point3(0.0,0.0,1.0);
	if(mati == -123) // VACUUM (GREY)
		return point3(0.3,0.3,0.3);
	if(mati == -122) // MIRROR (YELLOW)
		return point3(0.5,0.5,0.0);

	return point3(0.1,0.1,0.1);
};

sample_viewer::sample_viewer(const std::vector<triangle>* tri_vec_p)
	: m_control(), m_camera(m_control)
{
	load(tri_vec_p);
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

void sample_viewer::load(const std::vector<triangle>* tri_vec_p) {
	m_tri_vec_p = tri_vec_p;
	std::clog << "Loaded sample with " << m_tri_vec_p->size() << " interfaces" << std::endl;
}
void sample_viewer::update_voxel_mesh() {
	vertices.clear();
	normals.clear();
	colors.clear();

	point3 max;
	for(const triangle& tri : *m_tri_vec_p) {
		point3 tri_normal = tri.normal()/tri.normal().norm();
		for(const point3& v : {tri.A, tri.B, tri.C}) {
			vertices.push_back(v.x);
			vertices.push_back(v.y);
			vertices.push_back(v.z);
			normals.push_back(tri_normal.x);
			normals.push_back(tri_normal.y);
			normals.push_back(tri_normal.z);
			colors.push_back(mati2color(tri.in).x);
			colors.push_back(mati2color(tri.in).y);
			colors.push_back(mati2color(tri.in).z);

			max.x = std::max({max.x, v.x, -(v.x)});
			max.y = std::max({max.y, v.y, -(v.y)});
			max.z = std::max({max.z, v.z, -(v.z)});
		}
		for(const point3& v : {tri.C, tri.B, tri.A}) {
			vertices.push_back(v.x);
			vertices.push_back(v.y);
			vertices.push_back(v.z);
			normals.push_back(-tri_normal.x);
			normals.push_back(-tri_normal.y);
			normals.push_back(-tri_normal.z);
			colors.push_back(mati2color(tri.out).x);
			colors.push_back(mati2color(tri.out).y);
			colors.push_back(mati2color(tri.out).z);
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
	SDL_Event e;
	while( !m_control.quit_flag )
	{
		double t_prev = t;
		t = SDL_GetTicks()/1000.0f;
		double dt = t - t_prev;

		while( SDL_PollEvent(&e) != 0 )
			m_control.handleEvent(e);
		m_camera.update(dt);

		if (m_control.reload_flag) {
			load(m_tri_vec_p);
			update_voxel_mesh();
			m_control.reload_flag = false;
		}

		View = glm::lookAt(
			conv(m_camera.position()),
			conv(m_camera.direction()+m_camera.position()),
			conv(point3(0,0,1))
		);

		MVP = Projection * View * Model;

		render();
		SDL_GL_SwapWindow(gWindow);
	}
	std::clog << std::endl;

	close();
}
