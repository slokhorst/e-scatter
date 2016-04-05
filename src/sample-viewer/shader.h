/**
 * @file src/sample-viewer/shader.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SAMPLE_VIEWER__SHADER__HEADER_INCLUDED
#define SAMPLE_VIEWER__SHADER__HEADER_INCLUDED

GLuint loadShaders(const std::string& vertex_file_path, const std::string& fragment_file_path){
	std::string vertexCode = cpl::text::file_to_string(vertex_file_path);
	std::string fragmentCode = cpl::text::file_to_string(fragment_file_path);

	GLint result = GL_FALSE;
	GLint infoLogLength;

	GLuint vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	const GLchar* vertexCode_p = vertexCode.c_str();
	glShaderSource(vertexShaderID, 1, &vertexCode_p, NULL);
	glCompileShader(vertexShaderID);
	glGetShaderiv(vertexShaderID, GL_COMPILE_STATUS, &result);
	if(result != GL_TRUE) {
		glGetShaderiv(vertexShaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
		std::string infoLogString(infoLogLength, ' ');
		glGetShaderInfoLog(vertexShaderID, infoLogLength, NULL, &infoLogString[0]);
		throw std::runtime_error(cpl::text::cat("error compiling vertex shader:\n",infoLogString));
	}

	GLuint fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	const GLchar* fragmentCode_p = fragmentCode.c_str();
	glShaderSource(fragmentShaderID, 1, &fragmentCode_p, NULL);
	glCompileShader(fragmentShaderID);
	glGetShaderiv(fragmentShaderID, GL_COMPILE_STATUS, &result);
	if(result != GL_TRUE) {
		glGetShaderiv(fragmentShaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
		std::string infoLogString(infoLogLength, ' ');
		glGetShaderInfoLog(fragmentShaderID, infoLogLength, NULL, &infoLogString[0]);
		throw std::runtime_error(cpl::text::cat("error compiling fragment shader:\n",infoLogString));
	}

	GLuint programID = glCreateProgram();
	glAttachShader(programID, vertexShaderID);
	glAttachShader(programID, fragmentShaderID);
	glLinkProgram(programID);
	glGetProgramiv(programID, GL_LINK_STATUS, &result);
	if(result != GL_TRUE) {
		glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &infoLogLength);
		std::string infoLogString(infoLogLength, ' ');
		glGetProgramInfoLog(programID, infoLogLength, NULL, &infoLogString[0]);
		throw std::runtime_error(cpl::text::cat("error linking shader program:\n",infoLogString));
	}
	glDeleteShader(vertexShaderID);
	glDeleteShader(fragmentShaderID);

	return programID;
}

#endif
