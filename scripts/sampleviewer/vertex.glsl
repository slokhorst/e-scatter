#version 400

in vec3 vertexPosition_ms;
in vec3 vertexNormal_ms;
in vec3 vertexColor;

out vec3 fragmentColor;
out float lightIntensity;

uniform mat4 MVP;
uniform vec3 lightDirection_ms;

void main() {
	gl_Position = MVP * vec4(vertexPosition_ms,1);

	fragmentColor = vertexColor;

	if(lightDirection_ms.x == 0 && lightDirection_ms.y == 0 && lightDirection_ms.z == 0)
		lightIntensity = 1.0;
	else
		lightIntensity = clamp(dot(vertexNormal_ms,lightDirection_ms),0,1);
}
