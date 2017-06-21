import array
import os.path
import math
import warnings
from .geometry import Triangle

from PyQt5.QtGui import (
    QOpenGLShader,
    QOpenGLShaderProgram,
    QOpenGLVersionProfile,
    QSurfaceFormat,
    QMatrix4x4,
    QVector3D
)
from PyQt5.QtWidgets import QOpenGLWidget


mat_to_color_map = {
    0: [0.9, 0.9, 0.9],  # light grey
    1: [0.9, 0.5, 0.5],  # light red
    2: [0.5, 0.9, 0.5],  # light green
    3: [0.5, 0.5, 0.9],  # light blue
 -128: [0.0, 0.0, 0.0],  # NOP (black)
 -127: [1.0, 0.0, 0.0],  # TERMINATE (red)
 -126: [0.0, 1.0, 1.0],  # DETECTOR (cyan)
 -125: [0.0, 1.0, 0.0],  # SE DETECTOR (green)
 -124: [0.0, 0.0, 1.0],  # BSE DETECTOR (blue)
 -123: [1.0, 1.0, 1.0],  # VACUUM (white)
 -122: [0.5, 0.5, 0.0],  # MIRROR (yellow)
}


def mat_to_color(mati):
    if mati in mat_to_color_map:
        return array.array('f', mat_to_color_map[mati])
    else:
        return array.array('f', [0.1, 0.1, 0.1])


def parseFile(filename):
    triangles = []
    with open(filename, 'r') as f, warnings.catch_warnings():
        warnings.filterwarnings('error')
        ln = 0
        for line in f.readlines():
            try:
                ln += 1
                line = line.strip()
                if line == '' or line[0] == '#':
                    continue
                parts = line.split()
                if len(parts) != 11:
                    raise RuntimeError('malformed line')

                t = Triangle(
                    [float(parts[2]), float(parts[3]), float(parts[4])],
                    [float(parts[5]), float(parts[6]), float(parts[7])],
                    [float(parts[8]), float(parts[9]), float(parts[10])],
                    int(parts[0]), int(parts[1])
                )
                triangles.append(t)
            except Warning as e:
                print('warning: line {}: {}'.format(ln, e))
            except Exception as e:
                print('error: line {}: {}'.format(ln, e))
                exit()
    return triangles


class SampleViewerWidget(QOpenGLWidget):
    def __init__(self):
        super().__init__()
        self.num_triangles = 0
        self.wireframe = False
        self.max_dim = 1
        self.lightDir = QVector3D(0, 0, 1)
        self.resetCamera()

    def initializeGL(self):
        vp = QOpenGLVersionProfile()
        vp.setVersion(4, 1)  # PyQt supports up to OpenGL 4.1
        vp.setProfile(QSurfaceFormat.CoreProfile)
        self.gl = self.context().versionFunctions(vp)
        if not self.gl:
            raise RuntimeError("unable to set OpenGL version profile")

        self.gl.initializeOpenGLFunctions()
        self.gl.glClearColor(0.1, 0.1, 0.1, 1.0)
        self.gl.glEnable(self.gl.GL_DEPTH_TEST)
        self.gl.glDepthFunc(self.gl.GL_LESS)
        self.gl.glEnable(self.gl.GL_CULL_FACE)

        with open(os.path.join(os.path.dirname(__file__), 'vertex.glsl')) as f:
            vertexShaderSource = f.read()

        with open(os.path.join(os.path.dirname(__file__), 'fragment.glsl')) as f:
            fragmentShaderSource = f.read()

        self.createShaders(vertexShaderSource, fragmentShaderSource)
        self.loadTriangles([])

    def createShaders(self, vertexCode, fragmentCode):
        self.shader_prog = QOpenGLShaderProgram(self)

        self.shader_prog.addShaderFromSourceCode(QOpenGLShader.Vertex, vertexCode)
        self.shader_prog.addShaderFromSourceCode(QOpenGLShader.Fragment, fragmentCode)
        self.shader_prog.link()

        self.attr_vPosition = self.shader_prog.attributeLocation('vertexPosition_ms')
        self.attr_vNormal = self.shader_prog.attributeLocation('vertexNormal_ms')
        self.attr_vColor = self.shader_prog.attributeLocation('vertexColor')
        self.unif_MVP = self.shader_prog.uniformLocation('MVP')
        self.unif_lightDir = self.shader_prog.uniformLocation('lightDirection_ms')

    def setData(self, vertexArray, normalArray, colorArray):
        if len(vertexArray) != len(colorArray) or len(vertexArray) != len(normalArray):
            raise RuntimeError("vertexArray, normalArray and colorArray must be same length")
        self.num_triangles = int(len(vertexArray)/3/3)

        self.vertexArray = vertexArray
        self.colorArray = colorArray
        self.normalArray = normalArray

    def loadTriangles(self, triangles):
        vertexArray = array.array('f')
        normalArray = array.array('f')
        colorArray = array.array('f')
        for t in triangles:
            color_i = mat_to_color(t.mat_i)
            color_o = mat_to_color(t.mat_o)
            for v in t.v1, t.v2, t.v3:
                vertexArray += array.array('f', v)
                normalArray += array.array('f', t.n)
                colorArray += color_i
                self.max_dim = max(self.max_dim, v[0], -v[0])
                self.max_dim = max(self.max_dim, v[1], -v[1])
                self.max_dim = max(self.max_dim, v[2], -v[2])
            for v in t.v3, t.v2, t.v1:
                vertexArray += array.array('f', v)
                normalArray += array.array('f', -t.n)
                colorArray += color_o
        self.setData(vertexArray, normalArray, colorArray)

    def loadFile(self, filename):
        print("loading {}".format(filename))
        triangles = parseFile(filename)
        print("loaded {} triangles".format(len(triangles)))
        self.loadTriangles(triangles)
        self.update()

    def paintGL(self):
        if self.wireframe:
            self.gl.glPolygonMode(self.gl.GL_FRONT_AND_BACK, self.gl.GL_LINE)
        else:
            self.gl.glPolygonMode(self.gl.GL_FRONT_AND_BACK, self.gl.GL_FILL)

        self.gl.glClear(self.gl.GL_COLOR_BUFFER_BIT | self.gl.GL_DEPTH_BUFFER_BIT)

        self.shader_prog.bind()

        camera_pos = QVector3D(
            self.camera_pos_sph.x() * math.sin(self.camera_pos_sph.y())
                                    * math.cos(self.camera_pos_sph.z()),
            self.camera_pos_sph.x() * math.sin(self.camera_pos_sph.y())
                                    * math.sin(self.camera_pos_sph.z()),
            self.camera_pos_sph.x() * math.cos(self.camera_pos_sph.y())
        )

        MVP = QMatrix4x4()
        MVP.perspective(45, self.w/self.h, 1e-3, 1e4)
        MVP.lookAt(camera_pos, QVector3D(0, 0, 0), QVector3D(0, 0, 1))
        MVP.scale(QVector3D(1, 1, 1)/self.max_dim)

        self.shader_prog.setUniformValue(self.unif_MVP, MVP)
        self.shader_prog.setUniformValue(self.unif_lightDir, self.lightDir)

        self.gl.glEnableVertexAttribArray(self.attr_vPosition)
        self.gl.glEnableVertexAttribArray(self.attr_vNormal)
        self.gl.glEnableVertexAttribArray(self.attr_vColor)

        self.gl.glVertexAttribPointer(self.attr_vPosition, 3, self.gl.GL_FLOAT, False, 0, self.vertexArray)
        self.gl.glVertexAttribPointer(self.attr_vNormal, 3, self.gl.GL_FLOAT, False, 0, self.normalArray)
        self.gl.glVertexAttribPointer(self.attr_vColor, 3, self.gl.GL_FLOAT, False, 0, self.colorArray)

        self.gl.glDrawArrays(self.gl.GL_TRIANGLES, 0, 3*self.num_triangles)

        self.gl.glDisableVertexAttribArray(self.attr_vPosition)
        self.gl.glDisableVertexAttribArray(self.attr_vNormal)
        self.gl.glDisableVertexAttribArray(self.attr_vColor)

        self.shader_prog.release()

    def resizeGL(self, w, h):
        self.w = w
        self.h = h

    def setWireframe(self, value):
        self.wireframe = bool(value)

    def resetCamera(self):
        self.camera_pos_sph = QVector3D(3, 0.7, 0.7)

    def setLighting(self, value):
        self.lighting = bool(value)
        if self.lighting:
            self.lightDir = QVector3D(0, 0, 1)
        else:
            self.lightDir = QVector3D(0, 0, 0)

    def mouseMoveEvent(self, e):
        delta = e.pos() - self.mousePos
        self.camera_pos_sph -= QVector3D(0, delta.y(), delta.x())/200
        self.camera_pos_sph.setY(min(max(self.camera_pos_sph.y(), 0.01), 3.14))
        self.update()
        self.mousePos = e.pos()

    def mousePressEvent(self, e):
        self.mousePos = e.pos()

    def mouseReleaseEvent(self, e):
        self.mousePos = None

    def wheelEvent(self, e):
        self.camera_pos_sph -= QVector3D(e.angleDelta().y(), 0, 0)/1000
        self.update()
