//  Created by Antonio J. Rueda on 21/3/2024.
//  Copyright © 2024 Antonio J. Rueda.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

#include "util.h"
#include "particulas.h"

// Vertex shader para visualización de partículas
static const char* codVisVS = R"(
    #version 330
              
    layout(location = 0) in vec2 coord;
    layout(location = 1) in vec3 color;

    out vec3 colorEntrada;
                                               
    void main() {
        gl_Position = vec4(coord, 0.0, 1.0);
        gl_PointSize = 10.0;

        colorEntrada = color;
    }
)";


// Fragment shader para visualización de partículas
static const char* codVisFS = R"(
    #version 330
    in vec3 colorEntrada;                                           
    out vec4 color;
                                               
    void main() {
        // Dibujar los puntos como círculos
        vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
        if (dot(circCoord, circCoord) > 1.0) {
            discard;
        }
                                               
       color = vec4(colorEntrada, 0.75);
    }
)";

// Vertex shader para simulación de movimiento de partículas
static const char* codSimVS = R"(
    #version 330
    
    layout(location = 0) in vec2 coord;
    
    void main() {
        gl_Position = vec4(coord, 0.0, 1.0);
    }
)";


// Fragment shader para simulación de movimiento de partículas
static const char* codSimFS = R"(
    #version 330

    uniform sampler2D coordsIn;
    uniform sampler2D coordsPrevious;

    uniform float pasoT;
    uniform int numParticulas;
    uniform int anchoTextura;
    uniform float G;

    layout(location = 0) out vec2 coordOut;
    layout(location = 1) out float speedOut;

    vec2 compute(vec2 coordsActuales, vec2 coordsPrevias, sampler2D coordsIn, float pasoT, int numParticulas, int anchoTextura, float G) {
        // CALCULAMOS DISTANCIAS
        int countParticulas = 0;
        int countFila = 0;
        int countColumna = 0;
        vec2 sum = vec2(0, 0);
        while(countParticulas < numParticulas) {
            if(countColumna == anchoTextura) { 
                countColumna = 0;
                countFila; // Nos bajamos a una nueva fila
            }

            ivec2 texel = ivec2 (countColumna, countFila);
            vec2 otherParticle = texelFetch(coordsIn, texel, 0).xy; 

            if(coordsActuales != otherParticle) { 
                float d = length(otherParticle - coordsActuales);
                if(d >= 0.01f) sum += (G/(d*d*d))*(otherParticle - coordsActuales); // En caso de que la distancia salga 0 por ser un valor pequeño
            }

            ++countColumna;
            ++countParticulas;
        }

        // CALCULAMOS EL RESULTADO
        vec2 res = 2*coordsActuales - coordsPrevias + (pasoT*pasoT) * sum;

        return res;
    }

    void main() {
        ivec2 pos = ivec2(int(gl_FragCoord.x), int(gl_FragCoord.y));

        if (pos.y * anchoTextura + pos.x >= numParticulas)
            discard;

        vec2 coordsActuales = texelFetch(coordsIn, pos, 0).xy;
        vec2 coordsPrevias = texelFetch(coordsPrevious, pos, 0).xy;

        vec2 result = compute(coordsActuales, coordsPrevias, coordsIn, pasoT, numParticulas, anchoTextura, G);

        coordOut = result;
        speedOut = length((coordsActuales-coordsPrevias) / pasoT);
    }
)";


ParticulasApp::ParticulasApp(unsigned numParticulas) :
    GLFWWindow(512, 512, "Partículas"),
    numParticulas(numParticulas) {
    glfwSwapInterval(0);
    makeContextCurrent();
}


void ParticulasApp::mostrarFR() {
    double tiempoActual = glfwGetTime();
    double delta = tiempoActual - tVisUltimoFR;
    numVisFrames++;

    if (delta >= 1.0) {
        double fps = double(numVisFrames) / delta;
        std::stringstream ss;
        ss << "Partículas" << " [" << fps << " FPS]";

        setTitle(ss.str().c_str());

        numVisFrames = 0;
        tVisUltimoFR = tiempoActual;
    }
}


void ParticulasApp::prepararVisualizacion() {
    progVis = crearProgramaShader(codVisVS, codVisFS);

    // Activación de viewport
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Crear y activar VAO de partículas
    glGenVertexArrays(1, &vaoVis);
    glBindVertexArray(vaoVis);

    // Crear y activar VBO para coordenadas de las partículas, asociado al VAO de las partículas
    glGenBuffers(1, &vboVisCoords);
    glBindBuffer(GL_ARRAY_BUFFER, vboVisCoords);

    // Asociar VBO de coordenadas al parámetro vertex del shader
    glEnableVertexAttribArray(attribVisCoord);
    glVertexAttribPointer(attribVisCoord, 2, GL_FLOAT, GL_FALSE, 0, NULL);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glGenBuffers(1, &vboColorParticulas);
    glBindBuffer(GL_ARRAY_BUFFER, vboColorParticulas);
    glBufferData(GL_ARRAY_BUFFER, colorParticulas.size() * sizeof(vec3), &colorParticulas[0], GL_STREAM_DRAW);
    glEnableVertexAttribArray(colorParam);
    glVertexAttribPointer(colorParam, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    tVisUltimoFR = glfwGetTime();
    numVisFrames = 0;
}


void ParticulasApp::visualizar() {
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);

    // Activar shaders
    glUseProgram(progVis);

    // Activar array de vértices
    glBindVertexArray(vaoVis);
    glBindBuffer(GL_ARRAY_BUFFER, vboVisCoords);
    // Transferir coordenadas actualizadas al vbo
    glBufferData(GL_ARRAY_BUFFER, numParticulas * sizeof(Coord), &coords[0], GL_DYNAMIC_DRAW);

    glDrawArrays(GL_POINTS, 0, (GLsizei)numParticulas);

    glBindVertexArray(0);
    glUseProgram(0);

    glDisable(GL_BLEND);
    glDisable(GL_PROGRAM_POINT_SIZE);
    swapBuffers();
}


void ParticulasApp::ejecutar() {
#ifdef VISUALIZAR
    // Bucle de simulación/renderizado
    while (!shouldClose()) {
        simular();
        visualizar();
        mostrarFR();

        glfwPollEvents();
    }
#else
    // Bucle de solo simulación
    while (t < SIM_DURACION) {
        simular();
    }
#endif
    std::cout << "La partícula " << maxParticle << " posee la velocidad máxima (" << maxSpeed << ")\n";
}

void ParticulasApp::generarPuntos() {
    coords.clear();

    for (unsigned c = 0; c < numParticulas; ++c) {
        coords.push_back(Coord(aleatorio(-1.0f, 1.0f), aleatorio(-1.0f, 1.0f)));
        colorParticulas.push_back(vec3(aleatorio(0.0f, 1.0f), aleatorio(0.0f, 1.0f), aleatorio(0.0, 1.0f)));
    }

    /* Introducir coordenadas nulas adicionales para completar
     * un número entero de líneas de la textura
     */
    unsigned restoLineaTextura = numParticulas % SIM_ANCHO_TEXTURA_DATOS;
    if (restoLineaTextura > 0) {
        coords.insert(coords.end(), SIM_ANCHO_TEXTURA_DATOS - restoLineaTextura, Coord(0.0f, 0.0f));
    }
}

vector<ParticulasApp::Coord> ParticulasApp::generarPerturbacion() {
    vector<ParticulasApp::Coord> perturbacion = coords;

    for (unsigned c = 0; c < numParticulas; ++c) {
        perturbacion[c] = (coords[c] + Coord(aleatorio(-0.001f, 0.001f), aleatorio(-0.001f, 0.001f)));
        speeds.push_back(length((coords[c] - perturbacion[c]) / SIM_PASO_TIEMPO));
    }

    return perturbacion;
}

void ParticulasApp::prepararSimulacion() {
    generarPuntos();

    // Crear programa de shading
    progSim = crearProgramaShader(codSimVS, codSimFS);

    glUseProgram(progSim);
    // Obtener uniforme a sample de textura de entrada
    uniformSimCoordsIn = glGetUniformLocation(progSim, "coordsIn");
    uniformSimCoordsPrevious = glGetUniformLocation(progSim, "coordsPrevious");

    // Establecer uniformes constantes con número de partículas,
    // ancho textura y paso de tiempo
    glUniform1i(glGetUniformLocation(progSim, "numParticulas"), numParticulas);
    glUniform1i(glGetUniformLocation(progSim, "anchoTextura"), SIM_ANCHO_TEXTURA_DATOS);
    glUniform1f(glGetUniformLocation(progSim, "pasoT"), SIM_PASO_TIEMPO);
    glUniform1f(glGetUniformLocation(progSim, "G"), CONSTANTE_G);

    glUseProgram(0);

    // Crear texturas con coordenadas de entrada y salida
    texSimCoords[0] = crearTextura2DVec2(numParticulas, (float*)&coords[0]);
    texSimCoords[1] = crearTextura2DVec2(numParticulas, (float*)&generarPerturbacion()[0]); // POSICIONES EN INSTANTE ANTERIOR (MISMAS COORDS = SIN MOVIMIENTO INICIAL)
    texSimCoords[2] = crearTextura2DVec2(numParticulas, 0);
    texSimSpeeds[0] = crearTextura2DVec1(numParticulas, (float*)&speeds[0]);

    // Crear FBO para escribir en las coordenadas de salida
    glGenFramebuffers(1, &fboSimCoords);
    glGenFramebuffers(1, &fboSimSpeeds);

    //    glViewport(0, 0, tamTextura, tamTextura);

        // Preparar geometría para cubrir todo el viewport (dos triángulos)
    float simCoords[] = {
        -1.0f, -1.0f,
        -1.0f, 1.0f,
        1.0f, -1.0f,
        1.0f, 1.0f
    };

    // Generar vao y vbo para geometría
    glGenVertexArrays(1, &vaoSim);
    glBindVertexArray(vaoSim);

    glGenBuffers(1, &vboSimCoords);
    glBindBuffer(GL_ARRAY_BUFFER, vboSimCoords);
    glBufferData(GL_ARRAY_BUFFER, sizeof(simCoords), simCoords, GL_STATIC_DRAW);

    glEnableVertexAttribArray(attribSimCoord);
    glVertexAttribPointer(attribSimCoord, 2, GL_FLOAT, GL_FALSE, 0, 0);

    // Iniciar contador de tiempo
    t = 0;
    maxSpeed = 0;
    maxParticle = 0;
}


void ParticulasApp::simular() {
    // Activar programa de simulación
    glUseProgram(progSim);

    glBindVertexArray(vaoSim);

    // Cargar textura de entrada en unidad 0 y pasar a uniforme
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texSimCoords[0]);
    glUniform1i(uniformSimCoordsIn, 0);

    // Cargar textura de entrada (del instante anterior) en unidad 1 y pasar a uniforme
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, texSimCoords[1]);
    glUniform1i(uniformSimCoordsPrevious, 1);

    // Activar fbo de escritura a textura
    glBindFramebuffer(GL_FRAMEBUFFER, fboSimCoords);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texSimCoords[2], 0);

    // Lanzar computación
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    glActiveTexture(GL_TEXTURE2);
    glBindFramebuffer(GL_FRAMEBUFFER, fboSimSpeeds);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, texSimSpeeds[0], 0);

    glUseProgram(0);

    // Extraer coordenadas de partículas y pasarlas a vector
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, texSimCoords[2]);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RG, GL_FLOAT, &coords[0]);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, texSimSpeeds[0]);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, &speeds[0]);

    // Intercambiar textura de lectura y escritura
    swap(texSimCoords[1], texSimCoords[0]); // Pasamos nuestro instante T al instante anterior
    swap(texSimCoords[0], texSimCoords[2]); // Cargamos las nuevas posiciones al instante T

    // Avanzar variable de tiempo
    t += SIM_PASO_TIEMPO;
    for (int i = 0; i < NUM_PARTICULAS; i++) {
        if (speeds[i] > maxSpeed) {
            maxSpeed = speeds[i];
            maxParticle = i;
            std::cout << "Nueva velocidad máxima: Partícula " << i << "con velocidad " << maxSpeed <<"\n";
        }
    }
}


void ParticulasApp::liberarRecursos() {
#ifdef VISUALIZAR
    glDeleteBuffers(1, &vboVisCoords);
    glDeleteVertexArrays(1, &vaoVis);
    glDeleteBuffers(1, &vboVisCoords);
    glDeleteBuffers(1, &vboColorParticulas);
    glDeleteProgram(progVis);
#endif

    glDeleteBuffers(1, &vboSimCoords);
    glDeleteVertexArrays(1, &vaoSim);
    glDeleteBuffers(1, &vboSimCoords);
    glDeleteTextures(2, texSimCoords);
    glDeleteFramebuffers(1, &fboSimCoords);
    glDeleteProgram(progSim);
}


int main(int argc, char** argv) {
    // Iniciación glfw
    if (glfwInit() != GL_TRUE) {
        std::cerr << "Fallo en inicialización de GLFW" << std::endl;
        return 1;
    }

    // Application & window init
    ParticulasApp::hint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    ParticulasApp::hint(GLFW_CONTEXT_VERSION_MINOR, 1);
    ParticulasApp::hint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    ParticulasApp::hint(GLFW_RESIZABLE, GLFW_FALSE);

#ifndef VISUALIZAR
    ParticulasApp::hint(GLFW_VISIBLE, GLFW_FALSE);
#endif

#ifdef __APPLE__
    ParticulasApp::hint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    ParticulasApp particulasApp(NUM_PARTICULAS);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Fallo al cargar OpenGL" << std::endl;
        return 1;
    }

    particulasApp.prepararSimulacion();
#ifdef VISUALIZAR
    particulasApp.prepararVisualizacion();
#endif
    particulasApp.ejecutar();

    // Terminar GLFW
    glfwTerminate();

    return 0;
}