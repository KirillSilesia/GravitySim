#include <windows.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>

using namespace std;

float zoomLevel = 45.0f; // Initial zoom level
const float G = 6.67430e-11f; // Gravitational constant
float deltaTime = 0.1f; // Time step
double lastX = 320.0, lastY = 240.0; // Initial mouse position (center of window)
double lastRightX = 320.0, lastRightY = 240.0; // Right mouse position for dragging
float pitch = 0.0f, yaw = -90.0f;     // Camera rotation angles
bool isDragging = false;               // Flag for left mouse button hold
bool isRightDragging = false; // Flag for right mouse button hold
float cameraX = 0.0f, cameraY = 0.0f, cameraZ = -15.0f; // Camera position
int starCounter;

struct Vector3 {
    float x, y, z;
};

// Structure to represent celestial bodies
struct CelestialBody {
    float x, y, z;  // Position
    float size;     // Size
    float mass;     // Mass
    float color[3]; // Color
    float vx, vy, vz; // Velocity
    float initialVx, initialVz;
    bool isStar;
};

// Function to draw a celestial body (sphere)
void drawSphere(float radius, float x, float y, float z, float* color) {
    glPushMatrix();
    glTranslatef(x, y, z);  // Move the sphere to the position of the body
    glColor3f(color[0], color[1], color[2]);  // Set the color of the sphere
    GLUquadric* quadric = gluNewQuadric();  // Create a new quadratic object for the sphere
    gluSphere(quadric, radius, 20, 20);  // Draw the sphere with the specified radius
    glPopMatrix();  // Restore the matrix
}

// Function to draw the gravity grid (for visualization purposes)
void drawGravityGrid(int gridSize, float lineSpacing, CelestialBody* bodies, int numBodies) {
    int numLines = gridSize / lineSpacing;
    float gravityScale = 5.0f;

    glColor3f(0.5f, 0.5f, 0.5f);  // Set the color for the grid
    glBegin(GL_LINES);

    std::vector<std::vector<float>> heights(numLines + 1, std::vector<float>(numLines + 1, 0.0f));

    // Calculate the height values based on the gravitational force at each point
    for (int i = 0; i <= numLines; ++i) {
        for (int j = 0; j <= numLines; ++j) {
            float x = (i - numLines / 2) * lineSpacing;
            float z = (j - numLines / 2) * lineSpacing;
            float y = 0.0f;

            for (int b = 0; b < numBodies; ++b) {
                CelestialBody& body = bodies[b];

                float dx = x - body.x;
                float dz = z - body.z;
                float distanceSquared = dx * dx + dz * dz + 1.0f;  // Add small value to avoid division by zero
                float distance = sqrt(distanceSquared);

                // Calculate gravitational force and apply it to the height at the point
                float force = (G * body.mass) / distanceSquared;
                y -= force * gravityScale / (1.0f + distance * 0.1f);

            }

            heights[i][j] = y;
        }
    }

    // Draw the grid lines in the X direction
    for (int i = 0; i < numLines; ++i) {
        for (int j = 0; j <= numLines; ++j) {
            float x1 = (i - numLines / 2) * lineSpacing;
            float x2 = ((i + 1) - numLines / 2) * lineSpacing;
            float z = (j - numLines / 2) * lineSpacing;

            glVertex3f(x1, heights[i][j], z);
            glVertex3f(x2, heights[i + 1][j], z);
        }
    }

    // Draw the grid lines in the Z direction
    for (int i = 0; i <= numLines; ++i) {
        for (int j = 0; j < numLines; ++j) {
            float x = (i - numLines / 2) * lineSpacing;
            float z1 = (j - numLines / 2) * lineSpacing;
            float z2 = ((j + 1) - numLines / 2) * lineSpacing;

            glVertex3f(x, heights[i][j], z1);
            glVertex3f(x, heights[i][j + 1], z2);
        }
    }

    glEnd();
}

// Callback for zoom
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) return;

    float zoomSpeed = 5.0f;
    zoomLevel -= yoffset * zoomSpeed;
}


// Mouse callback function
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {

    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) {
        return; // Check if the mouse is in GUI
    }
    
        int leftButtonState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        int rightButtonState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);

        if (leftButtonState == GLFW_PRESS) {
            // Rotate the camera
            float xOffset = xpos - lastX;
            float yOffset = lastY - ypos;
            lastX = xpos;
            lastY = ypos;

            float sensitivity = 0.1f;
            xOffset *= sensitivity;
            yOffset *= sensitivity;

            yaw += xOffset;
            pitch += yOffset;

            if (pitch > 89.0f) pitch = 89.0f;
            if (pitch < -89.0f) pitch = -89.0f;
        }

        if (rightButtonState == GLFW_PRESS) {
            // Move the camera
            float xOffset = xpos - lastRightX;
            float yOffset = lastRightY - ypos;
            lastRightX = xpos;
            lastRightY = ypos;

            float sensitivity = 0.05f;
            cameraX -= xOffset * sensitivity;
            cameraY += yOffset * sensitivity;
        }
}



void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) return;

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        glfwGetCursorPos(window, &lastX, &lastY);
    }

    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        glfwGetCursorPos(window, &lastRightX, &lastRightY);
    }
}

// Calculating stars center of mass
glm::vec3 calculateCenterOfMass(CelestialBody* bodies, int numBodies) {
    glm::vec3 centerOfMass(0.0f);
    float totalMass = 0.0f;

    for (int i = 0; i < numBodies; ++i) {
        if (bodies[i].isStar) {  // Only stars contribute to the COM
            totalMass += bodies[i].mass;
            centerOfMass += bodies[i].mass * glm::vec3(bodies[i].x, bodies[i].y, bodies[i].z);
        }
    }

    if (totalMass > 0.0f) {
        centerOfMass /= totalMass;
    }

    return centerOfMass;
}


// Function to handle window resizing
void reshape(GLFWwindow* window, int width, int height) {
    if (height == 0) height = 1;
    float aspectRatio = (float)width / (float)height;

    glViewport(0, 0, width, height);  // Set the viewport size
    glMatrixMode(GL_PROJECTION);  // Set projection matrix mode
    glLoadIdentity();
    gluPerspective(zoomLevel, aspectRatio, 0.1f, 100.0f);  // Set the perspective projection
    glMatrixMode(GL_MODELVIEW);  // Set modelview matrix mode
    glLoadIdentity();
}

// Function to update star velocities
void updateStarVelocities(CelestialBody* bodies, int numBodies, float deltaTime) {
    for (int i = 0; i < numBodies; ++i) {
        if (!bodies[i].isStar) continue;  // Skip non-star bodies

        glm::vec3 totalForce(0.0f);

        for (int j = 0; j < numBodies; ++j) {
            if (i == j || !bodies[j].isStar) continue;  // Don't calculate force on itself or non-stars

            glm::vec3 direction = glm::vec3(bodies[j].x, bodies[j].y, bodies[j].z) - glm::vec3(bodies[i].x, bodies[i].y, bodies[i].z);
            float distance = glm::length(direction);
            direction = glm::normalize(direction);

            // Gravitational force magnitude
            float forceMagnitude = (G * bodies[i].mass * bodies[j].mass) / (distance * distance);

            // Force vector
            glm::vec3 force = forceMagnitude * direction;
            totalForce += force;
        }

        // Update velocity (velocity += force * time step / mass)
        bodies[i].vx += totalForce.x / bodies[i].mass * deltaTime;
        bodies[i].vy += totalForce.y / bodies[i].mass * deltaTime;
        bodies[i].vz += totalForce.z / bodies[i].mass * deltaTime;
    }
}


// Function to update the positions and velocities of celestial bodies based on gravity
void updatePositionsAndVelocities(CelestialBody* bodies, int numBodies, float deltaTime) {
    // Temporary arrays to store acceleration values
    float* ax = new float[numBodies]();
    float* az = new float[numBodies]();

    // Compute gravitational forces between all bodies
    for (int i = 0; i < numBodies; ++i) {
        for (int j = 0; j < numBodies; ++j) {
            if (i == j) continue;  // Skip self-interaction

            CelestialBody& body1 = bodies[i];
            CelestialBody& body2 = bodies[j];

            float dx = body2.x - body1.x;
            float dz = body2.z - body1.z;
            float r2 = dx * dx + dz * dz;  // Squared distance

            if (r2 < 1e-6f) continue;  // Prevent division by zero or extremely small values

            float r = sqrt(r2);
            float force = G * body1.mass * body2.mass / r2;  // Gravitational force magnitude

            // Acceleration components
            ax[i] += force * dx / (r * body1.mass);
            az[i] += force * dz / (r * body1.mass);
        }
    }

    // Update velocities and positions
    for (int i = 0; i < numBodies; ++i) {
        CelestialBody& body = bodies[i];

        body.vx += ax[i] * deltaTime;
        body.vz += az[i] * deltaTime;

        body.x += body.vx * deltaTime;
        body.z += body.vz * deltaTime;
    }

    // Free dynamically allocated arrays
    delete[] ax;
    delete[] az;
}

void initializeOrbitalVelocities(CelestialBody* bodies, int numBodies) {
    // Calculate the center of mass of the star system
    glm::vec3 CoM = calculateCenterOfMass(bodies, numBodies);

    for (int i = 0; i < numBodies; ++i) {
        CelestialBody& body = bodies[i];

        // Skip center of mass calculations for objects already at the CoM
        float dx = body.x - CoM.x;
        float dz = body.z - CoM.z;
        float r2 = dx * dx + dz * dz;  // Squared distance
        if (r2 < 1e-6f) continue;  // Avoid division by zero if already at the center

        float r = sqrt(r2);  // Distance from the center of mass

        // Determine which bodies contribute to gravity (only stars)
        float starMassSum = 0.0f;
        for (int j = 0; j < numBodies; ++j) {
            if (bodies[j].isStar) {
                starMassSum += bodies[j].mass;
            }
        }

        // Compute orbital velocity for a circular orbit around the CoM
        float velocityMagnitude = sqrt(G * starMassSum / r);

        // Set velocity perpendicular to the radial vector
        body.vx = -dz * velocityMagnitude / r;
        body.vz = dx * velocityMagnitude / r;
    }
}


// Function to initialize the OpenGL window
GLFWwindow* initializeOpenGL() {
    if (!glfwInit()) {
        return nullptr;  // Initialization failed
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(800, 600, "Gravity Simulator", nullptr, nullptr);

    if (!window) {
        glfwTerminate();  // Window creation failed
        return nullptr;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    
    return window;
}

// Test objects (celestial bodies)
CelestialBody star1 = { 0.0f, 0.0f, 0.0f, 5.0f, 100.0e10f, {1.0f, 1.0f, 0.0f}, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, true };
CelestialBody star2 = { 10.0f, 0.0f, 0.0f, 0.5f, 2.0e10f, {1.0f, 0.0f, 0.0f}, -0.2f, 0.0f, 0.1f, 0.0f, 0.0f, true };
CelestialBody planet = { 5.0f, 0.0f, 0.0f, 0.3f, 5.0e9f, {0.0f, 0.0f, 1.0f}, 0.0f, 0.0f, 0.0f, false };


CelestialBody bodies[] = { star1, star2, planet };

Vector3 initialPositions[] = {
    { bodies[0].x, bodies[0].y, bodies[0].z },
    { bodies[1].x, bodies[1].y, bodies[1].z },
    { bodies[2].x, bodies[2].y, bodies[2].z }
};

Vector3 initialVelocities[] = {
    { bodies[0].vx, bodies[0].vy, bodies[0].vz },
    { bodies[1].vx, bodies[1].vy, bodies[1].vz },
    { bodies[2].vx, bodies[2].vy, bodies[2].vz }
};

void resetSimulation() {
    for (int i = 0; i < 3; i++) {
        bodies[i].x = initialPositions[i].x;
        bodies[i].y = initialPositions[i].y;
        bodies[i].z = initialPositions[i].z;

        bodies[i].vx = initialVelocities[i].x;
        bodies[i].vy = initialVelocities[i].y;
        bodies[i].vz = initialVelocities[i].z;
    }
}

int main(void) {
    GLFWwindow* window;

    // Initialize GLFW
    if (!glfwInit())
        return -1;

    // Create a window
    window = glfwCreateWindow(640, 480, "Gravity Simulator", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  // Enable V-Sync
    glfwSetFramebufferSizeCallback(window, reshape);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);  // For mouse click handling
    glfwSetCursorPosCallback(window, mouse_callback);  // For rotation
    glEnable(GL_DEPTH_TEST); // Enable depth testing for proper 3D rendering
    reshape(window, 640, 480);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 640.0 / 480.0, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    initializeOrbitalVelocities(bodies, 3);
    float lastTime = glfwGetTime();
    float deltaTime;


    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
       
        float currentTime = glfwGetTime();
        deltaTime = currentTime - lastTime;
        lastTime = currentTime;

        // Update ImGui
        ImGui_ImplGlfw_NewFrame();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui::NewFrame();
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x, io.DisplaySize.y * 0.2f));
        ImGui::SetNextWindowPos(ImVec2(0.0f, io.DisplaySize.y * 0.8f), ImGuiCond_Always);

        // GUI Rendering
        ImGui::Begin("Simulation Controls", nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);

        // GUI elements
        ImGui::SliderFloat("Zoom", &zoomLevel, 10.0f, 90.0f);
        if (ImGui::Button("Reset Simulation")) {
            resetSimulation();

        }

        ImGui::End();

        // Cler the frame before rendering
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(zoomLevel, 640.0f / 480.0f, 0.1f, 100.0f);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glTranslatef(cameraX, cameraY, cameraZ); 
        glRotatef(pitch, 1.0f, 0.0f, 0.0f); 
        glRotatef(yaw, 0.0f, 1.0f, 0.0f);

        drawGravityGrid(100, 0.25f, bodies, 3);

        updatePositionsAndVelocities(bodies, 3, deltaTime);

        for (int i = 0; i < 3; ++i) {
            drawSphere(bodies[i].size, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].color);
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);  
    }



    // Clean up ImGui
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwTerminate();  // Clean up GLFW resources
    return 0;
}
