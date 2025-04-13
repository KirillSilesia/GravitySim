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
#include <string>
#include <limits>
#include <glm/glm.hpp>
#pragma managed(push, off)
#include <occi.h>
#include <oratypes.h>
#pragma managed(pop)
#define OCCI_DEMO_NEED_IMPORT
#define _CRT_SECURE_NO_DEPRECATE

using namespace std;
using namespace oracle::occi;

float zoomLevel = 45.0f;         // Initial zoom level
const float G = 0.66743f;        // Adjusted gravitational constant
double lastX = 320.0, lastY = 240.0;          // Initial mouse position (center of window)
double lastRightX = 320.0, lastRightY = 240.0;  // Right mouse position for dragging
float pitch = 0.0f, yaw = -90.0f; // Camera rotation angles
bool isDragging = false;         // Flag for left mouse button hold
bool isRightDragging = false;    // Flag for right mouse button hold
// Place camera far enough from center-of-mass – we will recenter on the system’s COM below.
float cameraX = 0.0f, cameraY = 0.0f, cameraZ = -30.0f;
int starCounter;

struct Vector3 {
    float x, y, z;
};

struct DatabaseParams {
    std::string user;
    std::string password;
    std::string connectString;
};

DatabaseParams getDatabaseParams() {
    DatabaseParams params;
    cout << "Enter database username: ";
    std::getline(std::cin, params.user);
    cout << "Enter database password: ";
    std::getline(std::cin, params.password);
    cout << "Enter database connect string (e.g., //111.111.111.11:1111/SID): ";
    std::getline(std::cin, params.connectString);
    return params;
}

// Structure to represent celestial bodies
struct CelestialBody {
    string BodyName;
    string BodyType;
    float x, y, z;     // Position
    float size;        // Size
    float mass;        // Mass
    float color[3];    // Color
    float vx, vy, vz;  // Velocity
    bool isStar;
};

// Function to draw a celestial body (sphere)
void drawSphere(float radius, float x, float y, float z, const float* color) {
    glPushMatrix();
    glTranslatef(x, y, z);  // Move the sphere to its position
    glColor3f(color[0], color[1], color[2]);
    GLUquadric* quadric = gluNewQuadric();
    gluSphere(quadric, radius, 20, 20);
    gluDeleteQuadric(quadric);
    glPopMatrix();
}

// Function to draw the gravity grid (for visualization purposes)
void drawGravityGrid(int gridSize, float lineSpacing, vector<CelestialBody>& bodies) {
    int numLines = gridSize / lineSpacing;
    float gravityScale = 50.f;

    glColor3f(0.9f, 0.9f, 0.9f);
    glBegin(GL_LINES);

    vector<vector<float>> heights(numLines + 1, vector<float>(numLines + 1, 0.0f));

    // Calculate the height values based on gravitational force at each point
    for (int i = 0; i <= numLines; ++i) {
        for (int j = 0; j <= numLines; ++j) {
            float x = (i - numLines / 2) * lineSpacing;
            float z = (j - numLines / 2) * lineSpacing;
            float y = 0.0f;

            for (size_t b = 0; b < bodies.size(); ++b) {
                CelestialBody& body = bodies[b];

                float dx = x - body.x;
                float dz = z - body.z;
                float distanceSquared = dx * dx + dz * dz + 1.0f;  // prevent division by zero
                float distance = sqrt(distanceSquared);

                float force = (G * body.mass) / distanceSquared;
                y -= force * gravityScale / (1.0f + distance * 0.1f);
            }
            heights[i][j] = y;
        }
    }

    // Draw grid lines along X
    for (int i = 0; i < numLines; ++i) {
        for (int j = 0; j <= numLines; ++j) {
            float x1 = (i - numLines / 2) * lineSpacing;
            float x2 = ((i + 1) - numLines / 2) * lineSpacing;
            float z = (j - numLines / 2) * lineSpacing;
            glVertex3f(x1, heights[i][j], z);
            glVertex3f(x2, heights[i + 1][j], z);
        }
    }

    // Draw grid lines along Z
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
    if (zoomLevel < 10.0f)
        zoomLevel = 10.0f;
    if (zoomLevel > 90.0f)
        zoomLevel = 90.0f;
}

// Mouse movement callback
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse)
        return;

    int leftButtonState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    int rightButtonState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);

    if (leftButtonState == GLFW_PRESS) {
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
        float xOffset = xpos - lastRightX;
        float yOffset = lastRightY - ypos;
        lastRightX = xpos;
        lastRightY = ypos;
        float sensitivity = 0.05f;
        cameraX -= xOffset * sensitivity;
        cameraY += yOffset * sensitivity;
    }
}

// Mouse button callback
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse)
        return;
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
        glfwGetCursorPos(window, &lastX, &lastY);
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
        glfwGetCursorPos(window, &lastRightX, &lastRightY);
}

// Calculate center-of-mass of all bodies
glm::vec3 calculateCenterOfMass(vector<CelestialBody>& bodies) {
    glm::vec3 centerOfMass(0.0f);
    float totalMass = 0.0f;
    for (size_t i = 0; i < bodies.size(); ++i) {
        totalMass += bodies[i].mass;
        centerOfMass += bodies[i].mass * glm::vec3(bodies[i].x, bodies[i].y, bodies[i].z);
    }
    if (totalMass > 0.0f)
        centerOfMass /= totalMass;
    return centerOfMass;
}

// Callback for window resizing
void reshape(GLFWwindow* window, int width, int height) {
    if (height == 0)
        height = 1;
    float aspectRatio = (float)width / height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(zoomLevel, aspectRatio, 0.1f, 1000.0f);  // increased far plane
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// Update the positions and velocities of celestial bodies based on gravity
void updatePositionsAndVelocities(vector<CelestialBody>& bodies, float dt) {
    vector<glm::vec3> accelerations(bodies.size(), glm::vec3(0.0f));
    // Compute gravitational forces between all bodies
    for (size_t i = 0; i < bodies.size(); ++i) {
        CelestialBody& body1 = bodies[i];
        for (size_t j = 0; j < bodies.size(); ++j) {
            if (i != j) {
                CelestialBody& body2 = bodies[j];
                glm::vec3 r_vec(body2.x - body1.x, body2.y - body1.y, body2.z - body1.z);
                float r2 = glm::dot(r_vec, r_vec);
                if (r2 < 1e-6f) continue;
                float r_dist = sqrt(r2);
                float force = (G * body1.mass * body2.mass) / r2;
                glm::vec3 acceleration = (force / body1.mass) * glm::normalize(r_vec);
                accelerations[i] += acceleration;
            }
        }
    }
    // Update velocities and positions
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].vx += accelerations[i].x * dt;
        bodies[i].vy += accelerations[i].y * dt;
        bodies[i].vz += accelerations[i].z * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }
}

// Initialize orbital velocities based on center-of-mass and system-dependent values
void initializeOrbitalVelocities(vector<CelestialBody>& bodies) {
    // Separate stars and non-stars
    vector<CelestialBody*> stars;
    vector<CelestialBody*> planets;
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].isStar)
            stars.push_back(&bodies[i]);
        else
            planets.push_back(&bodies[i]);
    }

    // Move these outside so they're visible to both sections
    glm::vec3 starsCOM(0.0f);
    float totalStarMass = 0.0f;

    // Compute stars COM and total mass
    for (auto star : stars) {
        starsCOM += star->mass * glm::vec3(star->x, star->y, star->z);
        totalStarMass += star->mass;
    }
    if (totalStarMass > 0.0f)
        starsCOM /= totalStarMass;

    // Planets orbiting the COM of the stars
    for (auto planet : planets) {
        float dx = planet->x - starsCOM.x;
        float dz = planet->z - starsCOM.z;
        float r = sqrt(dx * dx + dz * dz);
        if (r < 1e-6f) continue;
        float orbitalSpeed = sqrt(G * totalStarMass / r);
        planet->vx = -dz / r * orbitalSpeed;
        planet->vz = dx / r * orbitalSpeed;
    }

    // Stars orbiting the COM of the stars
    if (stars.size() >= 2) {
        for (auto star : stars) {
            float dx = star->x - starsCOM.x;
            float dz = star->z - starsCOM.z;
            float r = sqrt(dx * dx + dz * dz);
            if (r < 1e-6f) continue;
            float effectiveMass = totalStarMass - star->mass;
            float orbitalSpeed = sqrt(G * effectiveMass / r);
            star->vx = -dz / r * orbitalSpeed;
            star->vz = dx / r * orbitalSpeed;
        }
    } else if (stars.size() == 1) {
        stars[0]->vx = stars[0]->vy = stars[0]->vz = 0.0f;
    }
}


// Function to initialize the OpenGL window
GLFWwindow* initializeOpenGL() {
    if (!glfwInit())
        return nullptr;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(800, 600, "Gravity Simulator", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return nullptr;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    return window;
}

// Database fetch functions remain unchanged but add logging for debugging.
vector<string> fetchSolarSystemNames(Connection* conn) {
    vector<string> solarSystemNames;
    Statement* stmt = nullptr;
    ResultSet* rs = nullptr;
    try {
        stmt = conn->createStatement();
        stmt->setSQL("SELECT SolarSystemName FROM SolarSystems");
        rs = stmt->executeQuery();
        while (rs->next()) {
            solarSystemNames.push_back(rs->getString(1));
        }
        conn->terminateStatement(stmt);
    }
    catch (SQLException& ex) {
        cout << "Error fetching solar system names: " << ex.getMessage() << endl;
    }
    return solarSystemNames;
}

vector<CelestialBody> fetchCelestialBodies(Connection* conn, int solarSystemID) {
    vector<CelestialBody> celestialBodies;
    Statement* stmt = nullptr;
    ResultSet* rs = nullptr;
    try {
        stmt = conn->createStatement();
        stmt->setSQL("SELECT BodyName, BodyType, Mass, BodySize, Color, X, Y, Z, VX, VY, VZ, IsStar FROM CelestialBodies WHERE SolarSystemID = :1");
        stmt->setInt(1, solarSystemID);
        rs = stmt->executeQuery();
        while (rs->next()) {
            CelestialBody body;
            body.BodyName = rs->getString(1);
            body.BodyType = rs->getString(2);
            body.mass = rs->getDouble(3);
            body.size = rs->getFloat(4);
            string colorStr = rs->getString(5);
            float r, g, b;
            if (sscanf_s(colorStr.c_str(), "%f,%f,%f", &r, &g, &b) == 3) {
                body.color[0] = r;
                body.color[1] = g;
                body.color[2] = b;
            }
            else {
                body.color[0] = 1.0f;
                body.color[1] = 1.0f;
                body.color[2] = 1.0f;
                cout << "Failed to parse color string: " << colorStr << endl;
            }
            body.x = rs->getFloat(6);
            body.y = rs->getFloat(7);
            body.z = rs->getFloat(8);
            body.vx = rs->getFloat(9);
            body.vy = rs->getFloat(10);
            body.vz = rs->getFloat(11);
            body.isStar = rs->getInt(12) == 1;
            celestialBodies.push_back(body);
        }
        conn->terminateStatement(stmt);
    }
    catch (SQLException& ex) {
        cout << "Error fetching celestial bodies: " << ex.getMessage() << endl;
    }
    return celestialBodies;
}

// Update simulation bodies from database
void updateBodies(Connection* conn, int solarSystemID, vector<CelestialBody>& bodies) {
    vector<CelestialBody> fetchedBodies = fetchCelestialBodies(conn, solarSystemID);
    if (fetchedBodies.empty()) {
        cout << "No celestial bodies found for SolarSystemID " << solarSystemID << endl;
    }
    bodies = fetchedBodies;
}

int main(void) {
    DatabaseParams dbParams = getDatabaseParams();
    Environment* env = nullptr;
    Connection* conn = nullptr;
    Statement* stmt = nullptr;
    ResultSet* rs = nullptr;
    GLFWwindow* window = nullptr;
    vector<CelestialBody> bodiesVector;
    int numBodies = 0;
    CelestialBody* bodies = nullptr;

    try {
        string user = dbParams.user;
        string password = dbParams.password;
        string connectString = dbParams.connectString;
        env = Environment::createEnvironment(Environment::DEFAULT);
        conn = env->createConnection(user, password, connectString);
        cout << "Connected to the database!" << endl;

        stmt = conn->createStatement();
        stmt->setSQL("SELECT 1 FROM DUAL");
        rs = stmt->executeQuery();
        if (rs->next()) {
            cout << "Database test query successful: " << rs->getInt(1) << endl;
        }
        else {
            cout << "Database test query failed." << endl;
        }
        conn->terminateStatement(stmt);
        stmt = nullptr;
    }
    catch (SQLException& ex) {
        cout << "Error: " << ex.getMessage() << endl;
        return -1;
    }

    if (!glfwInit())
        return -1;

    window = glfwCreateWindow(640, 480, "Gravity Simulator", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glfwSetFramebufferSizeCallback(window, reshape);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glEnable(GL_DEPTH_TEST);
    reshape(window, 640, 480);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(zoomLevel, 640.0 / 480.0, 0.1, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    float lastTime = glfwGetTime();
    float dt = 0.01f;   // delta time for simulation (in seconds)

    static int selectedSolarSystemIndex = 0;
    vector<string> solarSystemNames = fetchSolarSystemNames(conn);

    // Load initial solar system (if available)
    if (!solarSystemNames.empty()) {
        updateBodies(conn, selectedSolarSystemIndex + 1, bodiesVector);
        numBodies = bodiesVector.size();
        cout << "Loaded " << numBodies << " bodies" << endl;
        initializeOrbitalVelocities(bodiesVector);
    }

    static int previousSolarSystemIndex = -1;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        float currentTime = glfwGetTime();
        dt = currentTime - lastTime;
        lastTime = currentTime;

        // Update ImGui
        ImGui_ImplGlfw_NewFrame();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui::NewFrame();
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x, io.DisplaySize.y * 0.2f));
        ImGui::SetNextWindowPos(ImVec2(0.0f, io.DisplaySize.y * 0.8f), ImGuiCond_Always);

        ImGui::Begin("Simulation Controls",
            nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);
        ImGui::SliderFloat("Zoom", &zoomLevel, 10.0f, 90.0f);

        if (!solarSystemNames.empty()) {
            const char* combo_preview_value = solarSystemNames[selectedSolarSystemIndex].c_str();
            if (ImGui::BeginCombo("Solar System", combo_preview_value)) {
                for (int n = 0; n < solarSystemNames.size(); n++) {
                    bool is_selected = (selectedSolarSystemIndex == n);
                    if (ImGui::Selectable(solarSystemNames[n].c_str(), is_selected)) {
                        selectedSolarSystemIndex = n;
                    }
                    if (is_selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        }
        ImGui::End();

        // Only update the bodies if a new system is selected.
        if (selectedSolarSystemIndex != previousSolarSystemIndex) {
            // SolarSystemID is assumed to be (index + 1)
            updateBodies(conn, selectedSolarSystemIndex + 1, bodiesVector);
            initializeOrbitalVelocities(bodiesVector);
            previousSolarSystemIndex = selectedSolarSystemIndex;
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(zoomLevel, 640.0f / 480.0f, 0.1f, 1000.0f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Recenter the scene based on the system’s center-of-mass.
        glm::vec3 systemCoM = calculateCenterOfMass(bodiesVector);
        glTranslatef(-systemCoM.x + cameraX, -systemCoM.y + cameraY, cameraZ);
        glRotatef(pitch, 1.0f, 0.0f, 0.0f);
        glRotatef(yaw, 0.0f, 1.0f, 0.0f);

        drawGravityGrid(100, 0.25f, bodiesVector);
        updatePositionsAndVelocities(bodiesVector, dt);

        // Draw each celestial body.
        for (const auto& body : bodiesVector) {
            drawSphere(body.size, body.x, body.y, body.z, body.color);
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    if (conn) {
        try {
            env->terminateConnection(conn);
            Environment::terminateEnvironment(env);
            cout << "Disconnected from the database." << endl;
        }
        catch (SQLException& ex) {
            cout << "Error during disconnection: " << ex.getMessage() << endl;
        }
    }
    return 0;
}
