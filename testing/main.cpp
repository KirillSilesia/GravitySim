#include <windows.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS
#define OCCI_DEMO_NEED_IMPORT
#include "GL/glew.h"
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <limits>
#include <glm/glm.hpp>
#pragma managed(push, off)
#include <occi.h>
#include <oratypes.h>
#pragma managed(pop)
#define STBI_NO_PKM
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using namespace std;
using namespace oracle::occi;

float zoomLevel = 45.0f;         // Initial zoom level
const float G = 0.66743f;        // Adjusted gravitational constant
double lastX = 320.0, lastY = 240.0;          // Initial mouse position (center of window)
double lastRightX = 320.0, lastRightY = 240.0;  // Right mouse position for dragging
float pitch = 0.0f, yaw = -90.0f; // Camera rotation angles
bool isDragging = false;         // Flag for left mouse button hold
bool isRightDragging = false;    // Flag for right mouse button hold
float cameraX = 0.0f, cameraY = 0.0f, cameraZ = -30.0f;
int starCounter;
enum CameraFocusType { FOCUS_COM = 0, FOCUS_OBJECT = 1 };
int cameraFocusIndex = 0; // 0 = COM, 1..N = celestial body index + 1

struct Vector3 {
    float x, y, z;
};

// Function to get the base directory of a project
std::string getExecutablePath() {
    char buffer[MAX_PATH];
    GetModuleFileNameA(NULL, buffer, MAX_PATH);
    std::string::size_type pos = std::string(buffer).find_last_of("\\/");
    return std::string(buffer).substr(0, pos);
}

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
    float x, y, z;
    float size;
    float mass;
    GLuint textureID;
    float vx, vy, vz;
    bool isStar;
};

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
    }
    else if (stars.size() == 1) {
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

std::string getTextureFileForBody(const std::string& bodyName, const std::string& bodyType) {
    std::string name = bodyName;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    const std::string folder = (bodyType == "Star")
        ? "Textures/Star_Textures/"
        : "Textures/Planet_Textures/";
    std::string candidate = folder + name + ".jpg";

    // Use the executable path for the existence check
    std::string exeFolder = getExecutablePath();
    std::string fullCandidate = exeFolder + "\\" + candidate;

    if (std::ifstream(fullCandidate, std::ios::binary)) {
        return candidate;
    }
    else {
        std::cerr << "Texture not found: " << candidate << " - using default\n";
        return (bodyType == "Star")
            ? "Textures/Star_Textures/default_star.jpg"
            : "Textures/Planet_Textures/default_planet.jpg";
    }
}

GLuint loadTexture(const std::string& fullPath) {
    int w, h, channels;
    unsigned char* data = stbi_load(fullPath.c_str(), &w, &h, &channels, 0);
    if (!data) {
        std::cerr << "Failed to load " << fullPath << "\n";
        return 0;
    }
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    GLenum fmt = (channels == 4 ? GL_RGBA :
        channels == 3 ? GL_RGB :
        GL_RED);
    glTexImage2D(GL_TEXTURE_2D, 0, fmt, w, h, 0, fmt, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);

    stbi_image_free(data);
    glBindTexture(GL_TEXTURE_2D, 0);
    return tex;
}

void drawSphere(float radius, float x, float y, float z, GLuint textureID) {
    glPushMatrix();
    glTranslatef(x, y, z);

    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);

    if (textureID != 0) {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, textureID);
    }
    else {
        glDisable(GL_TEXTURE_2D);
    }

    GLUquadric* quadric = gluNewQuadric();
    gluQuadricTexture(quadric, textureID != 0 ? GL_TRUE : GL_FALSE);
    gluQuadricNormals(quadric, GLU_SMOOTH);
    gluSphere(quadric, radius, 32, 32);
    gluDeleteQuadric(quadric);

    if (textureID != 0) {
        glBindTexture(GL_TEXTURE_2D, 0);
        glDisable(GL_TEXTURE_2D);
    }

    glPopMatrix();
}

// Database fetch functions remain unchanged but add logging for debugging.
vector<pair<int, string>> fetchSolarSystems(Connection* conn) {
    vector<pair<int, string>> systems;
    Statement* stmt = conn->createStatement();
    stmt->setSQL("SELECT SolarSystemID, SolarSystemName FROM SolarSystems ORDER BY SolarSystemID");
    ResultSet* rs = stmt->executeQuery();
    while (rs->next()) {
        int   id = rs->getInt(1);
        auto  name = rs->getString(2);
        systems.emplace_back(id, name);
    }
    conn->terminateStatement(stmt);
    return systems;
}

vector<CelestialBody> fetchCelestialBodies(Connection* conn, int solarSystemID) {
    vector<CelestialBody> celestialBodies;
    Statement* stmt = nullptr;
    ResultSet* rs = nullptr;
    // debug: how many rows *should* there be?
    {
        Statement* cstmt = conn->createStatement();
        cstmt->setSQL(
            "SELECT COUNT(*) FROM CelestialBodies WHERE SolarSystemID = :1"
        );
        cstmt->setInt(1, solarSystemID);
        ResultSet* crs = cstmt->executeQuery();
        if (crs->next()) {
            std::cout
                << "DEBUG: DB says "
                << crs->getInt(1)
                << " bodies for SolarSystemID="
                << solarSystemID
                << "\n";
        }
        conn->terminateStatement(cstmt);
    }
    try {
        stmt = conn->createStatement();
        stmt->setSQL("SELECT BodyName, BodyType, Mass, BodySize, X, Y, Z, VX, VY, VZ, IsStar FROM CelestialBodies WHERE SolarSystemID = :1");
        stmt->setInt(1, solarSystemID);
        rs = stmt->executeQuery();

        while (rs->next()) {
            CelestialBody body;
            body.BodyName = rs->getString(1);
            body.BodyType = rs->getString(2);
            body.mass = rs->getDouble(3);
            body.size = rs->getFloat(4);
            body.x = rs->getFloat(5);
            body.y = rs->getFloat(6);
            body.z = rs->getFloat(7);
            body.vx = rs->getFloat(8);
            body.vy = rs->getFloat(9);
            body.vz = rs->getFloat(10);
            body.isStar = rs->getInt(11) == 1;

            std::string name = body.BodyName;
            std::string type = body.BodyType;
            std::string relPath = getTextureFileForBody(name, type);
            std::string exeFolder = getExecutablePath();

            // Load the texture (assuming this code exists but was cut off)
            body.textureID = loadTexture(exeFolder + "\\" + relPath);

            // Add this line to store the body in the vector
            celestialBodies.push_back(body);
        }
        if (rs) {
            stmt->closeResultSet(rs);
            rs = nullptr;
        }
        if (stmt) {
            conn->terminateStatement(stmt);
            stmt = nullptr;
        }
    }
    catch (SQLException& e) {
        std::cerr << "Database error: " << e.getMessage() << std::endl;
        // Clean up if error occurred
        if (rs) {
            try { stmt->closeResultSet(rs); }
            catch (...) {}
            rs = nullptr;
        }
        if (stmt) {
            try { conn->terminateStatement(stmt); }
            catch (...) {}
            stmt = nullptr;
        }
    }
    catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        // Clean up if error occurred
        if (rs) {
            try { stmt->closeResultSet(rs); }
            catch (...) {}
            rs = nullptr;
        }
        if (stmt) {
            try { conn->terminateStatement(stmt); }
            catch (...) {}
            stmt = nullptr;
        }
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
    GLFWwindow* window = nullptr;
    int cameraFocusIndex = 0;

    try {
        string user = dbParams.user;
        string password = dbParams.password;
        string connectString = dbParams.connectString;
        env = Environment::createEnvironment(Environment::DEFAULT);
        conn = env->createConnection(user, password, connectString);
        cout << "Connected to the database!" << endl;

        // Test query
        Statement* stmt = conn->createStatement();
        stmt->setSQL("SELECT 1 FROM DUAL");
        ResultSet* rs = stmt->executeQuery();
        if (rs->next()) {
            cout << "Database test query successful: " << rs->getInt(1) << endl;
        }
        conn->terminateStatement(stmt);
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
    glewInit();
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
    float dt = 0.01f;

    // --- Solar System Selection Logic ---
    auto systems = fetchSolarSystems(conn);
    int selectedIndex = 0; // Current selection in the combo
    int previousIndex = -1; // Previous selection to detect changes

    // Load the initial system
    vector<CelestialBody> bodiesVector;
    if (!systems.empty()) {
        bodiesVector = fetchCelestialBodies(conn, systems[selectedIndex].first);
        initializeOrbitalVelocities(bodiesVector);
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        float currentTime = glfwGetTime();
        dt = currentTime - lastTime;
        lastTime = currentTime;

        // --- IMGUI UI ---
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

        if (!systems.empty()) {
            const char* comboPreview = systems[selectedIndex].second.c_str();
            if (ImGui::BeginCombo("Solar System", comboPreview)) {
                for (int n = 0; n < (int)systems.size(); n++) {
                    const bool isSelected = (selectedIndex == n);
                    if (ImGui::Selectable(systems[n].second.c_str(), isSelected)) {
                        selectedIndex = n; // Update selection immediately
                    }
                    if (isSelected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        }

        if (!bodiesVector.empty()) {
            static std::vector<std::string> focusOptions;
            focusOptions.clear();
            focusOptions.push_back("Center of Mass (COM)");
            for (const auto& body : bodiesVector)
                focusOptions.push_back(body.BodyName);

            // Directly use cameraFocusIndex for selection so it stays synced
            if (ImGui::Combo("Camera Center", &cameraFocusIndex,
                [](void* data, int idx, const char** out_text) {
                    const auto& opts = *static_cast<const std::vector<std::string>*>(data);
                    *out_text = opts[idx].c_str();
                    return true;
                }, static_cast<void*>(&focusOptions), (int)focusOptions.size()))
            {
                // cameraFocusIndex updated automatically by ImGui
            }

            // Calculate the current target position dynamically
            glm::vec3 target;
            if (cameraFocusIndex == 0) {
                target = calculateCenterOfMass(bodiesVector);
                ImGui::Text("Camera centered on: COM (%.2f, %.2f, %.2f)", target.x, target.y, target.z);
            }
            else {
                const auto& body = bodiesVector[cameraFocusIndex - 1];
                target = glm::vec3(body.x, body.y, body.z);
                ImGui::Text("Camera centered on: %s (%.2f, %.2f, %.2f)",
                    body.BodyName.c_str(), target.x, target.y, target.z);
            }
        }

        ImGui::End();

        // --- Handle Selection Change ---
        if (selectedIndex != previousIndex && !systems.empty()) {
            bodiesVector = fetchCelestialBodies(conn, systems[selectedIndex].first);
            initializeOrbitalVelocities(bodiesVector);
            previousIndex = selectedIndex;
        }

        // --- Simulation & Rendering ---
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(zoomLevel, 640.0f / 480.0f, 0.1f, 1000.0f);

        // --- Camera Centering Logic ---
        glm::vec3 focusPos;
        if (cameraFocusIndex == 0) {
            focusPos = calculateCenterOfMass(bodiesVector);
        }
        else {
            const auto& body = bodiesVector[cameraFocusIndex - 1];
            focusPos = glm::vec3(body.x, body.y, body.z);
        }

        float distance = zoomLevel;
        float pitchRad = glm::radians(pitch);
        float yawRad = glm::radians(yaw);

        glm::vec3 cameraOffset;
        cameraOffset.x = distance * cos(pitchRad) * cos(yawRad);
        cameraOffset.y = distance * sin(pitchRad);
        cameraOffset.z = distance * cos(pitchRad) * sin(yawRad);

        glm::vec3 cameraPos = focusPos + cameraOffset;

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(
            cameraPos.x, cameraPos.y, cameraPos.z,
            focusPos.x, focusPos.y, focusPos.z,
            0.0f, 1.0f, 0.0f
        );

        // Now draw your scene as usual:
        drawGravityGrid(100, 0.25f, bodiesVector);
        updatePositionsAndVelocities(bodiesVector, dt);

        for (const auto& body : bodiesVector) {
            drawSphere(body.size, body.x, body.y, body.z, body.textureID);
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // --- Cleanup ---
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
