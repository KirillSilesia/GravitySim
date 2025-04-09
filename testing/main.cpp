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
#pragma managed(push, off)
#include <occi.h>
#include <oratypes.h>
#pragma managed(pop)
#define OCCI_DEMO_NEED_IMPORT
#define _CRT_SECURE_NO_DEPRECATE

using namespace std;
using namespace oracle::occi;

float zoomLevel = 45.0f; // Initial zoom level
const float G = 0.66743f; // Adjusted gravitational constant
float deltaTime = 500.0f; // Smaller time step for accurate simulation
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
    string BodyName;
    string BodyType;
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
void drawGravityGrid(int gridSize, float lineSpacing, vector<CelestialBody>& bodies) {
    int numLines = gridSize / lineSpacing;
    float gravityScale = 50.f; // Increased scale for better visualization

    glColor3f(0.5f, 0.5f, 0.5f);  // Set the color for the grid
    glBegin(GL_LINES);

    std::vector<std::vector<float>> heights(numLines + 1, std::vector<float>(numLines + 1, 0.0f));

    // Calculate the height values based on the gravitational force at each point
    for (int i = 0; i <= numLines; ++i) {
        for (int j = 0; j <= numLines; ++j) {
            float x = (i - numLines / 2) * lineSpacing;
            float z = (j - numLines / 2) * lineSpacing;
            float y = 0.0f;

            for (size_t b = 0; b < bodies.size(); ++b) {
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

// Mouse button callback function
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
glm::vec3 calculateCenterOfMass(vector<CelestialBody>& bodies) {
    glm::vec3 centerOfMass(0.0f);
    float totalMass = 0.0f;

    // Loop through all bodies to accumulate their mass and position
    for (size_t i = 0; i < bodies.size(); ++i) {
        totalMass += bodies[i].mass;
        centerOfMass += bodies[i].mass * glm::vec3(bodies[i].x, bodies[i].y, bodies[i].z);
    }

    // Avoid division by zero if totalMass is zero
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

// Function to update the positions and velocities of celestial bodies based on gravity
void updatePositionsAndVelocities(vector<CelestialBody>& bodies, float deltaTime) {
    // Temporary arrays to store acceleration values
    vector<glm::vec3> accelerations(bodies.size(), glm::vec3(0.0f));

    // Compute gravitational forces between all pairs of bodies
    for (size_t i = 0; i < bodies.size(); ++i) {
        CelestialBody& body1 = bodies[i];

        for (size_t j = 0; j < bodies.size(); ++j) {
            if (i != j) {
                CelestialBody& body2 = bodies[j];

                // Calculate the distance vector from body1 to body2
                glm::vec3 r_vec = glm::vec3(body2.x - body1.x, 0, body2.z - body1.z); // 2D distance ignoring y-axis
                float r2 = glm::dot(r_vec, r_vec);  // r squared

                if (r2 < 1e-6f) continue;  // Avoid division by zero or too small values

                float r_dist = sqrt(r2);  // Scalar distance between the bodies
                float force = G * body1.mass * body2.mass / r2;  // Gravitational force magnitude

                // Calculate acceleration due to the gravitational force (Newton's 2nd law: F = ma)
                glm::vec3 acceleration = (force / body1.mass) * glm::normalize(r_vec); // Normalize the direction

                // Add the acceleration from body2 to body1 to the current body's acceleration
                accelerations[i] += acceleration;
            }
        }
    }

    // Update velocities and positions based on the calculated accelerations
    for (size_t i = 0; i < bodies.size(); ++i) {
        CelestialBody& body = bodies[i];

        // Update velocity
        body.vx += accelerations[i].x * deltaTime;
        body.vz += accelerations[i].z * deltaTime;

        // Update position
        body.x += body.vx * deltaTime;
        body.z += body.vz * deltaTime;
    }
}

void initializeOrbitalVelocities(vector<CelestialBody>& bodies) {
    glm::vec3 CoM = calculateCenterOfMass(bodies);

    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].isStar) continue; // Skip stars

        float dx = bodies[i].x - CoM.x;
        float dz = bodies[i].z - CoM.z;
        float r = sqrt(dx * dx + dz * dz);

        float orbitalSpeed = sqrt(G * 100.0e10f / r); //Approximation using mass of main star

        bodies[i].vx = -dz / r * orbitalSpeed;
        bodies[i].vz = dx / r * orbitalSpeed;
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

// Function to fetch solar system names from the database
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

// Function to fetch celestial bodies by SolarSystemID
vector<CelestialBody> fetchCelestialBodies(Connection* conn, int solarSystemID) {
    vector<CelestialBody> celestialBodies;
    Statement* stmt = nullptr;
    ResultSet* rs = nullptr;

    try {
        stmt = conn->createStatement();
        stmt->setSQL("SELECT BodyName, BodyType, Mass, Size, Color, X, Y, Z, VX, VY, VZ, IsStar FROM CelestialBodies WHERE SolarSystemID = :1");
        stmt->setInt(1, solarSystemID);
        rs = stmt->executeQuery();

        while (rs->next()) {
            CelestialBody body;
            body.BodyName = rs->getString(1);
            body.BodyType = rs->getString(2);
            body.mass = rs->getDouble(3);
            body.size = rs->getFloat(4);
            // Parse color string to float array
            string colorStr = rs->getString(5);
            sscanf_s(colorStr.c_str(), "%f,%f,%f", &body.color[0], &body.color[1], &body.color[2]);
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

// Update the bodies array with the selected solar system
void updateBodies(Connection* conn, int solarSystemID, vector<CelestialBody>& bodies) {
    vector<CelestialBody> fetchedBodies = fetchCelestialBodies(conn, solarSystemID);
    bodies = fetchedBodies;
}

int main(void) {
    Environment* env = nullptr;
    Connection* conn = nullptr;
    Statement* stmt = nullptr;
    ResultSet* rs = nullptr;
    GLFWwindow* window = nullptr;

    vector<CelestialBody> bodiesVector;
    CelestialBody* bodies = nullptr;
    int numBodies = 0;

    try {
        // Database Connection
        string user = "msbd15";
        string password = "haslo2025";
        string connectString = "//155.158.112.45:1521/oltpstud"; // Replace with your values

        env = Environment::createEnvironment(Environment::DEFAULT);
        conn = env->createConnection(user, password, connectString);
        cout << "Connected to the database!" << endl;

        // Test Query
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
    }

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

    float lastTime = glfwGetTime();
    float deltaTime;

    // Solar system selection
    static int selectedSolarSystemIndex = 0;
    vector<string> solarSystemNames = fetchSolarSystemNames(conn); // Fetch names from database
    if (!solarSystemNames.empty()) {
        updateBodies(conn, selectedSolarSystemIndex + 1, bodiesVector);
        numBodies = bodiesVector.size();
        bodies = new CelestialBody[numBodies];
        for (int i = 0; i < numBodies; ++i) {
            bodies[i] = bodiesVector[i];
        }
        initializeOrbitalVelocities(bodiesVector);
    }

    static int previousSolarSystemIndex = -1;

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        float currentTime = glfwGetTime();
        deltaTime = (currentTime - lastTime) * 10;
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
        // Solar system selection
        if (!solarSystemNames.empty()) {
            const char* combo_preview_value = solarSystemNames[selectedSolarSystemIndex].c_str();
            if (ImGui::BeginCombo("Solar System", combo_preview_value)) {
                for (int n = 0; n < solarSystemNames.size(); n++) {
                    const bool is_selected = (selectedSolarSystemIndex == n);
                    if (ImGui::Selectable(solarSystemNames[n].c_str(), is_selected)) {
                        selectedSolarSystemIndex = n;
                    }
                    if (is_selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        }

        if (selectedSolarSystemIndex != previousSolarSystemIndex) {
            updateBodies(conn, selectedSolarSystemIndex + 1, bodiesVector); // +1 because SolarSystemID starts from 1
            previousSolarSystemIndex = selectedSolarSystemIndex;

            // Convert the vector to a C-style array
            delete[] bodies;  // Free the old array
            numBodies = bodiesVector.size();
            bodies = new CelestialBody[numBodies];
            for (int i = 0; i < numBodies; ++i) {
                bodies[i] = bodiesVector[i];
            }
            initializeOrbitalVelocities(bodiesVector);
        }

        ImGui::End();

        // Clear the frame before rendering
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(zoomLevel, 640.0f / 480.0f, 0.1f, 100.0f);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glTranslatef(cameraX, cameraY, cameraZ);
        glRotatef(pitch, 1.0f, 0.0f, 0.0f);
        glRotatef(yaw, 0.0f, 1.0f, 0.0f);

        drawGravityGrid(100, 0.25f, bodiesVector);

        updatePositionsAndVelocities(bodiesVector, deltaTime);

        for (int i = 0; i < numBodies; ++i) {
            drawSphere(bodies[i].size, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].color);
        }

        // Render ImGui
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // Cleanup
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
