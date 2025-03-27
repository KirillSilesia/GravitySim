#include <windows.h>
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cmath>
#include <vector>

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

float zoomLevel = 45.0f; // Initial zoom level
const float G = 6.67430e-11f; // Gravitational constant
float deltaTime = 0.1f; // Time step

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

// Callback function for zooming using the scroll wheel
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    float zoomSpeed = 5.0f;
    zoomLevel -= yoffset * zoomSpeed;  // Adjust the zoom level based on scroll
    if (zoomLevel < 10.0f) zoomLevel = 10.0f;  // Set a minimum zoom level
    if (zoomLevel > 90.0f) zoomLevel = 90.0f;  // Set a maximum zoom level
}

// Variables to track mouse movement and camera rotation
double lastX = 320.0, lastY = 240.0; // Initial mouse position (center of window)
float pitch = 0.0f, yaw = -90.0f;     // Camera rotation angles
bool isDragging = false;               // Flag for left mouse button hold

// Function to handle mouse movement (camera rotation)
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (!isDragging) return;  // Only update rotation when dragging

    float xOffset = xpos - lastX;  // Calculate change in X position
    float yOffset = lastY - ypos;  // Calculate change in Y position (invert Y to match typical 3D controls)
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.1f;  // Rotation sensitivity
    xOffset *= sensitivity;
    yOffset *= sensitivity;

    // Update camera angles (yaw and pitch)
    yaw += xOffset;
    pitch += yOffset;

    // Constrain pitch to avoid gimbal lock
    if (pitch > 89.0f) pitch = 89.0f;
    if (pitch < -89.0f) pitch = -89.0f;
}

// Function to handle mouse button events
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            isDragging = true;
            glfwGetCursorPos(window, &lastX, &lastY);
        }
        if (action == GLFW_RELEASE) {
            isDragging = false;
        }
    }
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

// Function to update the positions of celestial bodies based on gravity
void updateBodies(CelestialBody* bodies, int numBodies) {
    int starCount = 0;

    // Count the number of stars
    for (int i = 0; i < numBodies; ++i) {
        if (bodies[i].isStar) {
            starCount++;
        }
    }

    for (int i = 0; i < numBodies; ++i) {
        CelestialBody& body = bodies[i];

        // Stars should only move if there are 2 or more stars
        if (body.isStar && starCount < 2) {
            continue;  // Single star remains static
        }

        float fx = 0.0f, fz = 0.0f; // Forces along X and Z axes

        for (int j = 0; j < numBodies; ++j) {
            if (i == j) continue;  // Skip self-interaction

            CelestialBody& other = bodies[j];

            // Planets should only be affected by stars
            if (!body.isStar && !other.isStar) continue;

            // Planets should not influence stars
            if (body.isStar && !other.isStar) continue;

            float dx = other.x - body.x;
            float dz = other.z - body.z;
            float distanceSquared = dx * dx + dz * dz + 1.0f;  // Avoid division by zero
            float distance = sqrt(distanceSquared);

            if (distance < 0.1f) continue;  // Prevent extreme forces at close range

            float force = (G * body.mass * other.mass) / distanceSquared;

            fx += force * dx / distance;
            fz += force * dz / distance;
        }

        // Apply acceleration
        float ax = fx / body.mass;
        float az = fz / body.mass;

        // Update velocity
        body.vx += ax * deltaTime;
        body.vz += az * deltaTime;

        // Update position (no update on Y-axis)
        body.x += body.vx * deltaTime;
        body.z += body.vz * deltaTime;
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
    glfwSetFramebufferSizeCallback(window, reshape);

    // Test objects (celestial bodies)
    CelestialBody star1 = { 0.0f, 0.0f, 0.0f, 5.0f, 100.0e10f, {1.0f, 1.0f, 0.0f}, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, true };
CelestialBody star2 = { 10.0f, 0.0f, 0.0f, 0.5f, 2.0e10f, {1.0f, 0.0f, 0.0f}, -0.2f, 0.0f, 0.1f, 0.0f, 0.0f, true };
CelestialBody planet = { 5.0f, 0.0f, 0.0f, 0.3f, 5.0e9f, {0.0f, 0.0f, 1.0f}, 0.0f, 0.0f, 0.0f, false };


    CelestialBody bodies[] = { star1, star2, planet };

    float cameraX = 0.0f, cameraY = 0.0f, cameraZ = -15.0f;

    glfwSetScrollCallback(window, scroll_callback);

    // Enable depth testing for proper 3D rendering
    glEnable(GL_DEPTH_TEST);

    reshape(window, 640, 480);

    // Set the perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 640.0 / 480.0, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);

    // Set GLFW mouse callback functions
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear the screen

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(zoomLevel, 640.0f / 480.0f, 0.1f, 100.0f);

        // Apply camera transformations (translation and rotation)
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        updateBodies(bodies, 3);  // Update the positions of celestial bodies

        glTranslatef(0.0f, 0.0f, cameraZ); // Apply camera translation
        glRotatef(pitch, 1.0f, 0.0f, 0.0f); // Apply pitch rotation
        glRotatef(yaw, 0.0f, 1.0f, 0.0f);   // Apply yaw rotation

        // Draw the gravity grid for visualization
        drawGravityGrid(50, 0.25f, bodies, 3);

        // Draw celestial bodies
        for (int i = 0; i < 3; ++i) {
            drawSphere(bodies[i].size, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].color);
        }


        updateBodies(bodies, 3);  // Update the positions of celestial bodies

        glfwSwapBuffers(window);  // Swap buffers
        glfwPollEvents();  // Poll for events
    }

    glfwTerminate();  // Clean up GLFW resources
    return 0;
}
