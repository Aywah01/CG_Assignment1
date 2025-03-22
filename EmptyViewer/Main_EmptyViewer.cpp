//#include <Windows.h>
//#include <iostream>
//#include <GL/glew.h>
//#include <GL/GL.h>
//#include <GL/freeglut.h>
//
//#define GLFW_INCLUDE_GLU
//#define GLFW_DLL
//#include <GLFW/glfw3.h>
//#include <vector>
//
//#define GLM_SWIZZLE
//#include <glm/glm.hpp>
//#include <glm/gtc/constants.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtx/string_cast.hpp>
//
//using namespace glm;
//
//// -------------------------------------------------
//// Global Variables
//// -------------------------------------------------
//int Width = 1280;
//int Height = 720;
//std::vector<float> OutputImage;
//// -------------------------------------------------
//
//
//void render()
//{
//	//Create our image. We don't want to do this in 
//	//the main loop since this may be too slow and we 
//	//want a responsive display of our beautiful image.
//	//Instead we draw to another buffer and copy this to the 
//	//framebuffer using glDrawPixels(...) every refresh
//	OutputImage.clear();
//	for (int j = 0; j < Height; ++j)
//	{
//		for (int i = 0; i < Width; ++i)
//		{
//			// ---------------------------------------------------
//			// --- Implement your code here to generate the image
//			// ---------------------------------------------------
//
//			// draw a red rectangle in the center of the image
//			vec3 color = glm::vec3(0.5f, 0.5f, 0.5f); // grey color [0,1] in RGB channel
//
//			if (i > Width / 4 && i < 3 * Width / 4
//				&& j > Height / 4 && j < 3 * Height / 4)
//			{
//				color = glm::vec3(1.0f, 1.0f, 0.0f); // red color [0,1] in RGB channel
//			}
//
//			// set the color
//			OutputImage.push_back(color.x); // R
//			OutputImage.push_back(color.y); // G
//			OutputImage.push_back(color.z); // B
//		}
//	}
//}
//
//void resize_callback(GLFWwindow*, int nw, int nh)
//{
//	//This is called in response to the window resizing.
//	//The new width and height are passed in so we make 
//	//any necessary changes:
//	Width = nw;
//	Height = nh;
//	//Tell the viewport to use all of our screen estate
//	glViewport(0, 0, nw, nh);
//
//	//This is not necessary, we're just working in 2d so
//	//why not let our spaces reflect it?
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//
//	glOrtho(0.0, static_cast<double>(Width)
//		, 0.0, static_cast<double>(Height)
//		, 1.0, -1.0);
//
//	//Reserve memory for our render so that we don't do 
//	//excessive allocations and render the image
//	OutputImage.reserve(Width * Height * 3);
//	render();
//}
//
//
//int main(int argc, char* argv[])
//{
//	// -------------------------------------------------
//	// Initialize Window
//	// -------------------------------------------------
//
//	GLFWwindow* window;
//
//	/* Initialize the library */
//	if (!glfwInit())
//		return -1;
//
//	/* Create a windowed mode window and its OpenGL context */
//	window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
//	if (!window)
//	{
//		glfwTerminate();
//		return -1;
//	}
//
//	/* Make the window's context current */
//	glfwMakeContextCurrent(window);
//
//	//We have an opengl context now. Everything from here on out 
//	//is just managing our window or opengl directly.
//
//	//Tell the opengl state machine we don't want it to make 
//	//any assumptions about how pixels are aligned in memory 
//	//during transfers between host and device (like glDrawPixels(...) )
//	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
//	glPixelStorei(GL_PACK_ALIGNMENT, 1);
//
//	//We call our resize function once to set everything up initially
//	//after registering it as a callback with glfw
//	glfwSetFramebufferSizeCallback(window, resize_callback);
//	resize_callback(NULL, Width, Height);
//
//	/* Loop until the user closes the window */
//	while (!glfwWindowShouldClose(window))
//	{
//		//Clear the screen
//		glClear(GL_COLOR_BUFFER_BIT);
//
//		// -------------------------------------------------------------
//		//Rendering begins!
//		glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
//		//and ends.
//		// -------------------------------------------------------------
//
//		/* Swap front and back buffers */
//		glfwSwapBuffers(window);
//
//		/* Poll for and process events */
//		glfwPollEvents();
//
//		//Close when the user hits 'q' or escape
//		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
//			|| glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
//		{
//			glfwSetWindowShouldClose(window, GL_TRUE);
//		}
//	}
//
//	glfwDestroyWindow(window);
//	glfwTerminate();
//	return 0;
//}

//#include <Windows.h>
//#include <iostream>
//#include <GL/glew.h>
//#include <GL/GL.h>
//#include <GL/freeglut.h>
//#include <GLFW/glfw3.h>
//#include <cmath>
//#include <vector>
//#include <memory>
//#include <limits>
//
//#define GLM_SWIZZLE
//#include <glm/glm.hpp>
//#include <glm/gtc/constants.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtx/string_cast.hpp>
//
//using namespace glm;
//
//// -------------------------------------------------
//// Global Variables
//// -------------------------------------------------
//int Width = 512; // Image width
//int Height = 512; // Image height
//std::vector<float> OutputImage;
//
//// Class Definitions
//class Ray {
//public:
//    vec3 origin; // Ray origin
//    vec3 direction; // Ray direction
//
//    Ray(const vec3& origin, const vec3& direction)
//        : origin(origin), direction(direction) {
//    }
//};
//
//class Camera {
//public:
//    vec3 eye; // Camera position
//    vec3 u, v, w; // Camera orientation (basis vectors)
//
//    Camera(const vec3& eye, const vec3& u, const vec3& v, const vec3& w)
//        : eye(eye), u(u), v(v), w(w) {
//    }
//
//    Ray getRay(int ix, int iy) {
//        // Image plane corners
//        float l = -0.1, r = 0.1, b = -0.1, t = 0.1, d = 0.1;
//
//        // Calculate image plane coordinates
//        float x = l + (r - l) * (ix + 0.5) / Width;
//        float y = b + (t - b) * (iy + 0.5) / Height;
//
//        // Projection direction
//        vec3 dir = normalize(x * u + y * v - d * w);
//        return Ray(eye, dir);
//    }
//};
//
//class Surface {
//public:
//    virtual bool intersect(const Ray& ray, float& t) = 0; // Pure virtual function for intersection
//};
//
//class Plane : public Surface {
//public:
//    float y; // height of the plane
//
//    Plane(float y) : y(y) {}
//
//    bool intersect(const Ray& ray, float& t) override {
//        if (abs(ray.direction.y) < 1e-6) { // Ray is parallel to the plane
//            return false;
//        }
//        t = (y - ray.origin.y) / ray.direction.y; // t is the distance to the intersection
//        return (t >= 0); // Intersection must be in the positive direction
//    }
//};
//
//class Sphere : public Surface {
//public:
//    vec3 center; // Center of the sphere
//    float radius; // Radius of the sphere
//
//    Sphere(const vec3& center, float radius) : center(center), radius(radius) {}
//
//    bool intersect(const Ray& ray, float& t) override {
//        vec3 oc = ray.origin - center;
//        float a = dot(ray.direction, ray.direction);
//        float b = 2.0f * dot(oc, ray.direction);
//        float c = dot(oc, oc) - radius * radius;
//        float discriminant = b * b - 4.0f * a * c;
//
//        if (discriminant < 0) {
//            return false; // No intersection
//        }
//
//        t = (-b - std::sqrt(discriminant)) / (2.0f * a); // Nearest intersection
//        return (t >= 0); // Intersection must be in the positive direction
//    }
//};
//
//class Scene {
//public:
//    std::vector<std::shared_ptr<Surface>> surfaces;
//
//    void addSurface(const std::shared_ptr<Surface>& surface) {
//        surfaces.push_back(surface);
//    }
//
//    vec3 trace(const Ray& ray) {
//        float closestT = std::numeric_limits<float>::infinity();
//        std::shared_ptr<Surface> closestSurface = nullptr;
//
//        for (const auto& surface : surfaces) {
//            float t;
//            if (surface->intersect(ray, t) && t < closestT) {
//                closestT = t;
//                closestSurface = surface;
//            }
//        }
//
//        // If any surface is hit, return white color; else return black
//        return (closestSurface != nullptr) ? vec3(1.0f) : vec3(0.0f);
//    }
//};
//
//// Render Function
//void render() {
//    OutputImage.clear();
//    Camera camera(vec3(0.0f, 0.0f, 0.0f), // Eye position
//        vec3(1.0f, 0.0f, 0.0f), // u direction
//        vec3(0.0f, 1.0f, 0.0f), // v direction
//        vec3(0.0f, 0.0f, -1.0f)); // w direction
//
//    // Creating the scene with a plane and three spheres
//    Scene scene;
//    scene.addSurface(std::make_shared<Plane>(-2.0f)); // Plane y = -2
//    scene.addSurface(std::make_shared<Sphere>(vec3(-4.0f, 0.0f, -7.0f), 1.0f)); // Sphere S1
//    scene.addSurface(std::make_shared<Sphere>(vec3(0.0f, 0.0f, -7.0f), 2.0f)); // Sphere S2
//    scene.addSurface(std::make_shared<Sphere>(vec3(4.0f, 0.0f, -7.0f), 1.0f)); // Sphere S3
//
//    for (int iy = 0; iy < Height; ++iy) {
//        for (int ix = 0; ix < Width; ++ix) {
//            Ray ray = camera.getRay(ix, iy);
//            vec3 color = scene.trace(ray);
//
//            // Set pixel color in the OutputImage
//            OutputImage.push_back(color.r);
//            OutputImage.push_back(color.g);
//            OutputImage.push_back(color.b);
//        }
//    }
//
//    // Render OutputImage as pixels
//    glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
//}
//
//// Resize Callback Function
//void resize_callback(GLFWwindow*, int nw, int nh) {
//    Width = nw;
//    Height = nh;
//    glViewport(0, 0, nw, nh);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//
//    // Use a perspective projection
//    gluPerspective(45.0, (double)Width / (double)Height, 0.1, 100.0);
//    glMatrixMode(GL_MODELVIEW);
//}
//
//// Main Function
//int main(int argc, char* argv[]) {
//    // Initialize Window
//    GLFWwindow* window;
//
//    // Initialize the library 
//    if (!glfwInit())
//        return -1;
//
//    // Create a windowed mode window and its OpenGL context
//    window = glfwCreateWindow(Width, Height, "Ray Tracer", NULL, NULL);
//    if (!window) {
//        glfwTerminate();
//        return -1;
//    }
//
//    // Make the window's context current
//    glfwMakeContextCurrent(window);
//
//    // Enable depth test
//    glEnable(GL_DEPTH_TEST);
//
//    // Set clear color for the window
//    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
//
//    // Set viewport and projection matrix
//    glfwSetFramebufferSizeCallback(window, resize_callback);
//    resize_callback(nullptr, Width, Height); // Initial setup
//
//    // Main Loop
//    while (!glfwWindowShouldClose(window)) {
//        render(); // Render the scene
//
//        // Swap front and back buffers
//        glfwSwapBuffers(window);
//
//        // Poll for and process events
//        glfwPollEvents();
//
//        // Close when the user hits 'q' or escape
//        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS ||
//            glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
//            glfwSetWindowShouldClose(window, GL_TRUE);
//        }
//    }
//
//    glfwDestroyWindow(window);
//    glfwTerminate();
//    return 0;
//}

//latest code 
#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
int Width = 1280;
int Height = 720;
std::vector<float> OutputImage;

// Function to draw a sphere with a given color
void drawSphere(float x, float y, float z, float radius, vec3 color) {
    const int latSteps = 20; // Number of latitude steps
    const int longSteps = 20; // Number of longitude steps

    for (int i = 0; i < latSteps; ++i) {
        float lat0 = glm::pi<float>() * (-0.5 + (float)(i) / latSteps); // getting latitude
        float z0 = radius * sin(lat0); // height of the bottom circle
        float r0 = radius * cos(lat0); // radius of the bottom circle

        float lat1 = glm::pi<float>() * (-0.5 + (float)(i + 1) / latSteps); // getting next latitude
        float z1 = radius * sin(lat1); // height of the top circle
        float r1 = radius * cos(lat1); // radius of the top circle

        glBegin(GL_TRIANGLE_STRIP); // Use TRIANGLE_STRIP to draw the sphere
        for (int j = 0; j <= longSteps; ++j) {
            float lng = 2 * glm::pi<float>() * (float)(j) / longSteps; // getting longitude
            float x0 = r0 * cos(lng); // X coordinate for current latitude
            float y0 = r0 * sin(lng); // Y coordinate for current latitude
            float x1 = r1 * cos(lng); // X coordinate for the next latitude
            float y1 = r1 * sin(lng); // Y coordinate for the next latitude

            // Set color for each vertex
            glColor3f(color.x, color.y, color.z);

            // First triangle vertex
            glVertex3f(x + x0, y + z0, z + y0);
            // Second triangle vertex
            glVertex3f(x + x1, y + z1, z + y1);
        }
        glEnd();
    }
}

// Function to draw a plane
void drawPlane(float width, float height, float z) {
    glBegin(GL_QUADS); // Draw a simple rectangle as a plane
    glColor3f(0.5f, 0.5f, 0.5f); // Grey color
    glVertex3f(-width / 2, z, -height / 2); // Bottom-left
    glVertex3f(width / 2, z, -height / 2);  // Bottom-right
    glVertex3f(width / 2, z, height / 2);   // Top-right
    glVertex3f(-width / 2, z, height / 2);  // Top-left
    glEnd();
}

void render() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // Move the camera back and slightly up
    gluLookAt(0.0, 2.0, 10.0,   // Eye position (x, y, z)
        0.0, 0.0, 0.0,    // Look-at position (x, y, z)
        0.0, 1.0, 0.0);   // Up direction (x, y, z)

    // Draw the plane
    drawPlane(10.0f, 10.0f, 0.0f); // Z-position of the plane

    // Draw three spheres at different depths
    drawSphere(-2.0f, 1.0f, -3.0f, 0.5f, vec3(1.0f, 0.0f, 0.0f)); // Red sphere
    drawSphere(0.0f, 1.0f, 0.0f, 0.5f, vec3(0.0f, 1.0f, 0.0f)); // Green sphere
    drawSphere(2.0f, 1.0f, 3.0f, 0.5f, vec3(0.0f, 0.0f, 1.0f)); // Blue sphere
}

void resize_callback(GLFWwindow*, int nw, int nh) {
    Width = nw;
    Height = nh;
    glViewport(0, 0, nw, nh);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Use a perspective projection
    gluPerspective(45.0, (double)Width / (double)Height, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char* argv[]) {
    // Initialize Window
    GLFWwindow* window;

    // Initialize the library 
    if (!glfwInit())
        return -1;

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);

    // Set clear color for the window
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background

    // Set viewport and projection matrix
    glfwSetFramebufferSizeCallback(window, resize_callback);
    resize_callback(nullptr, Width, Height); // Initial setup

    // Main Loop
    while (!glfwWindowShouldClose(window)) {
        render(); // Render the scene

        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();

        // Close when the user hits 'q' or escape
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <limits>
//#include <memory>
//#include <GLFW/glfw3.h>
//
//const int WIDTH = 512;
//const int HEIGHT = 512;
//const float INF = std::numeric_limits<float>::max();
//const float EPSILON = 0.001f;
//
//// Define a vector class
//struct Vec3 {
//    float x, y, z;
//
//    Vec3() : x(0), y(0), z(0) {}
//    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
//
//    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
//    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
//    Vec3 operator*(float t) const { return Vec3(x * t, y * t, z * t); }
//    Vec3 operator/(float t) const { return Vec3(x / t, y / t, z / t); }
//
//    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
//    Vec3 normalize() const {
//        float len = std::sqrt(dot(*this));
//        return Vec3(x / len, y / len, z / len);
//    }
//};
//
//class Ray {
//public:
//    Vec3 origin, direction;
//
//    Ray(const Vec3& origin, const Vec3& direction) : origin(origin), direction(direction.normalize()) {}
//};
//
//class Surface {
//public:
//    virtual ~Surface() = default;
//    virtual bool intersect(const Ray& ray, float& t) const = 0;
//};
//
//class Plane : public Surface {
//public:
//    float y;
//
//    Plane(float y) : y(y) {}
//
//    bool intersect(const Ray& ray, float& t) const override {
//        if (std::fabs(ray.direction.y) < EPSILON) return false; // Parallel to plane
//        t = (y - ray.origin.y) / ray.direction.y;
//        return (t >= 0);
//    }
//};
//
//class Sphere : public Surface {
//public:
//    Vec3 center;
//    float radius;
//
//    Sphere(const Vec3& center, float radius) : center(center), radius(radius) {}
//
//    bool intersect(const Ray& ray, float& t) const override {
//        Vec3 oc = ray.origin - center;
//        float a = ray.direction.dot(ray.direction);
//        float b = 2.0f * oc.dot(ray.direction);
//        float c = oc.dot(oc) - radius * radius;
//        float discriminant = b * b - 4 * a * c;
//
//        if (discriminant < 0) return false;
//        t = (-b - std::sqrt(discriminant)) / (2.0f * a);
//        return (t >= 0);
//    }
//};
//
//class Camera {
//public:
//    Vec3 position, u, v, w;
//    float l, r, b, t, d;
//
//    Camera(const Vec3& position) : position(position) {
//        u = Vec3(1, 0, 0);
//        v = Vec3(0, 1, 0);
//        w = Vec3(0, 0, -1);
//        l = -0.1f; r = 0.1f; b = -0.1f; t = 0.1f; d = 0.1f;
//    }
//
//    Ray getRay(int ix, int iy) const {
//        float x = l + (r - l) * (ix + 0.5f) / WIDTH;
//        float y = b + (t - b) * (iy + 0.5f) / HEIGHT;
//        Vec3 ray_dir = u * x + v * y - w * d;
//        return Ray(position, ray_dir);
//    }
//};
//
//class Scene {
//public:
//    std::vector<std::shared_ptr<Surface>> surfaces;
//
//    void addSurface(const std::shared_ptr<Surface>& surface) {
//        surfaces.push_back(surface);
//    }
//
//    Vec3 trace(const Ray& ray) const {
//        float closest_t = INF;
//        std::shared_ptr<Surface> closest_surface = nullptr;
//
//        for (const auto& surface : surfaces) {
//            float t;
//            if (surface->intersect(ray, t) && t < closest_t) {
//                closest_t = t;
//                closest_surface = surface;
//            }
//        }
//
//        if (closest_surface) {
//            return Vec3(1.0f, 1.0f, 1.0f); // Color white for intersections
//        }
//
//        return Vec3(0.0f, 0.0f, 0.0f); // Color black for no intersections
//    }
//};
//
//// Function to render the scene, filling the image buffer
//void renderScene(const Camera& camera, Scene& scene, Vec3* image) {
//    for (int iy = 0; iy < HEIGHT; iy++) {
//        for (int ix = 0; ix < WIDTH; ix++) {
//            Ray ray = camera.getRay(ix, iy);
//            image[iy * WIDTH + ix] = scene.trace(ray);
//        }
//    }
//}
//
//// Function to convert Vec3 to pixel value for OpenGL
//template <typename T>
//T clamp(T value, T min, T max) {
//    if (value < min) return min;
//    if (value > max) return max;
//    return value;
//}
//
//void setPixelColor(unsigned char* buffer, int x, int y, const Vec3& color) {
//    int idx = (y * WIDTH + x) * 3;
//    buffer[idx] = static_cast<unsigned char>(clamp(color.x * 255, 0.0f, 255.0f));
//    buffer[idx + 1] = static_cast<unsigned char>(clamp(color.y * 255, 0.0f, 255.0f));
//    buffer[idx + 2] = static_cast<unsigned char>(clamp(color.z * 255, 0.0f, 255.0f));
//}
//
//void displayImage(GLFWwindow* window, unsigned char* image) {
//    glClear(GL_COLOR_BUFFER_BIT);
//    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, image);
//    glfwSwapBuffers(window);
//}
//
//int main() {
//    // Initialize GLFW
//    if (!glfwInit()) {
//        return -1;
//    }
//
//    // Create a windowed mode window and its OpenGL context
//    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Ray Tracer", nullptr, nullptr);
//    if (!window) {
//        glfwTerminate();
//        return -1;
//    }
//
//    glfwMakeContextCurrent(window);
//
//    // Set up the scene
//    Scene scene;
//    scene.addSurface(std::make_shared<Plane>(-2.0f));
//    scene.addSurface(std::make_shared<Sphere>(Vec3(-4, 0, -7), 1.0f));
//    scene.addSurface(std::make_shared<Sphere>(Vec3(0, 0, -7), 2.0f));
//    scene.addSurface(std::make_shared<Sphere>(Vec3(4, 0, -7), 1.0f));
//
//    // Set up the camera
//    Camera camera(Vec3(0, 0, 0));
//
//    // Create an image buffer for rendering
//    unsigned char* image = new unsigned char[WIDTH * HEIGHT * 3];
//
//    // Render the scene
//    renderScene(camera, scene, reinterpret_cast<Vec3*>(image));
//
//    // Display the image in the window
//    while (!glfwWindowShouldClose(window)) {
//        displayImage(window, image);
//
//        // Poll for and process events
//        glfwPollEvents();
//    }
//
//    // Clean up memory
//    delete[] image;
//    glfwDestroyWindow(window);
//    glfwTerminate();
//    return 0;
//}
 
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <limits>
//#include <memory>
//#include <algorithm>
//#define _CRT_SECURE_NO_WARNINGS
//
//const int WIDTH = 512;
//const int HEIGHT = 512;
//const float INF = std::numeric_limits<float>::max();
//const float EPSILON = 0.001f;
//
//struct Vec3 {
//    float x, y, z;
//
//    Vec3() : x(0), y(0), z(0) {}
//    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
//
//    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
//    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
//    Vec3 operator*(float t) const { return Vec3(x * t, y * t, z * t); }
//    Vec3 operator/(float t) const { return Vec3(x / t, y / t, z / t); }
//
//    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
//    Vec3 normalize() const {
//        float len = sqrt(dot(*this));
//        return Vec3(x / len, y / len, z / len);
//    }
//};
//
//class Ray {
//public:
//    Vec3 origin, direction;
//
//    Ray(const Vec3& origin, const Vec3& direction) : origin(origin), direction(direction.normalize()) {}
//};
//
//class Surface {
//public:
//    virtual ~Surface() = default;
//    virtual bool intersect(const Ray& ray, float& t) const = 0;
//};
//
//class Plane : public Surface {
//public:
//    float y;
//
//    Plane(float y) : y(y) {}
//
//    bool intersect(const Ray& ray, float& t) const override {
//        if (fabs(ray.direction.y) < EPSILON) return false; // Parallel to plane
//        t = (y - ray.origin.y) / ray.direction.y;
//        return (t >= 0);
//    }
//};
//
//class Sphere : public Surface {
//public:
//    Vec3 center;
//    float radius;
//
//    Sphere(const Vec3& center, float radius) : center(center), radius(radius) {}
//
//    bool intersect(const Ray& ray, float& t) const override {
//        Vec3 oc = ray.origin - center;
//        float a = ray.direction.dot(ray.direction);
//        float b = 2.0f * oc.dot(ray.direction);
//        float c = oc.dot(oc) - radius * radius;
//        float discriminant = b * b - 4 * a * c;
//
//        if (discriminant < 0) return false;
//        t = (-b - sqrt(discriminant)) / (2.0f * a);
//        return (t >= 0);
//    }
//};
//
//class Camera {
//public:
//    Vec3 position, u, v, w;
//    float l, r, b, t, d;
//
//    Camera(const Vec3& position) : position(position) {
//        u = Vec3(1, 0, 0);
//        v = Vec3(0, 1, 0);
//        w = Vec3(0, 0, -1);
//        l = -0.1f; r = 0.1f; b = -0.1f; t = 0.1f; d = 0.1f;
//    }
//
//    Ray getRay(int ix, int iy) const {
//        float x = l + (r - l) * (ix + 0.5f) / WIDTH;
//        float y = b + (t - b) * (iy + 0.5f) / HEIGHT;
//        Vec3 ray_dir = u * x + v * y - w * d;
//        return Ray(position, ray_dir);
//    }
//};
//
//class Scene {
//public:
//    std::vector<std::shared_ptr<Surface>> surfaces;
//
//    void addSurface(const std::shared_ptr<Surface>& surface) {
//        surfaces.push_back(surface);
//    }
//
//    Vec3 trace(const Ray& ray) const {
//        float closest_t = INF;
//        std::shared_ptr<Surface> closest_surface = nullptr;
//
//        for (const auto& surface : surfaces) {
//            float t;
//            if (surface->intersect(ray, t) && t < closest_t) {
//                closest_t = t;
//                closest_surface = surface;
//            }
//        }
//
//        if (closest_surface) {
//            return Vec3(1.0f, 1.0f, 1.0f); //
//            // Color the pixel white if an intersection occurs
//        }
//
//        // Color the pixel black if no intersection
//        return Vec3(0.0f, 0.0f, 0.0f);
//    }
//};
//
//// Function to render the scene
//void renderScene(const Camera& camera, Scene& scene, Vec3 image[HEIGHT][WIDTH]) {
//    for (int iy = 0; iy < HEIGHT; iy++) {
//        for (int ix = 0; ix < WIDTH; ix++) {
//            Ray ray = camera.getRay(ix, iy);
//            image[iy][ix] = scene.trace(ray);
//        }
//    }
//}
//
//void saveImage(const Vec3 image[HEIGHT][WIDTH], const char* filename) {
//    FILE* file;
//    errno_t err = fopen_s(&file, filename, "w");
//    if (err != 0) {
//        std::cerr << "Error opening file: " << filename << std::endl;
//        return;
//    }
//    fprintf(file, "P3\n%d %d\n255\n", WIDTH, HEIGHT);
//    for (int iy = 0; iy < HEIGHT; iy++) {
//        for (int ix = 0; ix < WIDTH; ix++) {
//            int r = static_cast<int>(image[iy][ix].x * 255);
//            int g = static_cast<int>(image[iy][ix].y * 255);
//            int b = static_cast<int>(image[iy][ix].z * 255);
//            fprintf(file, "%d %d %d ", r, g, b);
//        }
//        fprintf(file, "\n");
//    }
//    fclose(file);
//}
//
//// Function to display the rendered image in a simple manner
//void displayImage(const Vec3 image[HEIGHT][WIDTH]) {
//    for (int iy = 0; iy < HEIGHT; iy++) {
//        for (int ix = 0; ix < WIDTH; ix++) {
//            // Print a character based on brightness
//            float brightness = (image[iy][ix].x + image[iy][ix].y + image[iy][ix].z) / 3.0f;
//            char pixel = (brightness > 0.5f) ? '#' : ' ';
//            std::cout << pixel;
//        }
//        std::cout << std::endl;
//    }
//}
//
//int main() {
//    // Initialize the scene
//    Scene scene;
//
//    // Create a plane at y = -2
//    scene.addSurface(std::make_shared<Plane>(-2.0f));
//
//    // Create spheres
//    scene.addSurface(std::make_shared<Sphere>(Vec3(-4, 0, -7), 1.0f));
//    scene.addSurface(std::make_shared<Sphere>(Vec3(0, 0, -7), 2.0f));
//    scene.addSurface(std::make_shared<Sphere>(Vec3(4, 0, -7), 1.0f));
//
//    // Set up the camera
//    Camera camera(Vec3(0, 0, 0));
//
//    // Create an image buffer to hold pixel colors
//    Vec3(*image)[WIDTH] = new Vec3[HEIGHT][WIDTH];
//
//    // Render the scene
//    renderScene(camera, scene, image);
//
//    // Save the rendered image to a file
//    saveImage(image, "output.ppm");
//
//    // Display the rendered image (in a simple text format)
//    displayImage(image);
//
//    // Free the allocated memory
//    delete[] image;
//
//    return 0;
//}