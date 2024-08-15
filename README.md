# CSCI534_RayTracing_Term_Project
Finished May 2023. My final project for CSCI 534, Computational Geometry. Implements a simple ray tracer and spline related topics. 
This is my computational geometry course project. 
Running the code in an IDE will bring up two windows, one containing a ray-tracer and another containing a canvas to draw and display Bezier Curves. 

To operate the Bezier Curve Canvas, select each point button and click on the canvas to place that control point. Once each control point has been placed, hit draw to generate the corresponding Bezier curve. The clear button resets the canvas, but not the control points. The user can reclick the point buttons to replace control points in new locations in order to draw a different Bezier Curve. 

It should be noted the ray-tracer is slow to render a scene. By clicking one of the four x++, x--, y++, or y-- buttons, the camera can be rotated about the scene about the corresponding axis in 1 degree intervals. Reset camera resets the camera's orientation. The subdivide button currently serves no function. New meshes can be switched out for rendering in line 316, using commented out WingedEdgeMesh object in the code immediately above. 
