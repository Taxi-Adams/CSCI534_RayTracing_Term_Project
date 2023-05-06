import javax.swing.JPanel;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.HashMap;

public class RayTracer{
        private double vW; // viewport width
        private double vH; // viewport height
        private int d; // distance to viewport
        private int cW; // canvas width
        private int cH; // canvas height
        private Light[] lights;

    public RayTracer() {
        vW = 1.77; // viewport width
        vH = 1.0; // viewport height
        d = 1; // distance to viewport
        cW = 880; // canvas width
        cH = 495; // canvas height - rendering image at 495p
        lights = new Light[3];
        lights[0] = null;
        lights[1] = null;
        lights[2] = null;
    }

    public double getVH() {
        return vH;
    }

    public double getVW() {
        return vW;
    }

    public int getD() {
        return d;
    }

    public int getCH() {
        return cH;
    }

    public int getCW() {
        return cW;
    }

    public Light[] getLights() {
        return lights;
    }

    public void traceScene(Camera camera, Sphere[] spheres, int recursionDepth, PixelPainter canvas, double[] xyzAngle, int point1, int point2, WingedEdgeMesh[] meshes) {
        // int i = -1*threadNumber*cW/numThreads - 10; i < threadNumber*cW/numThreads; i++
        for(int i = point1; i < point2+1; i++) {
            for(int j = -1*cH/2; j < cH/2; j++) {
                double[] dVector = canvasToViewport(i, j); // D vector
                dVector = rotateCamera(dVector, xyzAngle); // rotates camera by given degrees for each axis
                short[] color = traceRay(camera.position, dVector, 1.0, Double.MAX_VALUE, spheres, recursionDepth, meshes);
                Color canvColor = new Color(color[0], color[1], color[2]);
                // paint pixel to canvas
                canvas.putPixel(canvas.getGraphics(), convertedX(i), convertedY(j), canvColor);
            }
        }
    }

    // public void traceSceneInThread(Camera camera, Sphere[] spheres, int recursionDepth, PixelPainter canvas, double[] xyzAngle, int numThreads, int threadNumber) {
    //     for(int i = -1*cW/2-10; i < cW/2; i++) {
    //         for(int j = -1*cH/2; j < cH/2; j++) {
    //             double[] dVector = canvasToViewport(i, j); // D vector
    //             dVector = rotateCamera(dVector, xyzAngle);
    //             // System.out.println("D[0]: " + dVector[0] + ", D[1]: " + dVector[1] + ", D[2]: " + dVector[2]);
    //             int[] color = traceRay(camera.position, dVector, 1.0, Double.MAX_VALUE, spheres, recursionDepth);
    //             Color canvColor = new Color(color[0], color[1], color[2]);
    //             // paint pixel to canvas
    //             canvas.putPixel(canvas.getGraphics(), convertedX(i), convertedY(j), canvColor);
    //         }
    //     }
    // }

    // multiplies vector by a given scalar
    public double[] multVectorByScalar(double[] vector, double scalar) {
        double[] newVector = {vector[0]*scalar, vector[1]*scalar, vector[2]*scalar};
        return newVector;
    }

    public int convertedX(int x) {
        return cW/2 + x;
    }

    public int convertedY(int y) {
        return cH/2 - y;
    }

    public double calcVx(int cX) {
        return ((double)cX)*vW/cW;
    }

    public double calcVy(int cY) {
        return ((double)cY)*vH/cH;
    }

    public double[] reflectRay(double[] r, double[] n) {
        double[] ray = {2*n[0]*dot(n, r) - r[0], 2*n[1]*dot(n, r) - r[1], 2*n[2]*dot(n, r) - r[2]};
        return ray;
    }

    // converts location on canvas grid to location on the viewport
    public double[] canvasToViewport(int x, int y) {
        double[] arr = {calcVx(x), calcVy(y), d};
        return arr;
    }

    public short[] traceRay(double[] oVector, double[] dVector, double tMin, double tMax, Sphere[] spheres, int recursionDepth, WingedEdgeMesh[] meshes) {
        // Variables
        short[] sphereColor = {0, 0, 0};
        short[] meshColor = {0, 0, 0};
        HashMap<Sphere, Double> closestVals = closestIntersection(oVector, dVector, tMin, tMax, spheres);
        HashMap<WingedEdgeMesh, double[]> closestMeshVals = closestFVIntersection(oVector, dVector, tMin, tMax, meshes);
        double closest_t = Double.MAX_VALUE;
        double closest_St = Double.MAX_VALUE;
        double closest_Mt = Double.MAX_VALUE;
        Sphere closestSphere = null;
        WingedEdgeMesh closestMesh = null;
        boolean sphereIsClosest = true;
        double sphereDistance = Double.MAX_VALUE;
        double meshDistance = Double.MAX_VALUE;
        double[] pVector = new double[3];
        double[] nVector = {0.0, 0.0, 0.0};
        short[] localColor = {0, 0, 0};
        double[] negativeDVector = {-1*dVector[0], -1*dVector[1], -1*dVector[2]};
        double nLength = 0.0;
        double rflctv = 0.0;
        double iVal = 0.0;

        // closestVals and closestMeshVals should have only one entry, HashMaps used for simplicity of returning two different data types
        // determines whether the a sphere or a mesh is closer for the current ray, thus deciding what is handled by calculations for current iteration
        for(Sphere sphere : closestVals.keySet()) {
            closestSphere = sphere;
            closest_St = closestVals.get(sphere);
            double[] oDSphereT = {oVector[0] + closestVals.get(sphere)*dVector[0], oVector[1] + closestVals.get(sphere)*dVector[1], oVector[2] + closestVals.get(sphere)*dVector[2]};
            sphereDistance = calcLength(oDSphereT);
        }
        for(WingedEdgeMesh m : closestMeshVals.keySet()) {
            if(m != null) {
                closestMesh = m;
                closest_Mt = closestMeshVals.get(m)[0];
                double[] oDMeshT = {oVector[0]+closestMeshVals.get(m)[0]*dVector[0], oVector[1]+closestMeshVals.get(m)[0]*dVector[1], oVector[2]+closestMeshVals.get(m)[0]*dVector[2]};
                meshDistance = calcLength(oDMeshT);
            }
        }
        if(sphereDistance <= meshDistance) {
            closest_t = closest_St;
            sphereIsClosest = true;
        } else if(sphereDistance > meshDistance && meshes.length > 0) {
            closest_t = closest_Mt;
            sphereIsClosest = false;
        } else {
            closest_t = closest_St;
            sphereIsClosest = true;
        }
        // ensures that a generic background color is returned if no object is hit
        if(sphereIsClosest) {
            if(closestSphere == null) {
                sphereColor[0] = 0;
                sphereColor[1] = 0;
                sphereColor[2] = 0;
                return sphereColor;
            }
        } else {
            if(closestMesh == null) {
                sphereColor[0] = 0;
                sphereColor[1] = 0;
                sphereColor[2] = 0;
                return sphereColor;
            }
        }

        // gets the color of the object hit
        if(sphereIsClosest) {
            sphereColor = closestSphere.color;
        } else {
            meshColor = closestMesh.getColor();
        }

        // defining the point that is hit by the light particle
        pVector[0] = oVector[0] + closest_t*dVector[0];
        pVector[1] = oVector[1] + closest_t*dVector[1];
        pVector[2] = oVector[2] + closest_t*dVector[2];
        
        // calculates normal vector for the sphere if sphere is closest, gets the returned normal for the mesh face if mesh is closest
        if(sphereIsClosest) {
            nVector[0] = pVector[0] - closestSphere.center[0];
            nVector[1] = pVector[1] - closestSphere.center[1];
            nVector[2] = pVector[2] - closestSphere.center[2];
            nLength = calcLength(nVector);
            nVector[0] = nVector[0]/nLength;
            nVector[1] = nVector[1]/nLength;
            nVector[2] = nVector[2]/nLength;
        } else {
            nVector[0] = closestMeshVals.get(closestMesh)[1];
            nVector[1] = closestMeshVals.get(closestMesh)[2];
            nVector[2] = closestMeshVals.get(closestMesh)[3];
        }

        if(sphereIsClosest) {
            iVal = computeLighting(lights, pVector, nVector, negativeDVector, closestSphere.specular, spheres, meshes);
        } else {
            iVal = computeLighting(lights, pVector, nVector, negativeDVector, closestMesh.getSpecular(), spheres, meshes);
        }
        
        if(sphereIsClosest) {
            localColor[0] = (short)(sphereColor[0]*iVal);
            localColor[1] = (short)(sphereColor[1]*iVal);
            localColor[2] = (short)(sphereColor[2]*iVal);
        } else {
            localColor[0] = (short)(meshColor[0]*iVal);
            localColor[1] = (short)(meshColor[1]*iVal);
            localColor[2] = (short)(meshColor[2]*iVal);
        }
        if(localColor[0] > 255) {
            localColor[0] = 255;
        }
        if(localColor[1] > 255) {
            localColor[1] = 255;
        }
        if(localColor[2] > 255) {
            localColor[2] = 255;
        }

        // recursive part - called for rays to bounce a # of times = recursion depth, instead of hitting just once
        if(sphereIsClosest) {
            rflctv = closestSphere.reflective;
        } else {
            rflctv = closestMesh.getReflective();
        }
        if(recursionDepth <= 0 || rflctv <= 0) {
            return localColor;
        }

        double[] rVector = reflectRay(negativeDVector, nVector);
        short[] reflectedColor = traceRay(pVector, rVector, 0.001, Double.MAX_VALUE, spheres, recursionDepth - 1, meshes);

        localColor[0] = (short)(localColor[0]*(1-rflctv) + reflectedColor[0]*rflctv);
        localColor[1] = (short)(localColor[1]*(1-rflctv) + reflectedColor[1]*rflctv);
        localColor[2] = (short)(localColor[2]*(1-rflctv) + reflectedColor[2]*rflctv);
        if(localColor[0] > 255) {
            localColor[0] = 255;
        }
        if(localColor[1] > 255) {
            localColor[1] = 255;
        }
        if(localColor[2] > 255) {
            localColor[2] = 255;
        }
        return localColor;
    }

    public double[] intersectSphere(double[] oVector, double[] dVector, Sphere sphere) {
        double[] tVals = {0.0, 0.0};
        double r = sphere.radius;
        double[] co = {oVector[0] - sphere.center[0], oVector[1] - sphere.center[1], oVector[2] - sphere.center[2]};

        double a = dot(dVector, dVector); // dot D D
        double b = 2.0 * dot(co, dVector); // 2* dot co D
        double c = dot(co, co) - r*r; // dot co co - r^2

        double root = b*b - 4*a*c;
        if(root < 0) {
            tVals[0] = Double.MAX_VALUE;
            tVals[1] = Double.MAX_VALUE;
            return tVals;
        }
        tVals[0] = (-1*b + Math.sqrt(root)) / (2*a);
        tVals[1] = (-1*b - Math.sqrt(root)) / (2*a);

        return tVals;
    }

    public double[] intersectMesh(double[] oVector, double[] dVector, WingedEdgeMesh mesh) {
        // oVector = [i, j, k] is start of ray, dVector = [l, m, n] is its direction
        // plane eqn: ax + by + cz + d = 0. a, b, and c are from normal to plane, plug in pt from face to solve for d
        // intersection t val eqn = (ai + bj + ck + d)/(al + bm + cn)
        double[] meshCenter = mesh.getCenter();

        // calculate the eqns of each of the faces normals, planes, determine if intersection between each
        double[] faceTVals = new double[mesh.getFaceList().length];
        double smallestT = Double.MAX_VALUE;
        for(int i = 0; i < faceTVals.length; i++) {
            faceTVals[i] = Double.MAX_VALUE;
        }
        double[] tVals = {0.0, 0.0, 0.0, 0.0};
        double[][] normalVectors = new double[mesh.getFaceList().length][3];

        double[] vert0 = new double[3];
        double[] vert1 = new double[3];
        double[] vert2 = new double[3];

        double[] edge01 = new double[3];
        double[] edge12 = new double[3];
        double[] edge20 = new double[3];
        double edge01Length = 0.0;
        double edge12Length = 0.0;
        double edge20Length = 0.0;

        double[] normalVector = new double[3];

        double[] faceCenter = new double[3];
        double faceDistance = 0.0;
        double[] faceCenterPlusNormal = new double[3];
        double faceCenterPlusNormalDistance = 0.0;

        double numerator = 0.0;
        double denominator = 0.0;

        double[] intersectionPoint = new double[3];
        double[] ap = new double[3];
        double[] bp = new double[3];
        double[] cp = new double[3];
        double apLength = 0.0;
        double bpLength = 0.0;
        double cpLength = 0.0;

        double currentT = 0.0;

        double e01AndApAngle = 0.0;
        double e12AndBpAngle = 0.0;
        double e20AndCpAngle = 0.0;

        double angleABC = 0.0;
        double areaABC = 0.0;
        double areaABP = 0.0;
        double areaBCP = 0.0;
        double areaCAP = 0.0;

        double uWeight = 0.0;
        double vWeight = 0.0;
        double wWeight = 0.0;


        for(int i = 0; i < mesh.getFaceList().length; i++) {
            // face normal is cross product of two edges in triangle
            vert0 = mesh.getVertexList()[mesh.getFaceList()[i][0]][0];
            vert1 = mesh.getVertexList()[mesh.getFaceList()[i][1]][0];
            vert2 = mesh.getVertexList()[mesh.getFaceList()[i][2]][0];

            edge01 = subtractVectors(vert1, vert0);
            edge12 = subtractVectors(vert2, vert1);
            edge20 = subtractVectors(vert0, vert2);

            normalVector = crossProduct3D(edge01, edge12);
            double normalVectorLength = calcLength(normalVector);
            normalVector[0] = normalVector[0] / normalVectorLength;
            normalVector[1] = normalVector[1] / normalVectorLength;
            normalVector[2] = normalVector[2] / normalVectorLength;

            // ensure normal direction is correct
            faceCenter = multVectorByScalar(addVectors(addVectors(vert0, vert1), vert2), 0.33);
            faceDistance = calcLength(subtractVectors(faceCenter, meshCenter));
            faceCenterPlusNormal = addVectors(faceCenter, multVectorByScalar(normalVector,faceDistance));
            faceCenterPlusNormalDistance = calcLength(subtractVectors(faceCenterPlusNormal, meshCenter));
            // if dist to the point of the face center plus the normal is smaller than the dist to the face center, normal is backwards
            if(faceCenterPlusNormalDistance < faceDistance) {
                normalVector = multVectorByScalar(normalVector, -1.0);
            }

            // add normal vector to matrix of normal vectors for each face
            normalVectors[i][0] = normalVector[0];
            normalVectors[i][1] = normalVector[1];
            normalVectors[i][2] = normalVector[2];

            // next, determine if there is an intersection with the plane the face is in
            double[] planeEqnVector = {normalVector[0], normalVector[1], normalVector[2], -1*(normalVector[0]*vert1[0] + normalVector[1]*vert1[1] + normalVector[2]*vert1[2])};
            numerator = -1*(planeEqnVector[0] * oVector[0] + planeEqnVector[1] * oVector[1] + planeEqnVector[2] * oVector[2] + planeEqnVector[3]);
            denominator = planeEqnVector[0] * dVector[0] + planeEqnVector[1] * dVector[1] + planeEqnVector[2] * dVector[2];

            // if denominator of eqn is 0, then vector is parallel, doesn't need to be calculated for
            if(denominator >= 0.001 || denominator <= -0.001) {
            // if(denominator != 0) {
                // now check if point of intersection is actually within the face
                currentT = numerator/denominator; // can safely divide now that we know denom != 0
                if(currentT > 0) {
                    intersectionPoint = addVectors(oVector, multVectorByScalar(dVector, currentT)); // oVector + currentT*dVector

                    // use convex hull test to see if the intersection point is within triagle
                    // need vectors from each point defining face to the intersection point
                    ap = subtractVectors(vert0, intersectionPoint);
                    bp = subtractVectors(vert1, intersectionPoint);
                    cp = subtractVectors(vert2, intersectionPoint);

                    edge01Length = calcLength(edge01);
                    edge12Length = calcLength(edge12);
                    edge20Length = calcLength(edge20);

                    apLength = calcLength(ap);
                    bpLength = calcLength(bp);
                    cpLength = calcLength(cp);

                    e01AndApAngle = Math.asin(calcLength(crossProduct3D(edge01, ap)) / (edge01Length * apLength));
                    e12AndBpAngle = Math.asin(calcLength(crossProduct3D(edge12, bp)) / (edge12Length * bpLength));
                    e20AndCpAngle = Math.asin(calcLength(crossProduct3D(edge20, cp)) / (edge20Length * cpLength));

                    // Barycentric weight test
                    angleABC = Math.asin(calcLength(crossProduct3D(edge01, edge12)) / (edge01Length * edge12Length));
                    areaABC = edge01Length * edge12Length * Math.sin(angleABC) / 2.0;
                    areaABP = edge01Length * apLength * Math.sin(e01AndApAngle) / 2.0;
                    areaBCP = edge12Length * bpLength * Math.sin(e12AndBpAngle) / 2.0;
                    areaCAP = edge20Length * cpLength * Math.sin(e20AndCpAngle) / 2.0;

                    uWeight = areaABP/areaABC;
                    vWeight = areaBCP/areaABC;
                    wWeight = areaCAP/areaABC;

                    // <= 1.001 to allow a tiny bit of wiggle room - removes noise present in image when check is for <= 1.0
                    if((0 <= uWeight+vWeight+wWeight && uWeight+vWeight+wWeight <= 1.001)) {
                        faceTVals[i] = currentT;
                        if(currentT < smallestT) {
                            smallestT = currentT;
                        }
                    }
                }
            }
        }

        for(int i = 0; i < faceTVals.length; i++) {
            if(faceTVals[i] == smallestT) {
                tVals[0] = smallestT;
                tVals[1] = normalVectors[i][0];
                tVals[2] = normalVectors[i][1];
                tVals[3] = normalVectors[i][2];
            }
        }
        return tVals;
    }

    // add two vectors in R3
    public double[] addVectors(double[] a, double[] b) {
        double[] c = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
        return c;
    }

    // vector b - a
    public double[] subtractVectors(double[] a, double[] b) {
        double[] c = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        return c;
    }

    // calculates dot product of two vectors in R3
    public double dot(double[] a, double[] b) {
        double val = 0.0;
        for(int i = 0; i < 3; i++) {
            val += a[i]*b[i];
        }
        return val;
    } 

    // calculates length of a vector
    public double calcLength(double[] a) {
        double l = Math.sqrt(Math.pow(a[0], 2) + Math.pow(a[1], 2) + Math.pow(a[2], 2));

        return l;
    }

    // linear interpolation method used for de casteljaus algorithm
    public double[] lerp(double[] p1, double[] p2, double t) {
        double[] lerpP = {(1-t)*p1[0] + t*p2[0], (1-t)*p1[1] + t*p2[1]};

        return lerpP;
    }

    // calculates light intensity for pixel being rendered by iterating through light sources, determining what is occluded, calculating how bright
    // reflection is based on angle light is hitting surface at - loose description
    public double computeLighting(Light[] lights, double[] p, double[] n, double[] v, double specular, Sphere[] spheres, WingedEdgeMesh[] meshes) {
        double i = 0.0;
        double nDotL = 0.0;
        double rDotV = 0.0;

        for(Light light : lights) {
            if(light.type.equals("ambient")) {
                i += light.intensity;
            } else {
                double[] l = {0.0, 0.0, 0.0};
                double[] r = {0.0, 0.0, 0.0};
                if(light.type.equals("point")) {
                    l[0] = light.position[0] - p[0];
                    l[1] = light.position[1] - p[1];
                    l[2] = light.position[2] - p[2];
                } else {
                    // may need to change above else back to else if the if statement below
                    // if(light.type.equals("directional")) {
                    l[0] = light.direction[0];
                    l[1] = light.direction[1];
                    l[2] = light.direction[2];
                }

                boolean occluded = false;
                HashMap<Sphere, Double> shadowedSpheres = closestIntersection(p, l, 0.001, Double.MAX_VALUE, spheres);
                HashMap<WingedEdgeMesh, double[]> shadowedMeshes = closestFVIntersection(p, l, 0.001, Double.MAX_VALUE, meshes);
                for(Sphere sphere : spheres) {
                    if(shadowedSpheres.containsKey(sphere)) {
                        occluded = true;
                    }
                }
                for(WingedEdgeMesh mesh : meshes) {
                    if(shadowedMeshes.containsKey(mesh)) {
                        occluded = true;
                    }
                }

                if(!occluded) {
                    // diffuse
                    nDotL = dot(n, l);
                    if(nDotL > 0) {
                        i += light.intensity * nDotL/(calcLength(n) * calcLength(l));
                    }

                    // specular
                    if(specular != -1) {
                        // if dot product of normal n and light direction l is 0, then the vectors are perpendicular and there is no reflection
                        r[0] = 2*n[0]*dot(n,l) - l[0];
                        r[1] = 2*n[1]*dot(n,l) - l[1];
                        r[2] = 2*n[2]*dot(n,l) - l[2];
                        rDotV = dot(r,v);
                        if(rDotV > 0) {
                            i += light.intensity*Math.pow(rDotV/(calcLength(r)*calcLength(v)), specular);
                        }
                    }
                }
            }
        }
        return i;
    }

    // de Casteljaus algorithm
    public void deCasteljaus(BezierCurve bc, JPanel curveCanvas) {
        double[] pointA = {0.0, 0.0};
        double[] pointB = {0.0, 0.0};
        double[] pointC = {0.0, 0.0};
        double[] pointD = {0.0, 0.0};
        double[] pointE = {0.0, 0.0};
        double[] pointP = {0.0, 0.0};

        for(double t = 0.0; t <= 1.0; t += 0.0001) {
            pointA = lerp(bc.point1, bc.point2, t);
            pointB = lerp(bc.point2, bc.point3, t);
            pointC = lerp(bc.point3, bc.point4, t);
            pointD = lerp(pointA, pointB, t);
            pointE = lerp(pointB, pointC, t);
            pointP = lerp(pointD, pointE, t);

            putPixel(curveCanvas, (int)pointP[0], (int)pointP[1], Color.BLACK);
        }
    }

    // paints a pixes to a JPanel canvas
    public void putPixel(JPanel curveCanvas, int x, int y, Color color) {
        curveCanvas.getGraphics().setColor(color);
        curveCanvas.getGraphics().drawLine(x, y, x, y);
    }

    // defines a bezier curve object defined by control points, then runs de casteljaus algorithm on it to get curve
    public void drawBezCurves(JPanel curveCanvas, int cpx1, int cpy1, int cpx2, int cpy2, int cpx3, int cpy3, int cpx4, int cpy4) {
        double[] bcPoint1 = {cpx1, cpy1};
        double[] bcPoint2 = {cpx2, cpy2};
        double[] bcPoint3 = {cpx3, cpy3};
        double[] bcPoint4 = {cpx4, cpy4};
        BezierCurve bc = new BezierCurve(bcPoint1, bcPoint2, bcPoint3, bcPoint4);
        deCasteljaus(bc, curveCanvas);
    }

    // determines a sphere's t values the d vector can be multiplied by - how far along D direction is traveled, basically.
    public HashMap<Sphere, Double> closestIntersection(double[] oVector, double[] dVector, double t_min, double t_max, Sphere[] spheres) {
        double closest_t = Double.MAX_VALUE;
        Sphere closestSphere = null;
        HashMap<Sphere, Double> sphere_t_pair = new HashMap<>();

        for(Sphere sphere : spheres) {
            double[] tVals = intersectSphere(oVector, dVector, sphere);
            if(tVals[0] > t_min && tVals[0] <= t_max && tVals[0] < closest_t) {
                closest_t = tVals[0];
                closestSphere = sphere;
            }
            if(tVals[1] > t_min && tVals[1] <= t_max && tVals[1] < closest_t) {
                closest_t = tVals[1];
                closestSphere = sphere;
            }
        }

        sphere_t_pair.put(closestSphere, closest_t);
        return sphere_t_pair;
    }

    // same as above method, but for meshes, also includes normal for simplicity of calculations
    public HashMap<WingedEdgeMesh, double[]> closestFVIntersection(double[] oVector, double[] dVector, double t_min, double t_max, WingedEdgeMesh[] meshes) {
        double closest_t = Double.MAX_VALUE;
        WingedEdgeMesh closestMesh = null;
        HashMap<WingedEdgeMesh, double[]> mesh_t_pair = new HashMap<>();
        double[] newTVals = {0.0, 0.0, 0.0, 0.0};
        double[] tVals = new double[4];

        for(WingedEdgeMesh mesh : meshes) {
            tVals = intersectMesh(oVector, dVector, mesh);
            if(tVals[0] > t_min && tVals[0] <= t_max && tVals[0] < closest_t) {
                closest_t = tVals[0];
                closestMesh = mesh;
                newTVals[0] = closest_t;
                newTVals[1] = tVals[1];
                newTVals[2] = tVals[2];
                newTVals[3] = tVals[3];
            }
        }

        mesh_t_pair.put(closestMesh, newTVals);
        return mesh_t_pair;
    }

    // rotates camera by input degrees along x, y, and z axis using rotation matrices
    public double[] rotateCamera(double[] dVector, double[] xyzAngle) {
        double[] newDVector = {0.0, 0.0, 0.0};

        double xCos = Math.cos(Math.toRadians(xyzAngle[0]));
        double xSin = Math.sin(Math.toRadians(xyzAngle[0]));
        double yCos = Math.cos(Math.toRadians(xyzAngle[1]));
        double ySin = Math.sin(Math.toRadians(xyzAngle[1]));
        double zCos = Math.cos(Math.toRadians(xyzAngle[2]));
        double zSin = Math.sin(Math.toRadians(xyzAngle[2]));

        double[][] rx = {{1.0, 0.0, 0.0}, {0.0, xCos, -1.0*xSin}, {0.0, xSin, xCos}}; // row 1, 2, 3
        double[][] ry = {{yCos, 0.0, ySin}, {0.0, 1.0, 0.0}, {-1.0*ySin, 0.0, yCos}};
        double[][] rz = {{zCos, -1.0*zSin, 0.0}, {zSin, zCos, 0.0}, {0.0, 0.0, 1.0}};
        
        // rotating by x angle
        if(xyzAngle[0] != 0.0) {
            newDVector[0] = rx[0][0]*dVector[0] + rx[0][1]*dVector[1] + rx[0][2]*dVector[2];
            newDVector[1] = rx[1][0]*dVector[0] + rx[1][1]*dVector[1] + rx[1][2]*dVector[2];
            newDVector[2] = rx[2][0]*dVector[0] + rx[2][1]*dVector[1] + rx[2][2]*dVector[2];
        }

        // rotating by y angle
        if(xyzAngle[1] != 0.0) {
            newDVector[0] = ry[0][0]*dVector[0] + ry[0][1]*dVector[1] + ry[0][2]*dVector[2];
            newDVector[1] = ry[1][0]*dVector[0] + ry[1][1]*dVector[1] + ry[1][2]*dVector[2];
            newDVector[2] = ry[2][0]*dVector[0] + ry[2][1]*dVector[1] + ry[2][2]*dVector[2];
        }

        // rotating by z angle
        if(xyzAngle[2] != 0.0) {
            newDVector[0] = rz[0][0]*dVector[0] + rz[0][1]*dVector[1] + rz[0][2]*dVector[2];
            newDVector[1] = rz[1][0]*dVector[0] + rz[1][1]*dVector[1] + rz[1][2]*dVector[2];
            newDVector[2] = rz[2][0]*dVector[0] + rz[2][1]*dVector[1] + rz[2][2]*dVector[2];
        }

        if(xyzAngle[0] == 0.0 && xyzAngle[1] == 0.0 && xyzAngle[2] == 0.0) {
            newDVector[0] = dVector[0];
            newDVector[1] = dVector[1];
            newDVector[2] = dVector[2];
        }

        return newDVector;
    }

    // calculates cross product of two vectors in R3
    public double[] crossProduct3D(double[] a, double[] b) {
        double[] c = {(a[1]*b[2] - a[2]*b[1]), (a[2]*b[0] - a[0]*b[2]), (a[0]*b[1]- a[1]*b[0])};
        return c;
    }

    // Loop's Subdivision Surface algorithm for triangle meshes
    public WingedEdgeMesh loopSubdivision(WingedEdgeMesh mesh, int recursionDepth) {
        // number of new edges = 2*numEdges + 3*numbFaces
        // number of new faces = 4*numbFaces
        // number of new points = numbPoints + numbEdges
        double[][][] newVL = new double[mesh.getEdges().length + mesh.getVertexList().length][][];
        // int[][] newFL = new int[4*mesh.getFaceList().length][];
        ArrayList<int[]> tempNewFL = new ArrayList<>();
        Edge[] meshEdges = mesh.getEdges();
        double[][][] onlyNewVerts = new double[meshEdges.length][][]; // these HAVE to be the first |meshEdges| of newVL
        HashMap<Edge, double[][]> edgesAndCoords = new HashMap<>();
        HashMap<Edge, Integer> edgesAndIndices = new HashMap<>();


        for(int i = 0; i < meshEdges.length; i++) {
            double[] oldEdgePoint1 = mesh.getVertexList()[meshEdges[i].getVertexIndexList()[0]][0];
            double[] oldEdgePoint2 = mesh.getVertexList()[meshEdges[i].getVertexIndexList()[1]][0];
            // System.out.println(oldEdgePoint1[0] + ", " + oldEdgePoint1[1] + ", " + oldEdgePoint1[2]);
            int face0 = meshEdges[i].getFacePair().get(0);
            int face1 = 0;
            if(meshEdges[i].getFacePair().size() == 1) {
                face1 = face0;
            } else {
                face1 = meshEdges[i].getFacePair().get(1);
            }
            double[][] newVertex = {addVectors(multVectorByScalar(oldEdgePoint1, 0.5), multVectorByScalar(oldEdgePoint2, 0.5)), {face0, face1}}; 
            System.out.println("New Vertex " + i + ": " + newVertex[0][0] + ", " + newVertex[0][1] + ", " + newVertex[0][2]);
            
            newVL[i] = newVertex;
            edgesAndCoords.put(meshEdges[i], newVertex);
            edgesAndIndices.put(meshEdges[i], i);
        }
        for(int i = 0; i < mesh.getVertexList().length; i++) {
            newVL[i+meshEdges.length] = mesh.getVertexList()[i];
        }

        // regular division of each edge into 2 working correctly
        // for(int i = 0; i < newVL.length; i++) {
        //     System.out.println(i + ": " + newVL[i][0][0] + ", " + newVL[i][0][1] + ", " + newVL[i][0][2]);
        // }

        // go through all the faces of the old triangle. Inner-loop until all three new points are found. create four new faces.
        int faceSection = 0;
        // System.out.println(mesh.getFaceList().length);
        for(int i = 0; i < mesh.getFaceList().length; i++) {
            int facePoint1 = mesh.getFaceList()[i][0];
            int facePoint2 = mesh.getFaceList()[i][1];
            int facePoint3 = mesh.getFaceList()[i][2];

            // System.out.println("Face: " + i + ", Vertices: " + facePoint1 + ", " + facePoint2 + ", " + facePoint3);

            int[] newFace1 = new int[3];
            int[] newFace2 = new int[3];
            int[] newFace3 = new int[3];
            
            // cycle through the three points on the face. for each, determine the two edges connected. use those edges to get the pts needed to make a new face
            Edge e1 = mesh.getEdges()[mesh.getFaceEdgeList()[i][0]];
            Edge e2 = mesh.getEdges()[mesh.getFaceEdgeList()[i][1]];
            Edge e3 = mesh.getEdges()[mesh.getFaceEdgeList()[i][2]];

            // System.out.println("Face Edges: " + e1.getVertexIndexList()[0] + "-" + e1.getVertexIndexList()[1] + ", " + e2.getVertexIndexList()[0] + "-" + e2.getVertexIndexList()[1] + 
            // ", " + e3.getVertexIndexList()[0] + "-" + e3.getVertexIndexList()[1]);

            // point 1
            if(facePoint1 != e1.getVertexIndexList()[0] && facePoint1 != e1.getVertexIndexList()[1]) {
                newFace1[0] = edgesAndIndices.get(e2);
                newFace1[1] = edgesAndIndices.get(e3);
                newFace1[2] = facePoint1+meshEdges.length;
            } else if(facePoint1 != e2.getVertexIndexList()[0] && facePoint1 != e2.getVertexIndexList()[1]) {
                newFace1[0] = edgesAndIndices.get(e1);
                newFace1[1] = edgesAndIndices.get(e3);
                newFace1[2] = facePoint1+meshEdges.length;
            } else if(facePoint1 != e3.getVertexIndexList()[0] || facePoint1 != e3.getVertexIndexList()[1]) {
                newFace1[0] = edgesAndIndices.get(e1);
                newFace1[1] = edgesAndIndices.get(e2);
                newFace1[2] = facePoint1+meshEdges.length;
            }

            // point 2
            if(facePoint2 != e1.getVertexIndexList()[0] && facePoint2 != e1.getVertexIndexList()[1]) {
                newFace2[0] = edgesAndIndices.get(e2);
                newFace2[1] = edgesAndIndices.get(e3);
                newFace2[2] = facePoint2+meshEdges.length;
            } else if(facePoint2 != e2.getVertexIndexList()[0] && facePoint2 != e2.getVertexIndexList()[1]) {
                newFace2[0] = edgesAndIndices.get(e1);
                newFace2[1] = edgesAndIndices.get(e3);
                newFace2[2] = facePoint2+meshEdges.length;
            } else if(facePoint2 != e3.getVertexIndexList()[0] || facePoint2 != e3.getVertexIndexList()[1]) {
                newFace2[0] = edgesAndIndices.get(e1);
                newFace2[1] = edgesAndIndices.get(e2);
                newFace2[2] = facePoint2+meshEdges.length;
            }

            // point 3
            if(facePoint3 != e3.getVertexIndexList()[0] && facePoint3 != e1.getVertexIndexList()[1]) {
                newFace3[0] = edgesAndIndices.get(e2);
                newFace3[1] = edgesAndIndices.get(e3);
                newFace3[2] = facePoint3+meshEdges.length;
            } else if(facePoint3 != e2.getVertexIndexList()[0] && facePoint3 != e2.getVertexIndexList()[1]) {
                newFace3[0] = edgesAndIndices.get(e1);
                newFace3[1] = edgesAndIndices.get(e3);
                newFace3[2] = facePoint3+meshEdges.length;
            } else if(facePoint3!= e3.getVertexIndexList()[0] || facePoint3 != e3.getVertexIndexList()[1]) {
                newFace3[0] = edgesAndIndices.get(e1);
                newFace3[1] = edgesAndIndices.get(e2);
                newFace3[2] = facePoint3+meshEdges.length;
            }

            int[] newFace4 = {edgesAndIndices.get(e1), edgesAndIndices.get(e2), edgesAndIndices.get(e3)};

            // System.out.println("New Face 1: \n");
            // for(int k = 0; k < newFace1.length; k++) {
            //     System.out.print(newVL[newFace1[k]][0][0] + ", " + newVL[newFace1[k]][0][1] + ", " + newVL[newFace1[k]][0][2] + "\n");
            // }
            // System.out.println("\n");

            // System.out.println("New Face 2: \n");
            // for(int k = 0; k < newFace1.length; k++) {
            //     System.out.print(newVL[newFace2[k]][0][0] + ", " + newVL[newFace2[k]][0][1] + ", " + newVL[newFace2[k]][0][2] + "\n");
            // }
            // System.out.println("\n");

            // System.out.println("New Face 3: \n");
            // for(int k = 0; k < newFace1.length; k++) {
            //     System.out.print(newVL[newFace3[k]][0][0] + ", " + newVL[newFace3[k]][0][1] + ", " + newVL[newFace3[k]][0][2] + "\n");
            // }
            // System.out.println("\n");

            // System.out.println("New Face 4: \n");
            // for(int k = 0; k < newFace1.length; k++) {
            //     System.out.print(newVL[newFace4[k]][0][0] + ", " + newVL[newFace4[k]][0][1] + ", " + newVL[newFace4[k]][0][2] + "\n");
            // }
            // System.out.println("\n");


            tempNewFL.add(newFace1);
            tempNewFL.add(newFace2);
            tempNewFL.add(newFace3);
            tempNewFL.add(newFace4);
            faceSection += 4;
        }

        int[][] newFL = new int[tempNewFL.size()][];

        for(int i = 0; i < tempNewFL.size(); i++) {
            newFL[i] = tempNewFL.get(i);
            // System.out.println(i + ": " + newFL[i][0] + ", " + newFL[i][1] + ", " + newFL[i][2]); faces also working correctly
        }

        // updates the faces on the vertexlist matrix
        for(int i = 0; i < newVL.length; i++) {
            int numTimes = 0;
            for(int k = 0; k < newFL.length; k++) {
                if(newFL[k][0] == i) {
                    numTimes++;
                } else if(newFL[k][1] == i) {
                    numTimes++;
                } else if(newFL[k][2] == i) {
                    numTimes++;
                }
            }
            double[] vertsOnFace = new double[numTimes];
            int faceNumb = 0;
            for(int k = 0; k < newFL.length; k++) {
                for(int j = 0; j < newFL[k].length; j++) {
                    if(newFL[k][j] == i) {
                        vertsOnFace[faceNumb] = k;
                        faceNumb++;
                    }
                }
            }
            newVL[i][1] = vertsOnFace;
        }

        System.out.println("Debugging newVL faces for each vertex");
        for(int i = 0; i < newFL.length; i++) {
            System.out.println("Face " + i + ": ");
            for(int j = 0; j < newFL[i].length; j++) {
                System.out.println(j);
            }
        }
        for(int i = 0; i < newVL.length; i++) {
            System.out.println(i + ": " + newVL[i][0][0] + ", " + newVL[i][0][1] + ", " + newVL[i][0][2]);
            for(int j = 0; j < newVL[i][1].length; j++) {
                System.out.println("Face: " + newVL[i][1][j]);
            }
        }
        System.out.println();

        for(int i = 0; i < mesh.getVertexEdgeList().length; i++) {
            // determine the number of connected edges n
            double[] centerVertex = mesh.getVertexList()[i][0];
            // double n = 0.0;
            ArrayList<Integer> nCounter = new ArrayList<>();
            double n = mesh.getVertexEdgeList()[i].length;
            double beta = 0.0;
            if(n==3) {
                // beta = 3.0/16.0;
                beta = 3/16.0;
            } else {
                beta = 3.0/(8.0*n);
            }
            // double otherWeight = 1.0 - n*beta;
            double otherWeight = 1.0 - n*beta;
            double xUnitSum = 0.0;
            double yUnitSum = 0.0;
            double zUnitSum = 0.0;
        
            if(n >= 3) {
                // identify boundary vertices
                boolean isBounaryVert = false;
                for(int k = 0; k < n; k++) {
                    Edge connectedEdge = mesh.getEdges()[mesh.getVertexEdgeList()[i][k]];
                    if(connectedEdge.getIsBoundaryEdge()) {
                        isBounaryVert = true;
                    }
                }

                for(int k = 0; k < n; k++) {
                    Edge connectedEdge = mesh.getEdges()[mesh.getVertexEdgeList()[i][k]];
                    if(isBounaryVert) {
                        if(connectedEdge.getIsBoundaryEdge()) {
                            xUnitSum += (1.0/8.0)*edgesAndCoords.get(connectedEdge)[0][0];
                        }
                    } else {
                        xUnitSum += beta*edgesAndCoords.get(connectedEdge)[0][0];
                    }
                }
                if(isBounaryVert) {
                    xUnitSum += 0.75*centerVertex[0];
                } else {
                    xUnitSum += otherWeight*centerVertex[0];
                }

                for(int k = 0; k < n; k++) {
                    Edge connectedEdge = mesh.getEdges()[mesh.getVertexEdgeList()[i][k]];
                    if(isBounaryVert) {
                        if(connectedEdge.getIsBoundaryEdge()) {
                            yUnitSum += (1.0/8.0)*edgesAndCoords.get(connectedEdge)[0][1];
                        }
                    } else {
                        yUnitSum += beta*edgesAndCoords.get(connectedEdge)[0][1];
                    }
                }
                if(isBounaryVert) {
                    yUnitSum += 0.75*centerVertex[1];
                } else {
                    yUnitSum += otherWeight*centerVertex[1];
                }

                for(int k = 0; k < n; k++) {
                    Edge connectedEdge = mesh.getEdges()[mesh.getVertexEdgeList()[i][k]];
                    if(isBounaryVert) {
                        if(connectedEdge.getIsBoundaryEdge()) {
                            zUnitSum += (1.0/8.0)*edgesAndCoords.get(connectedEdge)[0][2];
                        }
                    } else {
                        zUnitSum += beta*edgesAndCoords.get(connectedEdge)[0][2];
                    }
                }
                if(isBounaryVert) {
                    zUnitSum += 0.75*centerVertex[2];
                } else {
                    zUnitSum += otherWeight*centerVertex[2];
                }
            } else if(n == 2) {
                // System.out.println("i: " + i);
                for(int k = 0; k < n; k++) {
                    Edge connectedEdge = mesh.getEdges()[mesh.getVertexEdgeList()[i][k]];
                    // System.out.println(edgesAndCoords.get(connectedEdge)[0][0]);
                    xUnitSum += (1.0/8.0)*edgesAndCoords.get(connectedEdge)[0][0];
                }
                // System.out.println("xUnitSum: " + xUnitSum + ", centerVertex[0]: " + centerVertex[0]);
                xUnitSum += 0.75*centerVertex[0];
                // System.out.println();
                for(int k = 0; k < n; k++) {
                    Edge connectedEdge = mesh.getEdges()[mesh.getVertexEdgeList()[i][k]];
                    // System.out.println(edgesAndCoords.get(connectedEdge)[0][1]);
                    yUnitSum += (1.0/8.0)*edgesAndCoords.get(connectedEdge)[0][1];
                }
                // System.out.println("yUnitSum: " + yUnitSum + ", centerVertex[1]: " + centerVertex[1]);
                yUnitSum += 0.75*centerVertex[1];
                // System.out.println();
                for(int k = 0; k < n; k++) {
                    Edge connectedEdge = mesh.getEdges()[mesh.getVertexEdgeList()[i][k]];
                    zUnitSum += (1.0/8.0)*edgesAndCoords.get(connectedEdge)[0][2];
                }
                zUnitSum += 0.75*centerVertex[2];
                // System.out.println();
            }
            
            double[] newCenterVertex = {xUnitSum, yUnitSum, zUnitSum};
            // System.out.println(i + ": " + newCenterVertex[0] + ", " + newCenterVertex[1] + ", " + newCenterVertex[2] + ", n = " + n);
            mesh.getVertexList()[i][0] = newCenterVertex;
        }

        WingedEdgeMesh newMesh = new WingedEdgeMesh(newFL, newVL, mesh.getColor(), mesh.getSpecular(), mesh.getReflective(), mesh.getCenter());

        return newMesh;
    }
}

class PixelPainter extends Canvas {
    public void putPixel(Graphics g, int x, int y, Color color) {
        g.setColor(color);
        g.drawLine(x, y, x, y);
    }
}
