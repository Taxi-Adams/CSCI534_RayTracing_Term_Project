import java.awt.Color;

public class RTContainer implements Runnable{
    Camera camera;
    Sphere[] spheres;
    int recursionDepth;
    PixelPainter canvas;
    double[] xyzAngle;
    int numThreads;
    int threadNumber;
    RayTracer rayTracer;
    int point1;
    int point2;
    WingedEdgeMesh[] meshes;


    public RTContainer(Camera camera, Sphere[] spheres, int recursionDepth, PixelPainter canvas, double[] xyzAngle, int numThreads, int threadNumber, RayTracer rayTracer, 
    int point1, int point2, WingedEdgeMesh[] meshes) {
        this.camera = camera;
        this.spheres = spheres;
        this.recursionDepth = recursionDepth;
        this.canvas = canvas;
        this.xyzAngle = xyzAngle;
        this.numThreads = numThreads;
        this.threadNumber = threadNumber;
        this.rayTracer = rayTracer;
        this.point1 = point1;
        this.point2 = point2;
        this.meshes = meshes;
    }

    public void run() {
        // rayTracer.traceSceneInThread(camera, spheres, recursionDepth, canvas, xyzAngle, numThreads, threadNumber);
        rayTracer.traceScene(camera, spheres, recursionDepth, canvas, xyzAngle, point1, point2, meshes);
    }
}
