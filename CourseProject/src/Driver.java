import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.event.MouseInputListener;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;

public class Driver {

    private static int cpx1 = 0;
    private static int cpy1 = 0;
    private static int cpx2 = 0;
    private static int cpy2 = 0;
    private static int cpx3 = 0;
    private static int cpy3 = 0;
    private static int cpx4 = 0;
    private static int cpy4 = 0;

    private static boolean settingP1s = false;
    private static boolean settingP2s = false;
    private static boolean settingP3s = false;
    private static boolean settingP4s = false;

    private static int p1Count = 0;
    private static int p2Count = 0;
    private static int p3Count = 0;
    private static int p4Count = 0;

    private static int xDegOffset = 0;
    private static int yDegOffset = 0;
    // private static int zDegOffset = 0;

    public static void main(String[] args) throws Exception {
        RayTracer rayTracer = new RayTracer();
        // raytracer/subdivision surfaces pane
        JFrame frame = new JFrame("Ray Tracer"); // create a frame to hold the canvas
        Dimension rtDim = new Dimension(880, 595);
        frame.setMinimumSize(rtDim);
        frame.setLayout(null);
        PixelPainter canvas = new PixelPainter();
        //Canvas canvas = new Canvas(); // create a canvas to display images on
        canvas.setSize(880, 495);
        frame.add(canvas);
        Color backgroundColor = new Color(255, 255, 255);
        canvas.setBackground(backgroundColor);
        // double[] cameraPosition = {0.0, 0.0, 0.0};
        double[] cameraPosition = {1.7, -0.18, 0.0};
        Camera camera = new Camera(cameraPosition, 1.0);

        JPanel rtButtons = new JPanel();
        FlowLayout rtFlowLyt = new FlowLayout(5, 10, 10);
        rtButtons.setLayout(rtFlowLyt);
        rtButtons.setLocation(0, 496);
        rtButtons.setSize(700, 100);

        JButton posXRotate = new JButton("x Deg++");
        JButton negXRotate = new JButton("x Deg--");
        JButton posYRotate = new JButton("y Deg++");
        JButton negYRotate = new JButton("y Deg--");
        JButton subdivideBtn = new JButton("Subdivide");
        JButton resetAngle = new JButton("Reset Camera");

        rtButtons.add(posXRotate);
        rtButtons.add(negXRotate);
        rtButtons.add(posYRotate);
        rtButtons.add(negYRotate);
        rtButtons.add(subdivideBtn);
        rtButtons.add(resetAngle);

        frame.add(rtButtons);
        frame.pack();
        frame.setVisible(true);


        JFrame curveFrame = new JFrame("Bezier Curves");
        Dimension curveSize = new Dimension(400, 515);
        curveFrame.setLayout(null);
        curveFrame.setMinimumSize(curveSize);
        JPanel curvePanel = new JPanel();
        curvePanel.setLayout(null);
        JPanel buttonPanel = new JPanel();
        FlowLayout btnsManager = new FlowLayout();
        buttonPanel.setLayout(btnsManager);

        JButton drawThing = new JButton("Draw");
        JButton clearPanel = new JButton("Clear");

        JButton setP1 = new JButton("Set Cntrl Pt 1");
        JButton setP2 = new JButton("Set Cntrl Pt 2");
        JButton setP3 = new JButton("Set Cntrl Pt 3");
        JButton setP4 = new JButton("Set Cntrl Pt 4");

        MouseInputListener pointsMouse = new MouseInputListener() {
            @Override
            public void mousePressed(MouseEvent e) {
            }

            @Override
            public void mouseExited(MouseEvent e) {
            }

            @Override
            public void mouseDragged(MouseEvent e) {
            }

            @Override
            public void mouseReleased(MouseEvent e) {
            }

            @Override
            public void mouseClicked(MouseEvent e) {
                if(e.getX() >= 0 && e.getX() <= curvePanel.getWidth() && e.getY() >=0 && e.getY() <= curvePanel.getHeight()) {
                    if(settingP1s) {
                        if(p1Count > 0) {
                            p1Count = 0;
                            curvePanel.getGraphics().clearRect(cpx3, cpy3, 9, 9);
                        }
                        cpx1 = e.getX();
                        cpy1 = e.getY();
                        p1Count++;
                        System.out.println("Setting Control Point 1s. x: " + cpx1 + ", y: " + cpy1);
                        curvePanel.getGraphics().drawOval(cpx1, cpy1, 7, 7);
                        settingP1s = false;
                    } else if(settingP2s) {
                        if(p2Count > 0) {
                            p2Count = 0;
                            curvePanel.getGraphics().clearRect(cpx3, cpy3, 9, 9);
                        }
                        cpx2 = e.getX();
                        cpy2 = e.getY();
                        p2Count++;
                        System.out.println("Setting Control Point 1s. x: " + cpx2 + ", y: " + cpy2);
                        curvePanel.getGraphics().drawOval(cpx2, cpy2, 7, 7);
                        settingP2s = false;
                    } else if(settingP3s) {
                        if(p3Count > 0) {
                            p3Count = 0;
                            curvePanel.getGraphics().clearRect(cpx3, cpy3, 9, 9);
                        }
                        cpx3 = e.getX();
                        cpy3 = e.getY();
                        p3Count++;
                        System.out.println("Setting Control Point 3s. x: " + cpx3 + ", y: " + cpy3);
                        curvePanel.getGraphics().drawOval(cpx3, cpy3, 7, 7);
                        settingP3s = false;
                    } else if(settingP4s) {
                        if(p4Count > 0) {
                            p4Count = 0;
                            curvePanel.getGraphics().clearRect(cpx3, cpy3, 9, 9);
                        }
                        cpx4 = e.getX();
                        cpy4 = e.getY();
                        p4Count++;
                        System.out.println("Setting Control Point 4s. x: " + cpx4 + ", y: " + cpy4);
                        curvePanel.getGraphics().drawOval(cpx4, cpy4, 7,7);
                        settingP4s = false;
                    }
                }
            }

            @Override
            public void mouseMoved(MouseEvent e) {
            }

            @Override
            public void mouseEntered(MouseEvent e) {
            }
        };
        curvePanel.addMouseListener(pointsMouse);
        
        drawThing.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                rayTracer.drawBezCurves(curvePanel, cpx1, cpy1, cpx2, cpy2, cpx3, cpy3, cpx4, cpy4);
            }
        });

        clearPanel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                curvePanel.getGraphics().clearRect(0, 0, curvePanel.getWidth(), curvePanel.getHeight());
            }
        });

        setP1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                settingP1s = true;
                settingP2s = false;
                settingP3s = false;
                settingP4s = false;
            }
        });

        setP2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                settingP1s = false;
                settingP2s = true;
                settingP3s = false;
                settingP4s = false;
            }
        });

        setP3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                settingP1s = false;
                settingP2s = false;
                settingP3s = true;
                settingP4s = false;
            }
        });

        setP4.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                settingP1s = false;
                settingP2s = false;
                settingP3s = false;
                settingP4s = true;
            }
        });


        curvePanel.setSize(400, 400);
        curvePanel.setBackground(Color.WHITE);

        buttonPanel.setSize(400, 115);
        buttonPanel.setLocation(0, 400);
        buttonPanel.setBackground(Color.WHITE);
        buttonPanel.add(setP1);
        buttonPanel.add(setP2);
        buttonPanel.add(setP3);
        buttonPanel.add(setP4);
        buttonPanel.add(drawThing);
        buttonPanel.add(clearPanel);

        curveFrame.setLocation(900, 200);
        curveFrame.add(curvePanel);
        curveFrame.add(buttonPanel);
        curveFrame.pack();
        curveFrame.setVisible(true);

        int recursionDepth = 3; // change back to 3

        double[] center1 = {0.0, -1.0, 3.0};
        double[] center2 = {2.0, 0.0, 4.0};
        double[] center3 = {-2.0, 0.0, 4.0};
        double[] center4 = {0, -5001, 0};

        short[] color1 = {255, 0, 0};
        short[] color2 = {0, 0, 255};
        short[] color3 = {0, 255, 0};
        short[] color4 = {255, 255, 0};

        double specular1 = 500.0;
        double specular2 = 500.0;
        double specular3 = 10.0;
        double specular4 = 1000.0;

        double reflective1 = 0.2;
        double reflective2 = 0.3;
        double reflective3 = 0.4;
        double reflective4 = 0.5;

        Sphere sphere1 = new Sphere(color1, center1, 1.0, specular1, reflective1);
        // Sphere sphere2 = new Sphere(color2, center2, 1.0, specular2, reflective2); // test sphere
        // Sphere sphere3 = new Sphere(color3, center3, 1.0, specular3, reflective3); // test sphere
        Sphere sphere4 = new Sphere(color4, center4, 5000.0, specular4, reflective4);
        // Sphere[] spheres = {sphere1, sphere2, sphere3, sphere4};
        Sphere[] spheres = {sphere1, sphere4};

        int[][] cube1FL = {{0,1,2}, {1,2,3}, {1,3,5}, {3,5,7}, {4,5,6}, {5,6,7}, {0,4,6}, {0,2,6}, {0,1,8}, {1,5,8}, {4,5,8}, {0,4,8}, {2,3,9}, {3,7,9}, {6,7,9}, {2,6,9}};

        double[][] c1Vert1 = {{2,1.0,4}, {0,6,7,8,11}};
        double[][] c1Vert2 = {{3,1.0,4}, {0,1,2,8,9}};
        double[][] c1Vert3 = {{2,0.0,4}, {0,1,7,12,15}};
        double[][] c1Vert4 = {{3,0.0,4}, {1,2,3,12,13}};
        double[][] c1Vert5 = {{2,1,5}, {4,6,10,11}};
        double[][] c1Vert6 = {{3,1,5}, {2,3,4,5,9,10}};
        double[][] c1Vert7 = {{2,0,5}, {4,5,6,7,14,15}};
        double[][] c1Vert8 = {{3,0,5}, {3,5,13,14}};
        double[][] c1Vert9 = {{2.5,1.0,4.5}, {12,13,14,15}};
        double[][] c1Vert10 = {{2.5,0.0,4.5}, {8,9,10,11}};

        double[] c1Center = {2.5,0.5,4.5};

        double[][][] cube1VL = {c1Vert1, c1Vert2, c1Vert3, c1Vert4, c1Vert5, c1Vert6, c1Vert7, c1Vert8, c1Vert9, c1Vert10};
        WingedEdgeMesh cube1 = new WingedEdgeMesh(cube1FL, cube1VL, color2, specular2, reflective3, c1Center);
        cube1.setVertexEdgeList();
        cube1.setBoundaryEdges();

        // WingedEdgeMesh subDivCube1 = rayTracer.loopSubdivision(cube1, recursionDepth);

        int[][] triExFL = {{0,2,4},{1,2,3},{1,4,5},{1,2,4}};
        double[][][] triExVL = {{{2.0,0.95,4.0},{0}}, {{2.0,0.05,4.0},{1,2,3}}, {{1.6,0.66,4.0},{0,1,3}}, {{1.6,0.33,4.0},{1}}, {{2.4,0.66,4.0},{0,2,3}}, {{2.4,0.33,4.0},{2}}};
        WingedEdgeMesh triEx = new WingedEdgeMesh(triExFL, triExVL, color3, specular4, reflective4, c1Center);
        triEx.setVertexEdgeList();
        triEx.setBoundaryEdges();

        // for(int j = 0; j < triEx.getEdges().length; j++) {
        //     System.out.println("Edge:" + triEx.getEdges()[j].getVertexIndexList()[0] + "-" + triEx.getEdges()[j].getVertexIndexList()[1]);
        //     if(triEx.getEdges()[j].getIsBoundaryEdge()) {
        //         System.out.println("is a boundary");
        //     }
        // }

        // WingedEdgeMesh triExSub = rayTracer.loopSubdivision(triEx, recursionDepth);
        // triExSub.setVertexEdgeList();
        // triExSub.setBoundaryEdges();
        // WingedEdgeMesh triExSub2 = rayTracer.loopSubdivision(triExSub, recursionDepth);
        // triExSub2.setVertexEdgeList();
        // triExSub2.setBoundaryEdges();
        // WingedEdgeMesh triExSub3 = rayTracer.loopSubdivision(triExSub2, recursionDepth);

        WingedEdgeMesh[] meshes = {cube1};

        Light ambientLight1 = new Light("ambient", 0.2, null, null);
        double[] lightPosition1 = {2.0, 1.0, 0.0};
        Light pointLight1 = new Light("point", 0.6, lightPosition1, null);
        double[] lightDirection1 = {1.0, 4.0, 4.0};
        Light directionalLight1 = new Light("directional", 0.2, null, lightDirection1);
        rayTracer.getLights()[0] = ambientLight1;
        rayTracer.getLights()[1] = pointLight1;
        rayTracer.getLights()[2] = directionalLight1;

        // double[] oVector = {0.0, 0.0, 0.0}; // O vector
        double[] xyzAngle = {0.0, 0.0, 0.0};
        int rtCWLower = -1 * rayTracer.getCW()/2 - 10;
        int rtCWHigher = rayTracer.getCW()/2;
        // older threading code. sadly not worth the effort in the end, but included for proof of work, more than anything

        // int numThreads = 5;
        // int[] xRangeVals = {-1*rayTracer.getCW()/2 + rayTracer.getCW()/5, -1*rayTracer.getCW()/2 + 2*rayTracer.getCW()/5, -1*rayTracer.getCW()/2 + 3*rayTracer.getCW()/5, 
        //     -1*rayTracer.getCW()/2 + 4*rayTracer.getCW()/5, -1*rayTracer.getCW()/2 + 5*rayTracer.getCW()/5};
        // Runnable thrd1 = new RTContainer(camera, spheres, recursionDepth, canvas, xyzAngle, numThreads, 1, rayTracer, -1*rayTracer.getCW()/2 - 10, xRangeVals[0]);
        // Runnable thrd2 = new RTContainer(camera, spheres, recursionDepth, canvas, xyzAngle, numThreads, 2, rayTracer, xRangeVals[0], xRangeVals[1]);
        // Runnable thrd3 = new RTContainer(camera, spheres, recursionDepth, canvas, xyzAngle, numThreads, 3, rayTracer, xRangeVals[1], xRangeVals[2]);
        // Runnable thrd4 = new RTContainer(camera, spheres, recursionDepth, canvas, xyzAngle, numThreads, 4, rayTracer, xRangeVals[2], xRangeVals[3]);
        // Runnable thrd5 = new RTContainer(camera, spheres, recursionDepth, canvas, xyzAngle, numThreads, 5, rayTracer, xRangeVals[3], xRangeVals[4]);

        

        posXRotate.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                xDegOffset++;
                xyzAngle[0] = xDegOffset;
                // new Thread(thrd1).start();
                // new Thread(thrd2).start();
                // new Thread(thrd3).start();
                // new Thread(thrd4).start();
                // new Thread(thrd5).start();
                // xyzAngle[0] = +1;
                rayTracer.traceScene(camera, spheres, recursionDepth, canvas, xyzAngle, rtCWLower, rtCWHigher, meshes);
            }
        });

        negXRotate.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                xDegOffset--;
                xyzAngle[0] = xDegOffset;
                // new Thread(thrd1).start();
                // new Thread(thrd2).start();
                // new Thread(thrd3).start();
                // new Thread(thrd4).start();
                // new Thread(thrd5).start();
                // xyzAngle[0] = -1;
                rayTracer.traceScene(camera, spheres, recursionDepth, canvas, xyzAngle, rtCWLower, rtCWHigher, meshes);
            }
        });

        posYRotate.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                yDegOffset++;
                xyzAngle[1] = yDegOffset;
                // new Thread(thrd1).start();
                // new Thread(thrd2).start();
                // new Thread(thrd3).start();
                // new Thread(thrd4).start();
                // new Thread(thrd5).start();
                // xyzAngle[1] = +1;
                rayTracer.traceScene(camera, spheres, recursionDepth, canvas, xyzAngle, rtCWLower, rtCWHigher, meshes);
            }
        });

        negYRotate.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                yDegOffset--;
                xyzAngle[1] = yDegOffset;
                // new Thread(thrd1).start();
                // new Thread(thrd2).start();
                // new Thread(thrd3).start();
                // new Thread(thrd4).start();
                // new Thread(thrd5).start();
                // xyzAngle[1] = -1;
                rayTracer.traceScene(camera, spheres, recursionDepth, canvas, xyzAngle, rtCWLower, rtCWHigher, meshes);
            }
        });

        resetAngle.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent a) {
                xDegOffset = 0;
                yDegOffset = 0;
                xyzAngle[0] = 0;
                xyzAngle[1] = 0;
                // new Thread(thrd1).start();
                // new Thread(thrd2).start();
                // new Thread(thrd3).start();
                // new Thread(thrd4).start();
                // new Thread(thrd5).start();
                rayTracer.traceScene(camera, spheres, recursionDepth, canvas, xyzAngle, rtCWLower, rtCWHigher, meshes);
            }
        });

        rayTracer.traceScene(camera, spheres, recursionDepth, canvas, xyzAngle, rtCWLower, rtCWHigher, meshes);

        // new Thread(thrd1).start();
        // new Thread(thrd2).start();
        // new Thread(thrd3).start();
        // new Thread(thrd4).start();
        // new Thread(thrd5).start();
        System.out.println("Finishing drawing pixels"); // debug helper line
    }
}
