import java.util.ArrayList;

public class WingedEdgeMesh {
    private int[][] faceList;
    private double[][][] vertexList;
    private short[] color;
    private double specular;
    private double reflective;
    private double[] center;
    private Edge[] edges;
    private int[][] faceEdgeList;
    private int[][] vertexEdgeList;

    public WingedEdgeMesh(int[][] faceList, double[][][] vertexList, short[] color, double specular, double reflective, double[] center) {
        this.faceList = faceList;
        this.vertexList = vertexList;
        this.color = color;
        this.specular = specular;
        this.reflective = reflective;
        this.center = center;
        // inefficient, but for my sanity
        // creating edges
        ArrayList<Edge> edgeList = new ArrayList<>();
        for(int i = 0; i < faceList.length; i++) {
            edgeList.add(new Edge(new int[] {faceList[i][0], faceList[i][1]}));
            edgeList.add(new Edge(new int[] {faceList[i][1], faceList[i][2]}));
            edgeList.add(new Edge(new int[] {faceList[i][0], faceList[i][2]}));
        }
        for(int i = 0; i < edgeList.size(); i++) {
            Edge comparisonEdge = edgeList.get(i);
            int[] compEdgeVerts = comparisonEdge.getVertexIndexList();
            double[] a = vertexList[compEdgeVerts[0]][0];
            double[] b = vertexList[compEdgeVerts[1]][0];
            for(int j = 0; j < edgeList.size(); j++) {
                if(i != j) {
                    int[] jEdgeVerts = edgeList.get(j).getVertexIndexList();
                    double[] c = vertexList[jEdgeVerts[0]][0];
                    double[] d = vertexList[jEdgeVerts[1]][0];
                    if(equalPoints(a, c) && equalPoints(b, d)) {
                        edgeList.remove(j);
                    }
                }
            }
        }
        edges = new Edge[edgeList.size()];
        for(int i = 0; i < edgeList.size(); i++) {
            edges[i] = edgeList.get(i);
        }

        // face edge list
        faceEdgeList = new int[faceList.length][3];
        int faceIndex = 0;
        for(int i = 0; i < faceList.length; i++) {
            for(int j = 0 ; j < edges.length; j++) {
                int count = 0;
                for(int k = 0; k < 3; k++) {
                    if(faceList[i][k] == edges[j].getVertexIndexList()[0] || faceList[i][k] == edges[j].getVertexIndexList()[1]) {
                        count++;
                    }
                }
                if(count == 2) {
                    faceEdgeList[i][faceIndex] = j;
                    edges[j].addFaceToPair(i);
                    edges[j].addVertsToFaceVertList(faceList[i]);
                    faceIndex++;
                }
            }
            faceIndex = 0;
        }
    }

    public int[][] getFaceList() {
        return faceList;
    }

    public double[][][] getVertexList() {
        return vertexList;
    }

    public double[] getCenter() {
        return center;
    }

    public short[] getColor() {
        return color;
    }

    public double getSpecular() {
        return specular;
    }

    public double getReflective() {
        return reflective;
    }

    public Edge[] getEdges() {
        return edges;
    }

    public boolean equalPoints(double[] a, double[] b) {
        boolean isEqual = false;
        if(a[0] == b[0] && a[1] == b[1] && a[2] == b[2]) {
            isEqual = true;
        }
        return isEqual;
    }

    public int[][] getFaceEdgeList() {
        return faceEdgeList;
    }

    public int[][] getVertexEdgeList() {
        return vertexEdgeList;
    }

    public void setVertexEdgeList() {
        // System.out.println(edges.length);

        vertexEdgeList = new int[vertexList.length][];
        for(int i = 0; i < vertexList.length; i++) {
            ArrayList<Integer> edgesWithVertex = new ArrayList<>();
            for(int j = 0; j < edges.length; j++) {
                if(equalPoints(vertexList[edges[j].getVertexIndexList()[0]][0], vertexList[i][0]) || equalPoints(vertexList[edges[j].getVertexIndexList()[1]][0], vertexList[i][0])) {
                    edgesWithVertex.add(j);
                }
            }
            vertexEdgeList[i] = new int[edgesWithVertex.size()];
            for(int j = 0; j < edgesWithVertex.size(); j++) {
                vertexEdgeList[i][j] = edgesWithVertex.get(j);
            }
        }
    }

    public void setBoundaryEdges() {
        for(int i = 0; i < edges.length; i++) {
            // if a vertex is only connected to two edges, then the edge is a boundary edge
            if(vertexEdgeList[edges[i].getVertexIndexList()[0]].length == 2 || vertexEdgeList[edges[i].getVertexIndexList()[1]].length == 2) {
                edges[i].setIsBoundaryEdge(true);
            }
        }
    }
}
