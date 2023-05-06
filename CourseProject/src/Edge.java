import java.util.ArrayList;

public class Edge {
    private int[] vertexIndexPair;
    private ArrayList<Integer> facePair;
    private ArrayList<int[]> faceVertList;
    private boolean isBoundaryEdge;

    public Edge(int[] vertexIndexPair) {
        this.vertexIndexPair = vertexIndexPair;
        facePair = new ArrayList<>();
        faceVertList = new ArrayList<>();
        isBoundaryEdge = false;
    }

    public int[] getVertexIndexList() {
        return vertexIndexPair;
    }

    public void addFaceToPair(int i) {
        facePair.add(i);
    }

    public void removeFaceFromPair(int i) {
        facePair.remove(Integer.valueOf(i));
    }

    public ArrayList<Integer> getFacePair() {
        return facePair;
    }

    public void addVertsToFaceVertList(int[] a) {
        faceVertList.add(a);
    }

    public ArrayList<int[]> getFaceVertList() {
        return faceVertList;
    }

    public int getFace1OtherVert() {
        int otherVert = -1;
        if(faceVertList.get(0)[0] != vertexIndexPair[0] && faceVertList.get(0)[0] != vertexIndexPair[1]) {
            otherVert = faceVertList.get(0)[0];
        } else if(faceVertList.get(0)[1] != vertexIndexPair[0] && faceVertList.get(0)[1] != vertexIndexPair[1]) {
            otherVert = faceVertList.get(0)[1];
        } else if(faceVertList.get(0)[2] != vertexIndexPair[0] && faceVertList.get(0)[2] != vertexIndexPair[1]) {
            otherVert = faceVertList.get(0)[2];
        }

        return otherVert;
    }

    public int getFace2OtherVert() {
        int otherVert = -1;
        if(faceVertList.get(1)[0] != vertexIndexPair[0] && faceVertList.get(1)[0] != vertexIndexPair[1]) {
            otherVert = faceVertList.get(1)[0];
        } else if(faceVertList.get(1)[1] != vertexIndexPair[0] && faceVertList.get(1)[1] != vertexIndexPair[1]) {
            otherVert = faceVertList.get(1)[1];
        } else if(faceVertList.get(1)[2] != vertexIndexPair[0] && faceVertList.get(1)[2] != vertexIndexPair[1]) {
            otherVert = faceVertList.get(1)[2];
        }

        return otherVert;
    }

    public boolean getIsBoundaryEdge() {
        return isBoundaryEdge;
    }

    public void setIsBoundaryEdge(boolean isOne) {
        isBoundaryEdge = isOne;
    }
}
