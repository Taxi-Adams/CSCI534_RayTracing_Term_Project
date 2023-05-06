public class Light {
    String type;
    double intensity;
    double[] position;
    double[] direction;

    public Light(String type, double intensity, double[] position, double[] direction) {
        this.type = type;
        this.intensity = intensity;
        if(position != null) {
            this.position = position;
        }
        if(direction != null) {
            this.direction = direction;
        }
    }
}
