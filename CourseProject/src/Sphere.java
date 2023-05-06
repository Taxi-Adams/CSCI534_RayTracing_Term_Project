public class Sphere {
    short[] color;
    double[] center;
    double radius;
    double specular;
    double reflective;

    public Sphere(short[] color, double[] center, double radius, double specular, double reflective) {
        this.color = color;
        this.center = center;
        this.radius = radius;
        this.specular = specular;
        this.reflective = reflective;
    }
}
