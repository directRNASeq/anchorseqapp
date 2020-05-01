package edu.nyit.sequencing.mass.analysis.model;

/**
 * @author Dr. Wenjia Li
 *
 *
 */

public class Pair {
    public TreeNode node;
    public String base;
    //added 6/18/2019
    public double baseMass;
    //added 6/8/2019
    public double retention_time;
    //added 1/14/2019
    public double ppm;
    //added 4/23/2019
    public int volume;
    //added 5/5/2019
    public double quality_score;

    public Pair(TreeNode node, String base) {
        this.node = node;
        this.base = base;
    }

    //added 1/14/2019
    public Pair(TreeNode node, String base, double ppm) {
        this.node = node;
        this.base = base;
        this.ppm = ppm;
    }

    //added 4/23/2019
    public Pair(TreeNode node, String base, double ppm, int volume) {
        this.node = node;
        this.base = base;
        this.ppm = ppm;
        this.volume = volume;
    }

    //added 5/5/2019
    public Pair(TreeNode node, String base, double ppm, int volume, double qs) {
        this.node = node;
        this.base = base;
        this.ppm = ppm;
        this.volume = volume;
        this.quality_score = qs;
    }

    //added 6/18/2019
    public Pair(TreeNode node, String base, double baseMass, double ppm, int volume, double qs) {
        this.node = node;
        this.base = base;
        this.baseMass = baseMass;
        this.ppm = ppm;
        this.volume = volume;
        this.quality_score = qs;
    }

    public void setPPM(double p){
        this.ppm = p;
    }

    //copy constructor, added 3/5/2019, updated 5/5/2019
    public Pair(Pair other){
        this.node = other.node;
        this.base = other.base;
        this.ppm = other.ppm;
        this.volume = other.volume;
        this.quality_score = other.quality_score;
    }

    @Override
    public String toString() {
        return node.toString();
    }

    @Override
    public int hashCode() {
        return this.node.mass_value.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        return (this.node.mass_value == ((Pair)obj).node.mass_value);
    }
}
