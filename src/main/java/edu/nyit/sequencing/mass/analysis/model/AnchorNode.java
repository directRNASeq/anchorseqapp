package edu.nyit.sequencing.mass.analysis.model;

/*
This class is used to record the information
for Anchor, including its name, theoretical
mass value, etc.
5/10/2019
 */

/**
 * @author Dr. Wenjia Li
 *
 *
 */

public class AnchorNode {
    private String name;
    private double mass;

    public AnchorNode(String name, double mass){
        this.name = name;
        this.mass = mass;
    }

    public String getName(){return this.name;}
    public void setName(String n){this.name = n;}
    public double getMass(){return this.mass;}
    public void setMass(double m){this.mass = m;}

    public String toString(){
        String output = "";
        output += "Name: " + this.name + " Mass: " + this.mass + "\n";
        return output;
    }
}
