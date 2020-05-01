package edu.nyit.sequencing.mass.analysis.model;

public class MassNode {
	public double lower;
	public double higher;
	public double mass;
	public String name;
	
	
	public MassNode(double mass, String name) {
		super();
		this.mass = mass;
		this.name = name;
	}
	
	public double getLower() {
		return lower;
	}
	public void setLower(double lower) {
		this.lower = lower;
	}
	public double getHigher() {
		return higher;
	}
	public void setHigher(double higher) {
		this.higher = higher;
	}
	public double getMass() {
		return mass;
	}
	public void setMass(double mass) {
		this.mass = mass;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}

	public String toString(){
		String output = "";
		output += "Base Name: " + this.name + ", Base Mass: " + this.mass + "\n";
		return output;
	}
	
	
}
