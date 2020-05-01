package edu.nyit.sequencing.mass.analysis.model;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TreeNode implements Comparable{
    public Double mass_value;
    public double theoretic_mass_value = -1.0;
    public Double retention_time;
    public int volume;
    public String Cpd;
    public Double width;
    public double quality_score;
    public boolean visited;
    public List<Pair> children;
    public Set<List<Pair>> seq_list;
    public int maxPathLength;

    public TreeNode(double mass_value, double retention_time) {
        this.children = new ArrayList<>();
        this.mass_value = mass_value;
        this.retention_time = retention_time;
        this.seq_list = new HashSet<>();
        this.visited = false;
        this.maxPathLength = 0;
    }

    public TreeNode(double mass_value, double retention_time, int volume) {
        this.children = new ArrayList<>();
        this.mass_value = mass_value;
        this.retention_time = retention_time;
        this.volume = volume;
        this.seq_list = new HashSet<>();
        this.visited = false;
        this.maxPathLength = 0;
    }

//    public TreeNode(double mass_value, double retention_time, int volume, String cpd, double width) {
//        this.children = new ArrayList<>();
//        this.mass_value = mass_value;
//        this.retention_time = retention_time;
//        this.volume = volume;
//        this.Cpd = cpd;
//        this.width = width;
//        this.seq_list = new HashSet<>();
//        this.visited = false;
//        this.maxPathLength = 0;
//    }

    public TreeNode(double mass_value, double retention_time, int volume, String cpd, double width, double qs) {
        this.children = new ArrayList<>();
        this.mass_value = mass_value;
        this.retention_time = retention_time;
        this.volume = volume;
        this.Cpd = cpd;
        this.width = width;
        this.quality_score = qs;
        this.seq_list = new HashSet<>();
        this.visited = false;
        this.maxPathLength = 0;
    }

    public TreeNode(double mass_value, double retention_time, int volume, double qs) {
        this.children = new ArrayList<>();
        this.mass_value = mass_value;
        this.retention_time = retention_time;
        this.volume = volume;
        this.quality_score = qs;
        this.seq_list = new HashSet<>();
        this.visited = false;
        this.maxPathLength = 0;
    }

    public TreeNode(double mass_value, double retention_time, int volume, String cpd, double qs) {
        this.children = new ArrayList<>();
        this.mass_value = mass_value;
        this.retention_time = retention_time;
        this.volume = volume;
        this.Cpd = cpd;
        this.quality_score = qs;
        this.seq_list = new HashSet<>();
        this.visited = false;
        this.maxPathLength = 0;
    }

    @Override
    public int compareTo(Object data1) {
        return this.mass_value.compareTo(((TreeNode) data1).mass_value);
    }

    @Override
    public String toString() {
        return this.Cpd;
    }

    @Override
    public int hashCode() {
        return this.Cpd.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        return this.Cpd == ((TreeNode)obj).Cpd;
    }

    public double getPPM(){
        if(this.theoretic_mass_value != -1){
            return Math.abs((this.theoretic_mass_value-this.mass_value)/this.theoretic_mass_value*1000000);
        }
        else
            return -1.0;
    }
}
