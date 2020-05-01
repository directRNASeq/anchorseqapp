package edu.nyit.sequencing.mass.analysis.model;

import java.util.ArrayList;
import java.util.List;

/*
To record a draft read (a path found by DFS)
 */

public class DraftRead implements Comparable{
    //instance variables
    private Integer length;
    private String bases;
    private String start_end; //start_end
    private Double average_ppm;
    private Double average_volume;
    private Double average_QS;
    private TreeNode head;
    private List<Pair> component;

    public DraftRead(int length, String bases, String start_end, double avg_ppm, TreeNode head, List<Pair> component){
        this.length = length;
        this.bases = bases;
        this.start_end = start_end;
        this.average_ppm = avg_ppm;
        this.head = head;
        this.component = component;
    }

    public DraftRead(int length, String bases, String start_end, double avg_ppm, double avg_volume, TreeNode head, List<Pair> component){
        this.length = length;
        this.bases = bases;
        this.start_end = start_end;
        this.average_ppm = avg_ppm;
        this.average_volume = avg_volume;
        this.head = head;
        this.component = component;
    }

    public DraftRead(int length, String bases, String start_end, double avg_ppm, double avg_volume, double avg_QS, TreeNode head, List<Pair> component){
        this.length = length;
        this.bases = bases;
        this.start_end = start_end;
        this.average_ppm = avg_ppm;
        this.average_volume = avg_volume;
        this.average_QS = avg_QS;
        this.head = head;
        this.component = component;
    }

    public String getStart_end() {
        return start_end;
    }

    public TreeNode getHead(){
        return this.head;
    }

    public int getLength(){return this.length;}

    public List<Pair> getComponent(){
        List<Pair> copy = new ArrayList<>();

        for (int i = 0; i < this.component.size(); i++){
            copy.add(this.component.get(i));
        }

        return copy;
    }

    public Double getAverage_ppm(){
        return this.average_ppm;
    }

    public String getBases(){return this.bases;}

    public Double getAverageVolume(){
        return this.average_volume;
    }

    public Double getAverageQS(){
        return this.average_QS;
    }

    @Override
    public int compareTo(Object other) {
        return this.average_volume.compareTo(((DraftRead) other).average_volume);
    }
}
