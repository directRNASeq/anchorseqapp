package edu.nyit.sequencing.mass.analysis.model;

/**
 * @author Dr. Wenjia Li
 *
 *
 */

public class MatchedPair implements Comparable{
    public Integer start;
    public Integer end;
    public String base;
    /**
     * @param start
     * @param end
     * @param base
     */
    public MatchedPair(int start, int end, String base) {
        super();
        this.start = start;
        this.end = end;
        this.base = base;
    }
    
    @Override
    public int compareTo(Object data1) {
    	//TODO: what if they have the same start, we need to order by end
        return this.start.compareTo(((MatchedPair)data1).start);   
    }
    
    @Override
    public String toString() {
//        return "[start = " + this.start + ", end = " + this.end + ", based = " + this.base + "]";
//        return "[" + this.start + ", " + this.end + ", " + this.base + "]";
        return this.base;
    }

    
    
}
