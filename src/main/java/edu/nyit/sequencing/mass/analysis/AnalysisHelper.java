package edu.nyit.sequencing.mass.analysis;

/**
 * @author Wenjia Li
 *
 */

public class AnalysisHelper {
	public static double ppmCal(double diff, double base){
        //System.out.println("diff: "+diff+", base: "+base);
        double res = (Math.abs(diff-base)/base)*1000000;
        //System.out.println(res);
        return res;
    }

    public static double ppm2dm(double mass, double ppm){
	    return (mass * ppm)/ 1.0E6;
    }

    public static double dm2ppm(double mass, double dm){
	    return (dm * 1.0E6) / mass;
    }

}
