package edu.nyit.sequencing.mass.analysis;


import edu.nyit.sequencing.mass.analysis.model.AnchorLoader;
import edu.nyit.sequencing.mass.analysis.model.MassNode;
import edu.nyit.sequencing.mass.analysis.model.AnchorNode;
import edu.nyit.sequencing.mass.analysis.model.Pair;
import edu.nyit.sequencing.mass.analysis.model.TreeNode;
import edu.nyit.sequencing.mass.analysis.util.ResultGenerator;

import java.io.*;
import java.util.*;

/**
 * @author Dr. Wenjia Li
 *
 *
 */

public class FindSequence {

    public static final double cluster_ppm = 10.0;
    public static final int rt_multiplier = 1;
    public static final int max_cluster_depth = 1;


    public double anchor_mass = -1.0;

    public AnchorNode anchor;


    public List<TreeNode> mass_data;

    public FindSequence() {
        mass_data = new ArrayList<>();
    }

    public List<TreeNode> loadData(String fileName){

        List<TreeNode> initialData = new ArrayList<>();

        //process and read in the input LC/MS MFE data if it is in the original csv file, added 2/15/2019

        if(fileName.endsWith(".csv")){
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));

                //try to read the first line, which is the column header
                String line = br.readLine();

                String[] columnHeader = line.replace("\n", "").split(",");

                //find the indices for columns that contain Mass, RT, and Vol data
                //can also add more columns of data for processing if needed
                int massIndex = -1, rtIndex = -1, volumeIndex = -1, cpdIndex = -1, widthIndex = -1;

                for(int i = 0; i < columnHeader.length; i++){
                    if(columnHeader[i].equalsIgnoreCase("Mass"))
                        massIndex = i;
                    else if(columnHeader[i].equalsIgnoreCase("RT"))
                        rtIndex = i;
                    else if(columnHeader[i].equalsIgnoreCase("Vol"))
                        volumeIndex = i;
                    else if(columnHeader[i].equalsIgnoreCase("Cpd"))
                        cpdIndex = i;
                    else if(columnHeader[i].equalsIgnoreCase("width"))
                        widthIndex = i;
                }

                //keep reading the columns for Mass, RT, and Vol and store them into the list of TreeNode
                while((line = br.readLine())!= null){
                    String[] splitArray = line.replace("\n", "").split(",");
                    //System.out.println(splitArray[massIndex] + " " + splitArray[rtIndex] + " " + splitArray[volumeIndex]);
                    TreeNode rawData = new TreeNode(Double.parseDouble(splitArray[massIndex]), Double.parseDouble(splitArray[rtIndex]), Integer.parseInt(splitArray[volumeIndex]), splitArray[cpdIndex], Double.parseDouble(splitArray[widthIndex]));
                    initialData.add(rawData);
                }
                br.close();

            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }


        //process and read in the input LC/MS MFE data if it has been stored as a .txt file, revised 2/15/2019

        if(fileName.endsWith(".txt")){
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
                String line = null;
                while((line = br.readLine())!= null){
                    String[] splitArray = line.replace("\n", "").split("\t");
//                    TreeNode rawData = new TreeNode(Double.parseDouble(splitArray[0]), Double.parseDouble(splitArray[1]), Integer.parseInt(splitArray[2]), splitArray[3], Double.parseDouble(splitArray[4]), Double.parseDouble(splitArray[5]));
                    TreeNode rawData = new TreeNode(Double.parseDouble(splitArray[0]), Double.parseDouble(splitArray[1]), Integer.parseInt(splitArray[2]), splitArray[3], Double.parseDouble(splitArray[4]));
//                    TreeNode rawData = new TreeNode(Double.parseDouble(splitArray[0]), Double.parseDouble(splitArray[1]), Integer.parseInt(splitArray[2]), splitArray[3], Double.parseDouble(splitArray[4]));
                    initialData.add(rawData);
                }
                br.close();

            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }


        Collections.sort(initialData);

        return initialData;
    }

    //Find LC/MS data points that are near Mass mass and RT rt within a mass ppm and RT rtMultiplier range
    //Added on 4/9/2019
    public List<TreeNode> findByMassRT(List<TreeNode> compounds, double mass, double rt, double ppm, int rtMultiplier){
        List<TreeNode> result = new ArrayList<>();

        double massDelta = AnalysisHelper.ppm2dm(mass, ppm);

        for (TreeNode t:compounds){
            if (((t.mass_value >= mass - massDelta) && (t.mass_value <= mass + massDelta)) && ((t.retention_time >= rt  - t.width*rtMultiplier) && (t.retention_time <= rt + t.width*rtMultiplier)))
                //System.out.println("Found some close data points");
                result.add(t);
        }

        return result;

    }

    public List<TreeNode> recursiveCluster(TreeNode t, List<TreeNode> compounds, List<MassNode> adducts, int depth){
        List<TreeNode> result = new ArrayList<>();
        List<TreeNode> candidates = new ArrayList<>();
        List<TreeNode> moreCandidates = new ArrayList<>();

        for(MassNode a:adducts){
            result = findByMassRT(compounds, t.mass_value.doubleValue()+a.mass,t.retention_time,cluster_ppm,rt_multiplier);

            for(int i = 0; i < result.size(); i ++){
                if(result.get(i).Cpd.equals(t.Cpd)){
                    result.remove(result.get(i));
                    continue;
                }
                else{
                    candidates.add(result.get(i));

                    List<TreeNode> altcompounds = new ArrayList<>();

                    for(TreeNode c:compounds){
                        if(!c.Cpd.equals(result.get(i).Cpd)){
                            altcompounds.add(c);
                        }
                    }

                    if (depth < max_cluster_depth){
                        moreCandidates = recursiveCluster(result.get(i),altcompounds,adducts,depth+1);
                        if(moreCandidates.size()>0){
                            candidates.addAll(moreCandidates);
                        }
                    }

                }
            }

        }

        return candidates;
    }

    public List<TreeNode> removeAdducts(List<TreeNode> compounds, List<MassNode> adducts){
        List<TreeNode> doneCompounds = new ArrayList<>();

        for(int i = 0; i < compounds.size();i++){
            if(doneCompounds.contains(compounds.get(i))){
                continue;
            }

            List<TreeNode> adductsFound = new ArrayList<>();
            adductsFound = recursiveCluster(compounds.get(i), compounds,adducts,0);

            if (!adductsFound.contains(compounds.get(i))){
                adductsFound.add(compounds.get(i));
            }

            int maxVolume = -1;

            for(TreeNode tn:adductsFound){
                if(tn.volume > maxVolume)
                    maxVolume = tn.volume;
            }

            TreeNode parent = null;

            for(TreeNode tn1:adductsFound){
                if(tn1.volume == maxVolume)
                    parent = new TreeNode(tn1.mass_value,tn1.retention_time,tn1.volume,tn1.Cpd,tn1.width,tn1.quality_score);
            }

            if(adductsFound.size() > 1 && parent != null){
                System.out.println("Found " + (adductsFound.size()-1) + " adducts for node with mass: " + compounds.get(i).mass_value);

                for(TreeNode tn2:adductsFound){
                    if(tn2.Cpd.equals(parent.Cpd)){
                        doneCompounds.add(tn2);
                        continue;
                    }

                    parent.volume += tn2.volume;
                    compounds.remove(tn2);
                    doneCompounds.add(tn2);
                }

                //compounds.add(parent);

                compounds.get(i).volume = parent.volume;


//                System.out.println("Mass node " + compounds.get(i).mass_value + " new volume is: " + compounds.get(i).volume);

            }

        }

        for(TreeNode t:compounds){
            this.mass_data.add(t);
        }

        Collections.sort(this.mass_data);

        return compounds;
    }

    public boolean isValidMatch(double ppmValue, int location, double baseValue){
        double ppmLimit = (10.0 * mass_data.get(location).mass_value)/baseValue;



        if (ppmValue <= ppmLimit){
            return true;
            //} else if (mass_data.get(location).mass_value > 11555.0 && ppmValue <= ppmLimit_large){
            //    return true;
        } else{
            return false;
        }

    }



    public Pair matchMassValue(int i, int j, double diff, double mass_value, List<MassNode> massNodes){
        for(MassNode massNode:massNodes) {
            //if (isValidMatch(Math.abs(AnalysisHelper.ppmCal(diff, massNode.getMass())), i, massNode.getMass())){

            double ppm_diff = (Math.abs(diff-massNode.getMass())/mass_value)*1000000;

            if(this.anchor_mass == 668.0992) {
                if (ppm_diff <= 10.9) {
                    return new Pair(mass_data.get(j), massNode.getName(), massNode.getMass(), ppm_diff, mass_data.get(j).volume, mass_data.get(j).quality_score);
                }
            }
            else{
                if (ppm_diff <= 13.0) {
                    return new Pair(mass_data.get(j), massNode.getName(), massNode.getMass(), ppm_diff, mass_data.get(j).volume, mass_data.get(j).quality_score);
                }
            }
        }
        return null;
    }


    //An alternative way to determine whether or not there is a matched pair
    //in the base calling step using the theoretic mass value to calculate the PPM
    //added 4/21/2019
    public Pair matchMassValueAlt(int i, int j, double starting_theoretic_mass, double observed, List<MassNode> massNodes){
        for(MassNode massNode:massNodes){

            double theoretic_mass = starting_theoretic_mass + massNode.mass;
            double ppm_diff = (Math.abs(theoretic_mass - observed)/theoretic_mass)*1000000;

            if(ppm_diff <= 13.0){
                return new Pair(mass_data.get(j), massNode.getName(), massNode.getMass(), ppm_diff, mass_data.get(j).volume, mass_data.get(j).quality_score);
            }
        }
        return null;
    }

    public static Boolean isBase(String base){
        Boolean result = false;
        if (base.equals("A") || base.equals("G") || base.equals("C") || base.equals("U")){
            result = true;
        }
        return result;
    }

    //Method for base calling
    public void pairDiff(List<MassNode> massNodes, List<AnchorNode> anchorNodes){
        int count = 0;

        //Try to first find the anchor node in the MFE dataset
        //Here we just assume that there is only one
        //anchor node in the MFE dataset

        for(AnchorNode anchor:anchorNodes) {
            for (TreeNode t : mass_data) {
                double delta_mass = AnalysisHelper.ppm2dm(anchor.getMass(), cluster_ppm);

                if ((t.mass_value >= anchor.getMass() - delta_mass) && (t.mass_value <= anchor.getMass() + delta_mass)) {

                    t.theoretic_mass_value = anchor.getMass();
                    this.anchor_mass = anchor.getMass();
                    this.anchor = new AnchorNode(anchor.getName(), anchor.getMass());

                    System.out.println("***Base Calling***\n\nWe found the anchor for this MFE data as: " + this.anchor.getName() + " and its mass is: " + this.anchor.getMass()
                                        + "\n\nThe result of base calling:");

                    break;
                }
            }
        }

        for (int i = 0; i < mass_data.size()-1; i++){
            for (int j = i+1; j < mass_data.size(); j++){

                double mass_diff = mass_data.get(j).mass_value - mass_data.get(i).mass_value;
                double time_diff = mass_data.get(j).retention_time - mass_data.get(i).retention_time;
                Pair matchedValue = null;


                if(this.anchor != null && this.anchor.getMass() == 694.2397) {
                    if (mass_data.get(i).theoretic_mass_value == this.anchor.getMass()) {
                        if (mass_data.get(i).theoretic_mass_value > 0.0 && mass_diff > 0.0) {
                            matchedValue = matchMassValueAlt(i, j, mass_data.get(i).theoretic_mass_value, mass_data.get(j).mass_value, massNodes);
                        }
                    } else {
                        if (mass_data.get(i).theoretic_mass_value > 0.0 && mass_diff > 0.0 && Math.abs(time_diff) <= 3.0) {
                            matchedValue = matchMassValueAlt(i, j, mass_data.get(i).theoretic_mass_value, mass_data.get(j).mass_value, massNodes);
                        }
                    }
                }
                else{
                    if (mass_data.get(i).theoretic_mass_value > 0.0 && mass_diff > 0.0) {
                        matchedValue = matchMassValueAlt(i, j, mass_data.get(i).theoretic_mass_value, mass_data.get(j).mass_value, massNodes);
                    }
                }

                if(matchedValue != null && this.anchor.getName().startsWith("3'") && matchedValue.base.contains("end")){
                    matchedValue = null;
                }

                if(matchedValue != null && this.anchor.getMass() == 694.2397 && !isBase(matchedValue.base)){
                    matchedValue = null;
                }


                if((matchedValue != null) && (this.anchor_mass != 1225.3243 && this.anchor_mass != 6398.1021) && (mass_data.get(i).theoretic_mass_value == this.anchor_mass) && !FindSequence.isBase(matchedValue.base)){
                    matchedValue = null;
                }

                if ((matchedValue != null) && (this.anchor_mass == 826.3184 || this.anchor_mass == 443.0243 || this.anchor_mass == 938.2216) && (matchedValue.base.equals("Y'") || matchedValue.base.equals("U+Cm") || matchedValue.base.equals("A+Gm"))){
                    matchedValue = null;
//                    mass_data.get(j).theoretic_mass_value = -1;
                }

                if((matchedValue != null) && (this.anchor_mass == 443.0243) && (matchedValue.base.equals("mC") || matchedValue.base.equals("T") || matchedValue.base.equals("Y") || matchedValue.base.equals("mA"))){
                    matchedValue = null;
//                    mass_data.get(j).theoretic_mass_value = -1;
                }

                if(matchedValue != null){
                    mass_data.get(j).theoretic_mass_value = mass_data.get(i).theoretic_mass_value + matchedValue.baseMass;
                }


                //Added 1/10/2019 to detect if there is already any child in the list
                //that has the same base with the newly found pair to be added
                //Updated 4/22/2019 and now we will replace the child with a higher volume if needed
                boolean replace = false;
                boolean containPair = false;


                //check to see if another pair with the same base has
                //already been added.

                if(matchedValue != null){
                    for(Pair pair:mass_data.get(i).children){
                        if(pair.base.equals(matchedValue.base)){
                            containPair = true;
                        }
                    }
                }

                //check to see if the new pair has the same base
                //with one of the children, but it has a higher
                //volume. If so, replace the old one with this
                //new pair
                if (matchedValue != null){
                    for (int x = 0; x < mass_data.get(i).children.size(); x++){
                        if ((mass_data.get(i).children.get(x).base.equals(matchedValue.base)) && (mass_data.get(i).children.get(x).node.volume < matchedValue.volume) ){

                            mass_data.get(i).children.remove(x);
                            replace = true;


                        }
                    }
                }


                if ((matchedValue != null) && ((replace == true) || (containPair == false)) && (mass_data.get(i).children.size() < 4)){

                    mass_data.get(i).children.add(matchedValue);

                    //print out the details for the base calling step
                    //added 11/23/2019
                    System.out.println("Base Calling #" + (count+1) + ": " + mass_data.get(i).mass_value + "-" + matchedValue.base + "-" + matchedValue.node.mass_value);

                    count++;
                }

            }
        }
        System.out.println("\nThe total number of base calling in the MFE data is: " + count + "\n");
    }

    public void findSequence(){
        int numOfPath = 0;
        HashSet<TreeNode> visited = new HashSet<>();
        for (TreeNode treeNode : mass_data) {

            double ppm_diff = (Math.abs(treeNode.mass_value-this.anchor_mass)/this.anchor_mass)*1000000;

            if (ppm_diff > 10.0){
                continue;
            }

            if (visited.contains(treeNode)) {
                System.out.println("Visited already for node with mass: "+treeNode.mass_value);
                continue;
            }

//            System.out.println("Started to process mass node: "+treeNode.mass_value);

            Set<List<Pair>> res = new HashSet<>();
            HashMap<String, Integer> seqDB = new HashMap<>();
            findSequence_DFS(treeNode, new ArrayList<>(), res, seqDB);
            int maxlen =  getMaxLength(res);

            for (List<Pair> r : res) {

                //If the dataset is for CMC conversion, then we will
                //allow draft sequence with short read length
                //revised 7/11/2019

                if(this.anchor_mass == 1225.3243 || this.anchor_mass == 6398.1021 || this.anchor_mass == 4348.7897 || this.anchor_mass == 668.0992 || this.anchor_mass == 877.1793) {

                    if (r.size() >= 4) {//try to limit the output sequence amount
                        treeNode.seq_list.add(r);
                    }
                }
                else if(this.anchor.getMass() == 443.0243 || this.anchor.getMass() == 938.2216 || this.anchor_mass == 668.0993){

                    if (r.size() >= 17) {//try to limit the output sequence amount
                        treeNode.seq_list.add(r);
                    }
                }
                else{
                    if (r.size() >= 18) {//try to limit the output sequence amount
                        treeNode.seq_list.add(r);
                    }
                }


                for (Pair item : r) {
                    visited.add(item.node);
                }
            }
            numOfPath += treeNode.seq_list.size();

        }
        System.out.println("***Draft Sequence Generation and Final Sequence(s) Identification***\n");

        ResultGenerator rg = new ResultGenerator();

        //for 3'-Biotin label, we may want to find
        //3 isoforms for tRNA
        if(this.anchor_mass == 826.3184){
            rg.generateResult(this.mass_data, this.anchor, 40000);
        }
        else if (this.anchor_mass == 877.1793){
            rg.generateResult(this.mass_data, this.anchor, 40000);
        }
        else if (this.anchor_mass == 612.1442){
            rg.generateResult(this.mass_data, this.anchor, 5);
        }
        else if(this.anchor_mass == -1.0){
            System.out.println("Could not find any anchor node in the MFE Data. STOP HERE!");
        }
        else{
            rg.generateResult(this.mass_data, this.anchor, 2);
        }


    }

    //find all possible path from input node to leaf nodes
    public void findSequence_DFS(TreeNode treeNode, List<Pair> localSeq, Set<List<Pair>> res, HashMap<String, Integer> seqDB) {
        //mark the current node
        treeNode.visited = true;


        List<Pair> children = treeNode.children;
        if (children.size() == 0) {

            if (localSeq.size() > 0) {
                List<Pair> newSeq = new ArrayList<>(localSeq);

                res.add(newSeq);
                treeNode.maxPathLength = newSeq.size();
                treeNode.visited = false;
            }
            return;
        }


        for (int i = 0; i < children.size(); i++) {
            Pair pair = children.get(i);
            if (!pair.node.visited) {
                localSeq.add(pair);
                findSequence_DFS(pair.node, localSeq, res, seqDB);
                localSeq.remove(pair);
            }
        }


        treeNode.visited = false;
    }


    public String seqToStr(List<Pair> seq) {
        StringBuilder sb = new StringBuilder();
        for (Pair s : seq) {
            sb.append(s.base);
        }
        return sb.toString();
    }


    public int getMaxLength(Set<List<Pair>> seq_list) {
        int max = 0;
        for (List<Pair> seq : seq_list) {
            if (max < seq.size()) {
                max = seq.size();
            }
        }
        return max;
    }

    public void testImage() {
        try {
            System.out.println("Start to generate image!");
            Process p = Runtime.getRuntime().exec("python webapps/plot_mass.py");
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("Failed to generate image!");
        }
    }

    public static void main(String[] args){
        FindSequence fs = new FindSequence();

        //Please note that a valid MFE data file in txt format with full file path is required
        //here to make the program run.

        String fileName = "data/input_data/TableS4_122718s07.txt";
        
        fs.mass_data.addAll(fs.loadData(fileName));

        System.out.println("***Data Pre-processing and Loading***\n\nThe size of MFE data is: " + fs.mass_data.size() + "\n");





        MassLoader massLoader = new MassLoader();
        massLoader.loadData();


        AnchorLoader anchorLoader = new AnchorLoader();
        anchorLoader.loadData();




        fs.pairDiff(massLoader.massNodes, anchorLoader.anchorNodes);

        fs.findSequence();
    }
}
