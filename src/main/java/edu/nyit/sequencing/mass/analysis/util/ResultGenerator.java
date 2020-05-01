package edu.nyit.sequencing.mass.analysis.util;

import edu.nyit.sequencing.mass.analysis.AnalysisHelper;
import edu.nyit.sequencing.mass.analysis.FindSequence;
import edu.nyit.sequencing.mass.analysis.model.AnchorNode;
import edu.nyit.sequencing.mass.analysis.model.DraftRead;
import edu.nyit.sequencing.mass.analysis.model.Pair;
import edu.nyit.sequencing.mass.analysis.model.TreeNode;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Dr. Wenjia Li
 *
 *
 */

public class ResultGenerator {
    List<DraftRead> finalReads;


    //This is the file path to output the sequencing result 
    
    public static final String directory = "Result/";
    
    public static final String fasta_directory = "fastaResult/";


    public ResultGenerator() {
        this.finalReads = new ArrayList<>();
    }

    public void writeToFile(AnchorNode anchor){
        HashMap<String, Integer> sizeMap = new HashMap<>();

        //First delete any existing file from previous
        //data analysis

        double target_mass1 = anchor.getMass();

        File folder = new File(directory);
        if (!folder.exists()) {
        	folder.mkdir();
        }
        
        for(File f: folder.listFiles()) {
            f.delete();
        }

        //Also delete any existing fasta file(s) from
        //previous data analysis, added 6/3/2019

        File fasta_folder = new File(fasta_directory);
        if (!fasta_folder.exists()) {
        	fasta_folder.mkdir();
        }
        
        for(File f1: fasta_folder.listFiles()) {
            f1.delete();
        }


        //Try to find the draft sequence with highest read length and also highest average volume

        int maxReadLength = -1;
        double maxVolume = -1.0;


        //Find the max read length

        for(DraftRead draftRead:finalReads){
            if (draftRead.getLength() > maxReadLength){
                maxReadLength = draftRead.getLength();
            }
        }

        for (DraftRead draftRead : finalReads) {
            if (draftRead.getLength() == maxReadLength) {

                //calculate the average volume for the draft reads with max read length

                if (draftRead.getAverageVolume() > maxVolume) {
                    maxVolume = draftRead.getAverageVolume();
                }
            }
        }


        //This list is used to store the best sequence(s)
        //The global hierarchical ranking strategy is
        //used to determine the final sequence from this list

        List<DraftRead> bestSequences = new ArrayList<>();

        System.out.println("Performing Global Hierarchical Ranking Strategy.\n");

        //Step 1: We first add the draft read(s) with the longest
        //read length and highest average volume to the list.

        for (DraftRead read : finalReads) {
            if (read.getLength() == maxReadLength && read.getAverageVolume() == maxVolume) {
                bestSequences.add(read);
            }
        }

        System.out.println("The num of best sequences after considering read length and avg volume is: " + bestSequences.size());

        //Step 2: If there are multiple draft reads that have the same
        //read length and same average volume, we will use the average QS
        //as the tie-breaker, and the draft read(s) with the highest average QS
        //will be the winner after Step 2.

        if(bestSequences.size() > 1){
            double maxQS = -1.0;

            for(DraftRead read:bestSequences){
                if (read.getAverageQS()> maxQS){
                    maxQS = read.getAverageQS();
                }
            }

            for(int i = 0; i < bestSequences.size(); i++){
                if (bestSequences.get(i).getAverageQS() != maxQS){
                    bestSequences.remove(i);
                }
            }

            System.out.println("The num of best sequences after further considering avg QS is: " + bestSequences.size());
        }



        //Step 3: If there are still multiple draft reads that have the same
        //read length, same average volume and same QS, we will use the average PPM
        //as the final tie-breaker, and the draft read(s) with the lowest average PPM
        //will be the winner after Step 3.

        if(bestSequences.size() > 1){
            double minPPM = 1000000.0;  //a large enough number to start with

            for(DraftRead read:bestSequences){
                if (read.getAverage_ppm() < minPPM){
                    minPPM = read.getAverage_ppm();
                }
            }

            for(int i = 0; i < bestSequences.size(); i++){
                if (bestSequences.get(i).getAverage_ppm() != minPPM){
                    bestSequences.remove(i);
                }
            }

            System.out.println("The num of best sequences after further considering avg PPM is: " + bestSequences.size());
        }

        System.out.println();


        Collections.sort(finalReads, Collections.reverseOrder());









        String fastaFile = "draftReads.fasta";

        int readCount = 1;

        //If we read from the 3'-Biotin labeled end, we should
        //put the reversed draft reads so that all of them are
        //consistently shown from the 5'-end.

        if(target_mass1 == 826.3184){
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fasta_directory + fastaFile)));

                for (DraftRead read : finalReads) {

                    bw.write(">SEQUENCE_" + readCount + "|Length: " + read.getLength() + "|Avg Volume: " + read.getAverageVolume());

                    String seq = "3'BiotinTag-" + read.getBases();
                    StringBuilder sb = new StringBuilder();
                    sb.append(seq);
                    //sb = sb.reverse();

                    bw.write("\n" + sb + "\n\n");
                    readCount += 1;
                }

                bw.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        else if(target_mass1 == 938.2216 || target_mass1 == 692.1105 || target_mass1 == 922.2269 || target_mass1 == 6398.1021 || target_mass1 == 668.0992 || target_mass1 == 877.1793) {
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fasta_directory + fastaFile)));

                for (DraftRead read : finalReads) {

                   bw.write(">SEQUENCE_" + readCount + "|Length: " + read.getLength()  + "|Avg Volume: " + read.getAverageVolume());

                   bw.write("\n" + read.getBases() + "\n\n");
                   readCount += 1;
                }

                bw.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        double similarityScore;


        String finalFasta = "finalReads.fasta";

        String isoFormFasta = "isoForms.fasta";

        int finalReadCount = 1;

        DraftRead model = null;

        if(finalReads.size() > 0){
            model = finalReads.get(0);
        }

        int match = 1;
        int mismatch = -1;
        int gap = -1;

        List<DraftRead> isoForms = new ArrayList<>();
        List<DraftRead> tRNAfinalISOForms = new ArrayList<>();




        if(target_mass1 == 826.3184 && model != null){
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fasta_directory + finalFasta)));

                isoForms.add(model);

                for (DraftRead read : finalReads) {

                    similarityScore = smithWatermanAlignment(model.getBases(), read.getBases(), match, mismatch, gap);

                    if (similarityScore >= 0.94) {

                        //System.out.println(read.getBases() + ": Similarity Score: " + similarityScore);

                        if(similarityScore < 1.0) {
                            isoForms.add(read);
                        }

                        bw.write(">ISOFORM_CANDIDATE_" + finalReadCount + "|Length: " + read.getComponent().size() + "|Avg Volume: " + read.getAverageVolume() + "|Smith-Waterman Alignment Similarity Score: " + similarityScore * 100 + "%");


                        String seq = "3'BiotinTag-" + read.getBases();
                        StringBuilder sb = new StringBuilder();
                        sb.append(seq);
                        //sb = sb.reverse();

                        bw.write("\n" + sb + "\n\n");
                        finalReadCount += 1;
                    }
                }

                bw.close();

            } catch (IOException e) {
                e.printStackTrace();
            }

            String modelBases = model.getBases();

            if(modelBases.startsWith("ACCACGC")){
                modelBases = modelBases.substring(3);
            }
            else if(modelBases.startsWith("CCACGC")){
                modelBases = modelBases.substring(2);
            }
            else if(modelBases.startsWith("CACGC")){
                modelBases = modelBases.substring(1);
            }
            else{
                System.out.println("ERROR!!! Model Sequence: " + modelBases + "does NOT belong to tRNA 3 iso forms!");
            }

            if(isoForms.size() > 0){
                for(int i = 0; i < isoForms.size(); i++){
                    //tRNA 76 form
                    if(isoForms.get(i).getBases().startsWith("ACCACGC")){

                        similarityScore = smithWatermanAlignment(modelBases, isoForms.get(i).getBases().substring(3), match, mismatch, gap);

                        if(similarityScore == 1.0){
                            tRNAfinalISOForms.add(isoForms.get(i));
                        }
                    }
                    //tRNA 75 form
                    else if(isoForms.get(i).getBases().startsWith("CCACGC")){

                        similarityScore = smithWatermanAlignment(modelBases, isoForms.get(i).getBases().substring(2), match, mismatch, gap);

                        if(similarityScore == 1.0){
                            tRNAfinalISOForms.add(isoForms.get(i));
                        }
                    }
                    //tRNA 74 form
                    else if(isoForms.get(i).getBases().startsWith("CACGC")){

                        similarityScore = smithWatermanAlignment(modelBases, isoForms.get(i).getBases().substring(1), match, mismatch, gap);

                        if(similarityScore == 1.0){
                            tRNAfinalISOForms.add(isoForms.get(i));
                        }
                    }
                    else{
                        System.out.println("ERROR!!! Sequence: " + isoForms.get(i).getBases() + "does NOT belong to tRNA 3 iso forms!");
                    }
                }
            }

            int isoFormCount = 1;

            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fasta_directory + isoFormFasta)));

                for(DraftRead read:tRNAfinalISOForms){
                    bw.write(">tRNA_isoform_" + isoFormCount + "|Length: " + read.getComponent().size() + "|Avg Volume: " + read.getAverageVolume());

                    String seq = "3'BiotinTag-" + read.getBases();
                    StringBuilder sb = new StringBuilder();
                    sb.append(seq);
                    //sb = sb.reverse();

                    bw.write("\n" + sb + "\n\n");
                    isoFormCount += 1;
                }

                bw.close();

            }
            catch (IOException e) {
                e.printStackTrace();
            }


            //We will print out the CCA truncated isoforms
            //only if all three are available and found.

            if(tRNAfinalISOForms.size() == 3){

                System.out.println("We have detected CCA truncated isoforms!");
                for(int i = 0; i < tRNAfinalISOForms.size(); i++){
                    System.out.println("CCA truncated isoform #" + (i+1) + ": " + tRNAfinalISOForms.get(i).getBases());
                }

                System.out.println();
            }



        }




        for(DraftRead draftRead:finalReads){

            double ppm_diff = (Math.abs(draftRead.getHead().mass_value-target_mass1)/target_mass1)*1000000;

            similarityScore = smithWatermanAlignment(model.getBases(), draftRead.getBases(), match, mismatch, gap);


            //Only write multiple draft sequences if it is
            //for 3'-Biotin labeled tRNA samples

            if (ppm_diff <= 10.0 && similarityScore >=0.94 && target_mass1 == 826.3184 && tRNAfinalISOForms.contains(draftRead))  {

                //print out ALL results to individual txt files
                sizeMap.put(String.valueOf(draftRead.getLength()), sizeMap.getOrDefault(String.valueOf(draftRead.getLength()), 0) + 1);

                String fname = String.valueOf(draftRead.getLength()) + "_" + draftRead.getAverage_ppm() + "_" + draftRead.getAverageVolume() + "_" + String.valueOf(sizeMap.get(String.valueOf(draftRead.getLength()))) + ".txt";
                try {
                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(directory + fname)));

                    //bw.write(draftRead.getHead().theoretic_mass_value + "\t" + draftRead.getHead().mass_value + "\t" + draftRead.getHead().retention_time + "\t\t" + draftRead.getHead().volume + "\n");
                    bw.write(draftRead.getHead().mass_value + "\t" + draftRead.getHead().retention_time + "\t\t" + draftRead.getHead().volume + "\n");

                    for (int i = 0; i < draftRead.getComponent().size(); i++) {
                        //bw.write(draftRead.getComponent().get(i).node.theoretic_mass_value + "\t" + draftRead.getComponent().get(i).node.mass_value + "\t" + draftRead.getComponent().get(i).node.retention_time + "\t" + draftRead.getComponent().get(i).base + "\t" + draftRead.getComponent().get(i).node.volume + "\t" + draftRead.getComponent().get(i).ppm + "\n");
                        bw.write(draftRead.getComponent().get(i).node.mass_value + "\t" + draftRead.getComponent().get(i).node.retention_time + "\t" + draftRead.getComponent().get(i).base + "\t" + draftRead.getComponent().get(i).node.volume + "\t" + String.format( "%.2f", draftRead.getComponent().get(i).ppm) + "\n");
                    }

                    bw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }

            }
            else if ((target_mass1 == 692.1105 || target_mass1 == 922.2269) && draftRead.getLength() == maxReadLength){
                sizeMap.put(String.valueOf(draftRead.getLength()), sizeMap.getOrDefault(String.valueOf(draftRead.getLength()), 0) + 1);

                String fname = "Length_" + String.valueOf(draftRead.getLength()) + "_AvgVolume_" + draftRead.getAverageVolume() + "_#" + String.valueOf(sizeMap.get(String.valueOf(draftRead.getLength()))) + ".txt";

                try {
                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(directory + fname)));

                    //bw.write("Mass of\t\tRT\t\tBase\tVolume\tPPM\nLadders\n");
                    bw.write(String.format("%8s\t%10s\t\t%6s\t%s\t%9s\t%4s\n","Fragment","Mass","RT","Base","Volume","PPM"));

                    if(target_mass1 == 443.0243){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tpG\t\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else if(target_mass1 == 826.3184){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tTag\t\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else if(target_mass1 == 922.2269){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tTag+A\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else if(target_mass1 == 938.2216){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tTag+G\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" +String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else if(target_mass1 == 692.1105 || target_mass1 == 612.1442){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tA+G\t\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else if(target_mass1 == 877.1793){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tA+C+C\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else if(target_mass1 == 6398.1021){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tMod+Psi\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else if(target_mass1 == 694.2397){
                        bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) + "\t3'Tag\t" + String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }
                    else {
                    	bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) + "\t" + anchor.getName() + "\t" + String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                    }

                    for (int i = 0; i < draftRead.getComponent().size(); i++) {

                        if(draftRead.getComponent().get(i).base.endsWith("end") || draftRead.getComponent().get(i).base.equals("U+Cm") || draftRead.getComponent().get(i).base.equals("A+Gm") || draftRead.getComponent().get(i).base.equals("ModPsi")){
                            bw.write(String.format("%8s", (i + 2)) + "\t" + String.format("%10s", String.format("%.4f", draftRead.getComponent().get(i).node.mass_value)) + "\t\t" + String.format("%6s", String.format("%.3f", draftRead.getComponent().get(i).node.retention_time)) + "\t" + draftRead.getComponent().get(i).base + "\t" + String.format("%9s", draftRead.getComponent().get(i).node.volume) + "\t" + String.format("%.2f", draftRead.getComponent().get(i).ppm) + "\n");
                        }else {
                            bw.write(String.format("%8s", (i + 2)) + "\t" + String.format("%10s", String.format("%.4f", draftRead.getComponent().get(i).node.mass_value)) + "\t\t" + String.format("%6s", String.format("%.3f", draftRead.getComponent().get(i).node.retention_time)) + "\t" + draftRead.getComponent().get(i).base + "\t\t" + String.format("%9s", draftRead.getComponent().get(i).node.volume) + "\t" + String.format("%.2f", draftRead.getComponent().get(i).ppm) + "\n");
                        }
                    }

                    bw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }


            //If we are reading from the 3'-Biotin labeled end and
            //all the three iso forms have been identified, we should
            //select and print out the final read from within the three
            //iso forms.

            if(target_mass1 == 826.3184 && tRNAfinalISOForms.size() == 3){
                double maxISOFormVolume = -1.0;

                String threeISOFormOutput = "Three CCA truncated isoforms\n\n";

                int isoFormCount = 1;

                for (DraftRead read : tRNAfinalISOForms) {
                    //calculate the average volume for the draft reads within three iso forms

                    if (read.getAverageVolume() > maxISOFormVolume) {
                        maxISOFormVolume = read.getAverageVolume();
                    }

                    threeISOFormOutput += "isoform #" + isoFormCount + "\n\n";

                    threeISOFormOutput += String.format("%8s\t%10s\t\t%6s\t%s\t%9s\t%4s\n","Fragment","Mass","RT","Base","Volume","PPM");

                    threeISOFormOutput += String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", read.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", read.getHead().retention_time)) +"\tTag\t\t"+ String.format("%9s",read.getHead().volume) + "\t" + String.format( "%.2f", read.getHead().getPPM()) + "\n";

                    for (int i = 0; i < read.getComponent().size(); i++) {
                        if(read.getComponent().get(i).base.endsWith("end") || read.getComponent().get(i).base.equals("U+Cm") || read.getComponent().get(i).base.equals("A+Gm")){
                            threeISOFormOutput += String.format("%8s", (i + 2)) + "\t" + String.format("%10s", String.format("%.4f", read.getComponent().get(i).node.mass_value)) + "\t\t" + String.format("%6s", String.format("%.3f", read.getComponent().get(i).node.retention_time)) + "\t" + read.getComponent().get(i).base + "\t" + String.format("%9s", read.getComponent().get(i).node.volume) + "\t" + String.format("%.2f", read.getComponent().get(i).ppm) + "\n";
                        }else {
                            threeISOFormOutput += String.format("%8s", (i + 2)) + "\t" + String.format("%10s", String.format("%.4f", read.getComponent().get(i).node.mass_value)) + "\t\t" + String.format("%6s", String.format("%.3f", read.getComponent().get(i).node.retention_time)) + "\t" + read.getComponent().get(i).base + "\t\t" + String.format("%9s", read.getComponent().get(i).node.volume) + "\t" + String.format("%.2f", read.getComponent().get(i).ppm) + "\n";
                        }
                    }

                    threeISOFormOutput += "\n\n";

                    isoFormCount++;

                }

                String finalOutput = "finalOutput.txt";

                try {
                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(directory + finalOutput)));

                    bw.write(threeISOFormOutput);

                    bw.close();

                }catch (IOException e) {
                    e.printStackTrace();
                }

                if (draftRead.getAverageVolume() == maxISOFormVolume){
                    String finalReadFile = "finalRead.txt";

                    String finalReadSeq = draftRead.getBases();

                    System.out.println("The best draft sequence is: " + finalReadSeq);


                    //Original code to produce the finalRead that can be used to draw the figure

                    try {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(directory+finalReadFile)));

                        bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tTag\t"+ draftRead.getHead().volume + "\n");

                        for (int i = 0; i < draftRead.getComponent().size(); i++) {
                            bw.write(draftRead.getComponent().get(i).node.mass_value+"\t"+draftRead.getComponent().get(i).node.retention_time+"\t"+draftRead.getComponent().get(i).base+"\t"+draftRead.getComponent().get(i).node.volume + "\t" + String.format( "%.2f", draftRead.getComponent().get(i).ppm)+"\n");
                        }

                        bw.close();

                    }catch (IOException e) {
                        e.printStackTrace();
                    }
                }


            }
            else{

                if (draftRead.getBases().equals(bestSequences.get(0).getBases()) && draftRead.getAverageVolume() == bestSequences.get(0).getAverageVolume()){

                    String finalReadFile = "finalRead.txt";

                    String finalReadSeq = draftRead.getBases();

                    System.out.println("The best draft sequence is: " + finalReadSeq);





                    //Original code to produce the finalRead that can be used to draw the figure

                    try {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(directory+finalReadFile)));



                        if(target_mass1 == 443.0243){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tpG\t"+ draftRead.getHead().volume + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 826.3184){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tTag\t"+ draftRead.getHead().volume + "\n");
                        }
                        else if(target_mass1 == 922.2269 || target_mass1 == 922.2268){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tTag+A\t"+ draftRead.getHead().volume + "\n");
                        }
                        else if(target_mass1 == 938.2216){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tTag+G\t"+ draftRead.getHead().volume + "\n");
                        }
                        else if(target_mass1 == 692.1105 || target_mass1 == 612.1442){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tA+G\t"+ draftRead.getHead().volume + "\n");
                        }
                        else if(target_mass1 == 877.1793){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tA+C+C\t"+ draftRead.getHead().volume + "\n");
                        }
                        else if(target_mass1 == 1225.3243 || target_mass1 == 6398.1021){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tCMC-Psi\t"+ draftRead.getHead().volume + "\n");
                        }
                        else if(target_mass1 == 668.0992 || target_mass1 == 668.0993){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\tG+C\t"+ draftRead.getHead().volume + "\n");
                        }
                        else if(target_mass1 == 694.2397){
                            bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\t3'Tag\t"+ draftRead.getHead().volume + "\n");
                        }
                        else {
                        	bw.write(draftRead.getHead().mass_value+"\t"+draftRead.getHead().retention_time +"\t" + anchor.getName() + "\t"+ draftRead.getHead().volume + "\n");
                        }

                        for (int i = 0; i < draftRead.getComponent().size(); i++) {
                            bw.write(draftRead.getComponent().get(i).node.mass_value+"\t"+draftRead.getComponent().get(i).node.retention_time+"\t"+draftRead.getComponent().get(i).base+"\t"+draftRead.getComponent().get(i).node.volume + "\t" + String.format( "%.2f", draftRead.getComponent().get(i).ppm)+"\n");
                        }

                        bw.close();

                    }catch (IOException e) {
                        e.printStackTrace();
                    }

                    //New code to produce the finalOutput file that is better aligned for human readers, added 5/28/2019
                    String finalOutput = "finalOutput.txt";

                    try {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(directory+finalOutput)));

                        //bw.write("Mass of\t\tRT\t\tBase\tVolume\tPPM\nLadders\n");
                        bw.write(String.format("%8s\t%10s\t\t%6s\t%s\t%9s\t%4s\n","Fragment","Mass","RT","Base","Volume","PPM"));

                        if(target_mass1 == 443.0243){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tpG\t\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 826.3184){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tTag\t\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 922.2269 || target_mass1 == 922.2268){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tTag+A\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 938.2216){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tTag+G\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" +String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 692.1105 || target_mass1 == 612.1442){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tA+G\t\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 877.1793){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tA+C+C\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 1225.3243 || target_mass1 == 6398.1021){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tMod-Psi\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 668.0992 || target_mass1 == 668.0993){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) +"\tG+C\t\t"+ String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 694.2397){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) + "\t3'Tag\t" + String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }
                        else if(target_mass1 == 898.2156){
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) + "\t5'Tag\t" + String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }

                        //For a new tag
                        else{
                            bw.write(String.format("%8s","1") + "\t" + String.format("%10s",String.format( "%.4f", draftRead.getHead().mass_value)) + "\t\t" + String.format("%6s", String.format( "%.3f", draftRead.getHead().retention_time)) + "\t" + anchor.getName() + "\t" + String.format("%9s",draftRead.getHead().volume) + "\t" + String.format( "%.2f", draftRead.getHead().getPPM()) + "\n");
                        }

                            for (int i = 0; i < draftRead.getComponent().size(); i++) {

                                if(draftRead.getComponent().get(i).base.endsWith("end") || draftRead.getComponent().get(i).base.equals("U+Cm") || draftRead.getComponent().get(i).base.equals("A+Gm")){
                                    bw.write(String.format("%8s", (i + 2)) + "\t" + String.format("%10s", String.format("%.4f", draftRead.getComponent().get(i).node.mass_value)) + "\t\t" + String.format("%6s", String.format("%.3f", draftRead.getComponent().get(i).node.retention_time)) + "\t" + draftRead.getComponent().get(i).base + "\t" + String.format("%9s", draftRead.getComponent().get(i).node.volume) + "\t" + String.format("%.2f", draftRead.getComponent().get(i).ppm) + "\n");
                                }else {
                                    bw.write(String.format("%8s", (i + 2)) + "\t" + String.format("%10s", String.format("%.4f", draftRead.getComponent().get(i).node.mass_value)) + "\t\t" + String.format("%6s", String.format("%.3f", draftRead.getComponent().get(i).node.retention_time)) + "\t" + draftRead.getComponent().get(i).base + "\t\t" + String.format("%9s", draftRead.getComponent().get(i).node.volume) + "\t" + String.format("%.2f", draftRead.getComponent().get(i).ppm) + "\n");
                                }
                            }

                            bw.close();

                    }catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }



        }
    }

    //This method implements the Smith-Waterman Alignment algorithm
    //added on 6/6/2019
    public double smithWatermanAlignment(String seq1, String seq2, int match, int mismatch, int gap){

        //First get the length of two sequences
        int seq1Len = seq1.length() + 1;
        int seq2Len = seq2.length() + 1;

        //Also convert the String to char[]
        char[] sequence1 = seq1.toCharArray();
        char[] sequence2 = seq2.toCharArray();

        //Store the resulting seq1 and seq2 after
        //alignment and before reverse
        String res1 = "";
        String res2 = "";

        String finalSeq1 = "";
        String finalSeq2 = "";

        String pipe = "";

        int[][] score = new int[seq2Len][seq1Len];

        HashMap<String, String> map = new HashMap<>();


        for(int i = 0; i < seq2Len; i++){
            for(int j = 0; j < seq1Len; j++) {
                score[i][j] = 0;
                //System.out.print(score[i][j] + " ");
            }
            //System.out.println();
        }



        int maxScore = 0;
        int sim,top,left,k=0,l=0;
        String key = "";
        String backTrackKey = "";

        //Filling the scoring matrix
        for(int m = 1; m < seq2Len; m++){
            for(int n = 1; n < seq1Len; n++){
                if(sequence1[n-1] == sequence2[m-1]){
                    sim = score[m-1][n-1] + match;
                }
                else {
                    sim = score[m-1][n-1] + mismatch;
                }
                top = score[m-1][n] + gap;
                left = score[m][n-1] + gap;

                score[m][n] = Math.max(sim, Math.max(top, Math.max(left, 0)));

                if(score[m][n] > maxScore){
                    maxScore = score[m][n];
                    //System.out.println("max score is: " + maxScore);
                    k = m;
                    l = n;
                }

                key = (m + "_" + n);

//                if(!map.containsKey(key)){
//                    map.
//                }

                if(score[m][n] == 0){
                    map.put(key,"STOP");
                }
                else if(score[m][n] == sim){
                    map.put(key,"sim");
                }
                else if(score[m][n] == top){
                    map.put(key,"top");
                }
                else if(score[m][n] == left){
                    map.put(key,"left");
                }

            }
        }

        //Backtracking
        int x = k;
        int y = l;

        while((x > 0) && (y > 0)){
            backTrackKey = (x + "_" + y);

            if(map.containsKey(backTrackKey)){
                if(map.get(backTrackKey).equals("STOP")){
                    break;
                }
                else if(map.get(backTrackKey).equals("sim")){
                    x -= 1;
                    y -= 1;
                    res1 += sequence1[y];
                    res2 += sequence2[x];
                }
                else if(map.get(backTrackKey).equals("top")){
                    x -= 1;
                    res1 += "-";
                    res2 += sequence2[x];
                }
                else if(map.get(backTrackKey).equals("left")){
                    y -= 1;
                    res1 += sequence1[y];
                    res2 += "-";
                }
            }
        }

        StringBuilder build1 = new StringBuilder();
        build1.append(res1);
        finalSeq1 += build1.reverse();

        StringBuilder build2 = new StringBuilder();
        build2.append(res2);
        finalSeq2 += build2.reverse();

        int comparison = Math.max(finalSeq1.length(), finalSeq2.length());

        for (int p = 0; p < comparison; p++){
            if (finalSeq1.charAt(p) == finalSeq2.charAt(p)){
                pipe += "|";
            }
            else {
                pipe += ".";
            }
        }

        //Print out the alignment result
//        System.out.println("Alignment Result:\n\n");
//        System.out.println(finalSeq1);
//        System.out.println(pipe);
//        System.out.println(finalSeq2);
//        System.out.println("\nAlignment Score: " + maxScore);

        double similarityScore = (double)maxScore / Math.max(seq1.length(), seq2.length());

//        System.out.println("\nSimilarity Score: " + similarityScore*100 + "%");

        return similarityScore;

    }

    public void generateResult(List<TreeNode> mass_data, AnchorNode anchor, int k) {
        Map<Integer, List<DraftRead>> lenGroups = new HashMap<>();
        FindSequence fs = new FindSequence();

        double anchor_mass = -1.0;

        if(anchor != null) {
            anchor_mass = anchor.getMass();
        }
        //step 1: for loop to populate this hashmap
        for (TreeNode treeNode : mass_data) {
            for (List<Pair> seq : treeNode.seq_list) {


                String start_end = "";

                Boolean addDraftRead = true;

                for(int i = 0; i < seq.size()-1;i++) {
                    if (seq.get(i).base.equals("U+Cm") && !seq.get(i+1).base.equals("A+Gm")) {
                        addDraftRead = false;
                    }

                    if (seq.get(i).base.equals("2mG") && seq.get(i+1).base.equals("Y'")) {
                        addDraftRead = false;
                    }

                    if (seq.get(i).base.equals("2mG") && seq.get(i+1).base.equals("G")) {
                        addDraftRead = false;
                    }

                    if (seq.get(i).base.equals("mA") && seq.get(i+1).base.equals("D")) {
                        addDraftRead = false;
                    }


                    if ((anchor_mass == 877.1793) && (seq.get(i).base.equals("Y'") || seq.get(i).base.equals("U+Cm") || seq.get(i).base.equals("A+Gm"))){
                        addDraftRead = false;
                    }

                    if(anchor_mass == 922.2269 && seq.get(i).base.equals("mA")){
                        addDraftRead = false;
                    }
                    
                    if(seq.get(i).base.endsWith("-end")) {
                    	addDraftRead = false;
                    }

                }




                //Remove the extra base after the correct sequence
                //in the case of RNA mixture reading
                //added 6/9/2019 and updated on 10/29/2019
                if(anchor_mass == 826.3184 || anchor_mass == 694.2397) {
                    double volumeRatio = (double) seq.get(seq.size() - 1).volume / seq.get(seq.size() - 2).volume;

                    if (volumeRatio < 1.0) {
                        //System.out.println("The old seq is: " + fs.seqToStr(seq));
                        seq.remove(seq.size() - 1);
                        //System.out.println("The new seq is: " + fs.seqToStr(seq));
                    }

//                    double newVolRatio = (double) seq.get(seq.size() - 1).volume / seq.get(seq.size() - 2).volume;
//
//                    if (volumeRatio < 10.0) {
//                        addDraftRead = false;
//                    }

                }



//                if(anchor_mass == 694.2397){
//                    if(seq.get(seq.size() - 1).volume < 1.0E7){
//                        addDraftRead = false;
//                    }
//                }

                //3'-Biotin labeled sequence should not use
                //any 5'-ending base mass value
                //added 11/6/2019

                if(anchor_mass == 694.2397){
                    if(seq.get(seq.size()-1).base.endsWith("end")){
                        addDraftRead = false;
                    }

                    if(seq.get(seq.size()-1).volume < 1.0E7){
                        addDraftRead = false;
                    }

                    if(seq.get(seq.size()-1).volume > 1.0E7 && seq.get(seq.size()-2).volume > 1.0E7){
                        addDraftRead = false;
                    }
                }

//                if(anchor_mass == 694.2397 && seq.size() == 21 && !seq.get(0).base.equals("A")){
//                    addDraftRead = false;
//                }

                if(anchor.getName().startsWith("5'Cy3")){
                    if(seq.get(seq.size()-1).base.endsWith("end") || seq.get(seq.size()-1).base.endsWith("end)")){
                        addDraftRead = true;
                    }else{
                        addDraftRead = false;
                    }

                    double volumeR = (double) seq.get(seq.size() - 1).volume / seq.get(seq.size() - 2).volume;

                    if(volumeR < 1.0 && seq.get(seq.size() - 2).volume > 1.0E7){
                        seq.remove(seq.size() - 1);
//                        System.out.println("Now the raw seq is: " + fs.seqToStr(seq));
                    }

                    for(int i = 0; i < seq.size()-1;i++){
                        if(seq.get(i).base.endsWith("end") || seq.get(i).base.endsWith("end)")){
                            addDraftRead = false;
                        }
                    }

                }

                if(anchor_mass == 922.2269){
                    for(int i = 0; i < seq.size()-1;i++){
                        if(seq.get(i).ppm > 10.0){
                            addDraftRead = false;
                        }
                    }
                }


                int length = seq.size();


                if (anchor_mass == 692.1105){
                    for(int i = 0; i < seq.size()-1;i++) {
                        if (seq.get(i).base.equals("U+Cm") && seq.get(i+1).base.equals("A+Gm") && (i >= 5)){
                            length += 2;
                        }
                    }
                }

                //calculate the average ppm for each draft read
                double total_ppm = 0.0;
                double avg_ppm = 0.0;

                //For recalculation of PPM, should be used
                //if using the adjacent mass diff to do base calling

                seq.get(0).node.theoretic_mass_value = anchor_mass + seq.get(0).baseMass;

                for (int i = 1; i < seq.size(); i++){
                    seq.get(i).node.theoretic_mass_value = seq.get(i-1).node.theoretic_mass_value + seq.get(i).baseMass;
                }

                for (int i = 0; i < seq.size(); i++){
                    double dm = Math.abs(seq.get(i).node.theoretic_mass_value - seq.get(i).node.mass_value);
                    double ppm = AnalysisHelper.dm2ppm(seq.get(i).node.theoretic_mass_value, dm);

                    if(ppm <= 10.0) {
                        seq.get(i).setPPM(ppm);
                    }
                    //Disabled 10/23/2019 for 2 mixture RNA readings, may need to enable for tRNA phenol
                    else{
                        if(anchor_mass == 694.2397 && seq.size() == 21) {
                            addDraftRead = false;
                        }
                    }
                }


                for (int i = 0; i < seq.size(); i++){
                    total_ppm += seq.get(i).ppm;
                }

                avg_ppm = total_ppm / seq.size();

                //calculate the average volume for each draft read

                int total_volume = treeNode.volume;
                double avg_volume = 0.0;

                for (int i = 0; i < seq.size(); i++){
                    total_volume += seq.get(i).volume;
                }

                avg_volume = total_volume / (seq.size()+1);

                //calculate the average QS for each draft read

                double total_QS = treeNode.quality_score;
                double avg_QS = 0.0;

                for (int i = 0; i < seq.size(); i++){
                    total_QS += seq.get(i).quality_score;
                }

                avg_QS = total_QS / (seq.size()+1);


                DraftRead draftRead = null;

                if(addDraftRead == true) {
                    draftRead = new DraftRead(length, fs.seqToStr(seq), start_end, avg_ppm, avg_volume, avg_QS, treeNode, seq); //has to create the object
                }
                //int len = length;

                boolean dubDraftRead = false;

                if(draftRead != null) {
                    if (!lenGroups.containsKey(length)) {
                        lenGroups.put(length, new ArrayList<>());
                    }

                    for(DraftRead d:lenGroups.get(length)){
                        if(d.getBases().equals(draftRead.getBases())){
                            if(d.getAverageVolume() >= draftRead.getAverageVolume()){
//                                System.out.println("There is already one draft read " + d.getBases() + " that has same or higher avg vol! Skip.");
                                dubDraftRead = true;
                            }
                        }
                    }
                    if(dubDraftRead == false) {
                        lenGroups.get(length).add(draftRead);
                    }
                }
            }
        }

        //step 2:
        for (Map.Entry<Integer, List<DraftRead>> entry : lenGroups.entrySet()) {
            List<DraftRead> draftReads = entry.getValue();
            List<DraftRead> afterFilter = simReadFilter(draftReads, k);
            finalReads.addAll(afterFilter);
        }

        writeToFile(anchor);
    }


    public List<DraftRead> simReadFilter(List<DraftRead> draftReads, int k) {
        List<DraftRead> afterFilter = new ArrayList<>();
        Map<String, List<DraftRead>> simGroups = new HashMap<>();
        //step 1: for loop to populate this hashmap to group by start_end
        for (DraftRead draftRead : draftReads) {

            if (!simGroups.containsKey(draftRead.getStart_end())) {
                simGroups.put(draftRead.getStart_end(), new ArrayList<>());
            }
            simGroups.get(draftRead.getStart_end()).add(draftRead);
        }

        //step 2: find top k min
        for (Map.Entry<String, List<DraftRead>> entry : simGroups.entrySet()) {
            PriorityQueue<DraftRead> tops = findTopKMin(entry.getValue(), k);
            //TODO: add to afterFilter
            for (DraftRead top : tops) {


                afterFilter.add(top);
            }
        }
        return afterFilter;
    }


    public  PriorityQueue<DraftRead> findTopKMin(List<DraftRead> draftReads, int k) {
        PriorityQueue<DraftRead> top = new PriorityQueue<>();
        for (DraftRead draftRead : draftReads) {
            top.add(draftRead);
            if (top.size() > k) {
                top.remove();
            }
        }
        //top stores top k min
        return top;
    }


}
