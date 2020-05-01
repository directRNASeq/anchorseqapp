package edu.nyit.sequencing.mass.analysis.model;

import edu.nyit.sequencing.mass.analysis.MassLoader;
import edu.nyit.sequencing.mass.analysis.model.AnchorNode;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Dr. Wenjia Li
 *
 *
 */

public class AnchorLoader {
    public List<AnchorNode> anchorNodes;

    public AnchorLoader() {
        this.anchorNodes = new ArrayList<AnchorNode>();
    }

    public void loadData(){

        //This is the file path for anchor_bank.csv

        String fileName = "config/anchor_bank.csv";

        try {
        	//URL resource = AnchorLoader.class.getClassLoader().getResource(fileName);
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String line = null;

            while((line = br.readLine())!= null){
                String[] splitArray =line.replace("\n", "").split("\t");
                AnchorNode anchorNode = new AnchorNode(splitArray[0], Double.parseDouble(splitArray[1]));
                this.anchorNodes.add(anchorNode);
            }


        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        //Arrays.sort(mass_data);
    }

    public static void main(String[] args) {
        AnchorLoader test = new AnchorLoader();
        test.loadData();

        for(AnchorNode anchor:test.anchorNodes){
            System.out.println(anchor.toString());
        }

        System.out.println("ok");

//        System.out.println(String.format("%s","Mass"));
    }
}
