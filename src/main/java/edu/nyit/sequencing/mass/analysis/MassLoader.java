package edu.nyit.sequencing.mass.analysis;
import edu.nyit.sequencing.mass.analysis.model.MassNode;

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


public class MassLoader {
	public List<MassNode> massNodes;
	
	public MassLoader() {
		this.massNodes = new ArrayList<MassNode>();
	}
	
	public void loadData(){

        //This is the file path for base_bank.csv. 

        String fileName = "config/base_bank.csv";

        try {
        	
           BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
        	
           String line = null;
         
           while((line = br.readLine())!= null){
                String[] splitArray =line.replace("\n", "").split("\t");
                MassNode massNode = new MassNode(Double.parseDouble(splitArray[1]), splitArray[0]);
                this.massNodes.add(massNode);
           }
            
            
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        //Arrays.sort(mass_data);
    }
	
	public static void main(String[] args) {
		MassLoader test = new MassLoader();
		test.loadData();

		for(MassNode m:test.massNodes){
		    System.out.println("Name: "+ m.name + " Mass: " + m.mass);
        }

		System.out.println("ok");
	}
}
