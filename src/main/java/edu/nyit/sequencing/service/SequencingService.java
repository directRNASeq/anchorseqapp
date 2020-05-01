package edu.nyit.sequencing.service;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import javax.ws.rs.Consumes;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.core.Context;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import javax.ws.rs.core.UriInfo;

import com.sun.jersey.core.header.FormDataContentDisposition;
import com.sun.jersey.multipart.FormDataParam;
import edu.nyit.sequencing.mass.analysis.FindSequence;
import edu.nyit.sequencing.mass.analysis.MassLoader;
import edu.nyit.sequencing.mass.analysis.model.AnchorLoader;



/**
 * @author Dr. Wenjia Li
 *
 *
 */

@Path("/upload")
public class SequencingService {

    /** The path to the folder where we want to store the uploaded files */
    private static final String UPLOAD_FOLDER = "uploadedFiles/";
    private static final String RESULT_PLOT = "Result/";
    private static final String FASTA_RESULT = "fastaResult/";

    public SequencingService() {}

    @Context
    private UriInfo context;

    /**
     * Returns text response to caller containing current time-stamp
     * @return error response in case of missing parameters an internal exception or
     * success response if file has been stored successfully
     */
    @POST
    @Consumes(MediaType.MULTIPART_FORM_DATA)
    public Response uploadFile(
            @FormDataParam("file") InputStream uploadedInputStream,
            @FormDataParam("file") FormDataContentDisposition fileDetail) {

        // check if all form parameters are provided
        if (uploadedInputStream == null || fileDetail == null)
            return Response.status(400).entity("Invalid form data").build();

        // create our destination folder, if it not exists
        try {
            createFolderIfNotExists(UPLOAD_FOLDER);
            createFolderIfNotExists(RESULT_PLOT);
            createFolderIfNotExists(FASTA_RESULT);
        } catch (SecurityException se) {
            return Response.status(500).entity("Can not create destination folder on server").build();
        }

        String uploadedFileLocation = UPLOAD_FOLDER + fileDetail.getFileName();

        //First delete any existing file from previous data analysis

        File folder = new File(UPLOAD_FOLDER);

        for(File f: folder.listFiles())
            f.delete();

        /*try {
            Process p = Runtime.getRuntime().exec("del uploadedFiles/*.*");
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("Failed to delete the previous data file(s)");
        }*/

        try {
            saveToFile(uploadedInputStream, uploadedFileLocation);
            //sequencing
            FindSequence fs = new FindSequence();
            fs.mass_data.addAll(fs.loadData(uploadedFileLocation));

            System.out.println("***Data Pre-processing and Loading***\n\nThe size of MFE data is: " + fs.mass_data.size() + "\n");

            MassLoader test = new MassLoader();
            test.loadData();
            AnchorLoader anchorLoader = new AnchorLoader();
            anchorLoader.loadData();
            fs.pairDiff(test.massNodes, anchorLoader.anchorNodes);
            fs.findSequence();
            fs.testImage();

        } catch (IOException e) {
            return Response.status(500).entity("Can not save file").build();
        }

        return Response.status(200).entity("File saved to " + uploadedFileLocation).build();

    }

    /**
     * Utility method to save InputStream data to target location/file
     * @param inStream - InputStream to be saved
     * @param target - full path to destination file
     */
    private void saveToFile(InputStream inStream, String target) throws IOException {
        OutputStream out = null;
        OutputStream out1 = null;
        
        String lcmsFile = UPLOAD_FOLDER + "LCMS_data.txt";

        int read = 0;
        byte[] bytes = new byte[1024];
        
        System.out.println("Target path: "+target);
        
        System.out.println("LC-MS path: "+ lcmsFile);
        
        // create our destination folder, if it not exists
//        try {
//            createFolderIfNotExists("resultFiles/data/");
//        } catch (SecurityException se) {
//            throw se;
//        }

        out = new FileOutputStream(new File(target));
        out1 = new FileOutputStream(new File(lcmsFile));

        while ((read = inStream.read(bytes)) != -1) {
            out.write(bytes, 0, read);
            out1.write(bytes, 0, read);
        }
        
        out.flush();
        out.close();
        
        out1.flush();
        out1.close();
    }

    /**
     * Creates a folder to desired location if it not already exists
     * @param dirName - full path to the folder
     * @throws SecurityException - in case you don't have permission to create the folder
     */
    private void createFolderIfNotExists(String dirName) throws SecurityException {
        File theDir = new File(dirName);
        if (!theDir.exists()) {
            theDir.mkdir();
        }
    }

}
