package ListMF.LPMF;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

import org.ejml.simple.SimpleMatrix;


public class test {
	public static void main(String[] args) throws Exception {
		//DataModel model = new FileDataModel(new File("/home/manni/Google Drive/research/code/ds/ml-100k/u.train"));
		//System.out.println(model.getItemIDsFromUser(1).toString());
		//System.out.println(model.getPreferencesFromUser(1).toString());
		//MF_SGDTest m = new MF_SGDTest();
		//m.MovielensTest();
		run();					
	}
	
	public static void run() throws Exception{
		
		//String trainSet = "/home/manni/Google Drive/research/code/ds/ml-100k/u1.base";  //Discarded 
		String trainSet = "src/main/java/u1.base";
		//String testSet = "/home/manni/Google Drive/research/code/ds/ml-100k/u1.test";
		//String testSet = "src/main/java/u1.test";
        //String resFile = "/home/manni/Google Drive/research/code/ds/MF_SGD_RESULT.txt";
        String resFile = "src/main/java/MF_SGD_RESULT.txt";
        
        BufferedWriter out = new BufferedWriter(new FileWriter(resFile));
        
        SimpleMatrix R,RP;  //RP is the matrix after normalization. 
        
        R  = MatrixMaker.createMatrix(trainSet, 943, 1682);
        RP=MatrixMaker.getPMatrix(R);
        RP.saveToFileCSV("RP2");
        System.exit(0);
        
        R  = MatrixMaker.createMatrix("src/main/java/RP2", 943, 1682);   // This is for supplying NORMALIZED matrix to the SGD  
        double alpha[] = new double[943]; 
        Arrays.fill(alpha,1.0);
        SGD mf = new SGD(R,10,alpha,0.1);  // all goes to SGD
        
               
        out.close();       
	}

}
