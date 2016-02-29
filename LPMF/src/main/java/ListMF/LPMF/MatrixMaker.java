package ListMF.LPMF;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


import org.ejml.simple.SimpleMatrix;


public class MatrixMaker {
	
	@SuppressWarnings("finally")
	public static SimpleMatrix createMatrix(String fileName, int numr, int numc) {
        SimpleMatrix R = new SimpleMatrix(numr+1, numc+1);
        File file = new File(fileName);
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(file));
            String str;
            while ((str = reader.readLine()) != null) {
                String[] ss = str.split("\t");
                R.set(Integer.parseInt(ss[0]),Integer.parseInt(ss[1]), Double.parseDouble(ss[2]));
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                }
            }
            return R;
        }
    }
	
	public static SimpleMatrix getPMatrix(SimpleMatrix R) throws Exception{
		int r_r = R.numRows();
		int r_c = R.numCols();
		SimpleMatrix RP = new SimpleMatrix(r_r, r_c);
		SimpleMatrix RPN = new SimpleMatrix(r_r, r_c);
		int rp_r = RP.numRows();
		int rp_c = RP.numCols();
		BufferedWriter out = new BufferedWriter(new FileWriter("src/main/java/RP"));
		
		for(int i=1;i<r_r;i++){
			for(int j=1;j<r_c;j++){						
				if(R.get(i,j) > 0){
					RP.set(i,j,Math.exp(R.get(i,j)));	
					continue;
				}
				RP.set(i,j,0);
			}
		}
		
		double[] rs = new double[944];				
		for(int row=1;row<rp_r;row++){
			double sum=0;
			for(int col=1;col<rp_c;col++){
				sum+=RP.get(row,col);						
			}
			rs[row]=sum;
		}
		
		for(int i=1;i<rp_r;i++){
			for(int j=1;j<rp_c;j++){						
				if(RP.get(i,j) > 0){
					RPN.set(i,j,RP.get(i,j)/rs[i]);	
					out.write(i+"\t"+j+"\t"+RPN.get(i,j)+"\t\n");
					continue;
				}
				RPN.set(i,j,0);
				out.write(i+"\t"+j+"\t"+0+"\n");
			}
		}				
		out.close();
		return RPN;
	}
	
}
