package ListMF.LPMF;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.mahout.math.DenseMatrix;
import org.apache.mahout.math.Matrix;
import org.apache.mahout.math.SparseMatrix;

public class MFTestCommon {
	@SuppressWarnings("finally")
			public static Matrix createMatrix(String fileName, int numu, int numi) {
		        Matrix R = new SparseMatrix(numu+1, numi+1);
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
			
	
			public static Matrix getPMatrix(Matrix R) throws Exception{
				Matrix RP = new DenseMatrix(R.rowSize(), R.columnSize());
				Matrix RPN = new DenseMatrix(R.rowSize(), R.columnSize());
				BufferedWriter out = new BufferedWriter(new FileWriter("src/main/java/RP"));
				
				for(int i=1;i<R.rowSize();i++){
					for(int j=1;j<R.columnSize();j++){						
						if(R.get(i,j) > 0){
							RP.set(i,j,Math.exp(R.get(i,j)));	
							continue;
						}
						RP.set(i,j,0);
					}
				}
				
				double[] rs = new double[944];				
				for(int row=1;row<RP.rowSize();row++){
					double sum=0;
					for(int col=1;col<RP.columnSize();col++){
						sum+=RP.get(row,col);						
					}
					rs[row]=sum;
				}
				
				for(int i=1;i<RP.rowSize();i++){
					for(int j=1;j<RP.columnSize();j++){						
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
			
			public static double calMse(Matrix S, Matrix SS) {
		        int cnt = 0;
		        double num = 0;
		        for(int i=0; i<SS.rowSize(); i++) for(int j=0; j<SS.columnSize(); j++)
		            if(SS.get(i, j) != 0) {
		                cnt ++;
		                num += Math.pow(S.get(i, j)-SS.get(i, j), 2);
		            }
		        return num/cnt;
		    }
		
		    public static double calDensity(Matrix R) {
		        double cnt = 0;
		        for(int i=0; i<R.rowSize(); i++)
		            for(int j=0; j<R.columnSize(); j++)
		                if(R.get(i, j) != 0) cnt++;
		        return cnt/(R.rowSize()*R.columnSize());
		    }	
		
}
