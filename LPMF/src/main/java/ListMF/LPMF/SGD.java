package ListMF.LPMF;


import java.util.Random;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.ejml.simple.SimpleMatrix;


public class SGD {
	
	public SimpleMatrix W,Wt,WtH,H;
	public static SimpleMatrix V;
	public double lambdaW,lambdaH,max_err;
	public double[] alpha;
	public int n,m,r,step,max_step,WtH_n,WtH_m,W_n,W_m,H_n,H_m,Wt_n,Wt_m;
	
	public double getCov(SimpleMatrix M){
		Covariance c = new Covariance();
		double[][] row = new double[M.numRows()][M.numCols()];
		double[] col = new double[M.numCols()];
		for(int ui=1;ui<M.numRows(); ui++){
    		for(int uj=1;uj<M.numCols();uj++){
    			col[uj-1]=M.get(ui,uj);
    		}
    		row[ui-1]=col;
		}
		double cov = 0;
		for(int i=0;i<row.length;i++){
			cov+=c.covariance(row[i], row[i+1]);
		}
		return cov/row.length;
	}
    
	public SGD(SimpleMatrix R, int r, double[] alpha,double max_err) {
        this.n = R.numRows();
        this.m = R.numCols();
        this.r = r; //Matrix factorization factor. 
        this.step = 0;
        //this.max_step = step; //donno what it was
        this.max_err = max_err;
        this.alpha = alpha;
        
        this.W = new SimpleMatrix(SimpleMatrix.random(n,r,0.1,0.9,new Random()));
        this.H = new SimpleMatrix(SimpleMatrix.random(r,m,0.1,0.9,new Random()));
        V = R;
        //this.lambdaW = getCov(W)/getCov(V);
        this.lambdaW = 1;   
        //this.lambdaH = getCov(H)/getCov(V);
        this.lambdaH = 1;
        this.Wt = W.transpose();
        this.W_n = W.numRows();
    	this.W_m = W.numCols();
    	this.H_n = H.numRows();
    	this.H_m = H.numCols();
    	this.Wt_n = Wt.numRows();
    	this.Wt_m = Wt.numCols();
    }
	
	public static double doMultiplyVectors(double[] V1,double[] V2){
		double sum=0;
		for(int i=0;i<V1.length;i++){
			sum+=V1[i]*V2[i];
		}
		return sum;
	}
	
	public double LF(){
		// this is the output
    	double E;
    	// sum of g where g(i) is not 0.
    	double[] H_sum = new double[H_n];
    	for(int r=1;r<H_n;r++){
    		double sum=0;
    		for(int c=1;c<H_m;c++){
    			if(V.get(r,c)>0){
    				sum+=H.get(r,c);
    			}  			
    		}
    		H_sum[r]=sum;
    	} 
    	//Main loop for part1
    	double P1 = 0;
    	double[] Wt_vectors = new double[Wt_n];
    	double[] H_vectors = new double[H_n];
    	for(int r=1;r<Wt_n;r++){
    		//Forming vectors
    		for(int c=1;c<Wt_m;c++){
    			Wt_vectors[c]=Wt.get(r,c);
    		}
    		for(int c=1;c<H_m;c++){
    			H_vectors[c]=H.get(r,c);
    		}
    		//main part 1 equation.. 
    		for(int i=1; i<n; i++){
        		for(int j=1; j<m; j++){
        			if(V.get(i, j)>0) {
        				P1+=V.get(i,j)*(Math.log(doMultiplyVectors(Wt_vectors,H_vectors))-Math.log(doMultiplyVectors(H_sum,Wt_vectors)));
        			}
        		}
        	}
    	}

    	P1 = -P1;
    	
    	//Determinant of W
    	double W_det = 0;	
    	for(int i=1;i<W_m;i++){
    		for(int j=1;j<W_n;j++){
    			W_det+=Math.pow(W.get(i,j),2);
    		}
    	}
    	//Determinant of H
    	double H_det = 0;	
    	for(int i=1;i<H_m;i++){
    		for(int j=1;j<H_n;j++){
    			H_det+=Math.pow(H.get(i,j),2);
    		}
    	}
    	double P2=((lambdaW/2)*W_det)+((lambdaH/2)*H_det);
    	
    	double P3=0;
    	for(int r=1;r<W_n;r++){
    		P3+=alpha[r]*(1-doMultiplyVectors(H_sum,Wt_vectors));
    	}
    	
    	E=P1+P2+P3;
    	return E;
    }
	
	public void GDA(){
		double P2=lambdaW*W.determinant()+getMatrixSum(H);
		for(int vi=1; vi<n; vi++){
			for(int vj=1;vj<m;vj++){
    	    	if(V.get(vi,vj)>0){
    	    	double P1=V.get(vi,vj)*((H.mult(WtH.invert())).determinant() - (getMatrixSumCons(H)/getMatrixSumCons(WtH)));
    	    	P1=-P1;
    	    	W.set(vi, vj, P1+P2);
    	    	}
			}
		}
    }
	
	public static double getMatrixSumCons(SimpleMatrix M){
		double sum = 0;
		for(int ui=1;ui<M.numRows(); ui++){
    		for(int uj=1;uj<M.numCols();uj++){
    			if(V.get(ui,uj)>0){
    				sum+=M.get(ui, uj);
    			}
    		}
		}
		return sum;
	}
	
	public static double getMatrixSum(SimpleMatrix M){
		double sum = 0;
		for(int ui=1;ui<M.numRows(); ui++){
    		for(int uj=1;uj<M.numCols();uj++){
    			sum+=M.get(ui, uj);
    		}
		}
		return sum;
	}
		
	public void GDB(){
		double P2=lambdaH*H.determinant()+getMatrixSum(W);
		for(int vi=1; vi<n; vi++){
			for(int vj=1;vj<m;vj++){
    	    	if(V.get(vi,vj)>0){
    	    	double P1=V.get(vi,vj)*((W.mult(WtH.invert())).determinant() - (W.determinant()/getMatrixSumCons(WtH)));
    	    	P1=-P1;
    	    	H.set(vi, vj, P1+P2);
    	    	}
			}
		}
    }
	
	public void GDC(){
	    	
	}
		
	}
