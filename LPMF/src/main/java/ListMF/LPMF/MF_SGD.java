package ListMF.LPMF;



import org.apache.mahout.math.DenseMatrix;
import org.apache.mahout.math.Matrix;
import org.apache.mahout.math.function.Functions;
import org.la4j.Matrices;
import org.la4j.inversion.GaussJordanInverter;

/**
 * Created by lmq on 2014/12/24.
 */
public class MF_SGD implements MF{

    /** Arrays for internal storage of V, W and H. */
    private final Matrix V;
    private final Matrix W,Wt;
    private final Matrix H;

    private final double lambda;
    private final double beta;

    /** Row and column dimensions. */
    private final int n;
    private final int m;

    /** number of features */
    private final int r;

    /** iteration step number before stop*/
    private int stepNum;

    /** max iteration steps */
    private final int stepMax;



    /** when fabs(oldObjectFunction-newObjectFunction) < errMax, stop iteration*/
    private final double errMax;
    
    public Matrix WH;


    public MF_SGD(Matrix arg, int r, double lambda, double beta, int stepNum, double errMax) {
        this.n = arg.numRows();
        this.m = arg.numCols();
        this.r = r; //Matrix factorization factor.
        this.stepNum = 0;
        this.stepMax = stepNum;
        this.errMax = errMax;
        this.lambda = lambda;
        this.beta = beta;

        this.W = new DenseMatrix(n, r).assign(Functions.random());
        this.H = new DenseMatrix(r, m).assign(Functions.random());
        this.V = arg;
        this.Wt = W.transpose();
    }

    //@Override
    public void solve() {
        for(; stepNum < stepMax; stepNum++) {
            for(int i=0; i<n; i++) for(int j=0; j<m; j++)
                if(V.get(i, j) != 0) {
                    double eij = V.get(i, j)-W.viewRow(i).dot(H.viewColumn(j));
                    for(int k=0; k<r; k++) {
                        W.set(i, k, W.get(i, k)+lambda*(2*eij*H.get(k, j)-beta*W.get(i, k)));
                        H.set(k, j, H.get(k, j)+lambda*(2*eij*W.get(i, k)-beta*H.get(k, j)));
                    }
                }

            double newObject = calObject(V, W, H);
            System.out.printf("%d %f %f\n", stepNum, newObject, Math.abs(object-newObject));
            if(Math.abs(object - newObject) < errMax) {
                object = newObject;
                break;
            }
            object = newObject;
        }

    }
    
    
    
    public double LF(){
    	this.WH = (W.transpose().times(H));
    	double E;
    	double[] rs = new double[944];
    	for(int c=1;c<WH.numRows();c++){
    		double sum=0;
    		for(int r=1;r<WH.numCols();r++){
    			if(V.get(c,r)>0){
    				sum+=WH.get(c,r);
    			}  			
    		}
    		rs[c]=sum;
    	}    	
    	double P1 = 0;    	    		
    	for(int i=1; i<n; i++){
    		for(int j=1; j<m; j++){
    			if(V.get(i, j)>0) {
    				P1+=V.get(i,j)*(Math.log(WH.get(i,j))-Math.log(rs[i]));
    			}
    		}
    	}
    	P1 = -P1;
    		
    	double Wt = 0;	
    	for(int i=1;i<W.numCols();i++){
    		for(int j=1;j<W.numRows();j++){
    			Wt+=Math.pow(W.get(i,j),2);
    		}
    	}
    	
    	double Ht = 0;	
    	for(int i=1;i<H.numCols();i++){
    		for(int j=1;j<H.numRows();j++){
    			Ht+=Math.pow(H.get(i,j),2);
    		}
    	}
    	double P2=(lambda/2)*(Wt+Ht);
    	
    	double P3=0;
    	double row[]=new double[944];
    	for(int i=1;i<WH.numRows();i++){
    		double sum=0;
    		for(int j=1;j<WH.numCols();j++){
    			sum+=WH.get(i,j);
    		}
    		row[i]=sum;
    	}
    	for(int i=1;i<WH.numRows();i++){
    		P3+=row[i];
    	}
    	E=P1+P2+P3;
    	return E;
    }

    public double GDA(){
    	
    	GaussJordanInverter M = new GaussJordanInverter(WH1);
    	
    	M.inverse();
    	for(int ui=1;ui<W.numRows(); ui++){
    		for(int uj=1;uj<W.numCols();uj++){
    			double P1=0;
    	    	for(int vi=1; vi<n; vi++){
    	    		for(int vj=1;vj<m;vj++){
    	    			if(V.get(vi,vj)>0){
    	    				P1+=V.get(vi,vj)*((H.times(WH)-(1/Wt.get(i,j)));
    	    			}
    	    		}
    	    	}
    	    	P1=-P1;
    		}
    	}
    	
    	double P2=lambda*W.determinant();
    	double P3=0;
    	for(int i=1; i<V.numRows(); i++){
    		for(int j=1;j<V.numCols();j++){
    			P3+=H.get(i,j);
    		}
    	}
    	return P1+P2+P3;
    }
    
    public double GDB(){
    	double P1=0;
    	double sum=0;
    	for(int i=1; i<H.numRows(); i++){
    		for(int j=1;j<H.numCols();j++){
    			sum+=H.get(i,j);
    		}
    	}
    	for(int i=1; i<W.numRows(); i++){
    		for(int j=1;j<W.numCols();j++){
    			if(W.get(i,j)>0){
    				P1+=V.get(i,j)*((W.get(i,j)/WH.get(i,j))-(W.get(i,j)/(Wt.get(i,j)*sum)));
    			}
    		}	
    	}
    	P1=-P1;
    	double P2=lambda*H.determinant();
    	double P3=0;
    	for(int i=1; i<W.numRows(); i++){
    		for(int j=1;j<W.numCols();j++){
    			P3+=W.get(i,j);
    		}
    	}
    	P3=-P3;
    	return P1+P2+P3;
    }
    
    public double GDC(){
    	double P3=0;
    	for(int i=1; i<H.numRows(); i++){
    		for(int j=1;j<H.numCols();j++){
    			P3+=H.get(i,j);
    		}
    	}
    	return 1-P3;
    }

    //@Override
    public double calObject(Matrix V, Matrix W, Matrix H) {
        Matrix WH = (W.times(H));
        double err = 0;
        for(int i=0; i<n; i++) for(int j=0; j<m; j++)
            if(V.get(i, j) != 0) {
                err += Math.pow(V.get(i, j)-WH.get(i, j), 2);
                for(int k=0; k<r; k++) {
                    err += (beta/2)*(Math.pow(W.get(i,k), 2)+Math.pow(H.get(k,j), 2));
                }
            }
        return err;
    }

    //@Override
    public Matrix getW() {
        return W;
    }

    //@Override
    public Matrix getH() {
        return H;
    }

    //@Override
    public int getStep() {
        return stepNum;
    }

    //@Override
    public double getObjectFunctionValue(){
        return object;
    }
}
