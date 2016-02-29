package ListMF.LPMF;



import org.apache.mahout.math.Matrix;

public interface MF {
	Matrix getW();

    Matrix getH();

    void solve();

    int getStep();

    double calObject(Matrix V, Matrix W, Matrix H);

    double getObjectFunctionValue();
}
