package MAYGEN;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class NewClass {
    private final MAYGEN maygen;

    public NewClass(MAYGEN maygen) {
        this.maygen = maygen;
    }

    public void run(int[] degree) {
        boolean[] learningFromConnectivity = new boolean[]{false};
        int[] nonCanonicalIndices = new int[2];
        ArrayList<ArrayList<Permutation>> formerPermutations = new ArrayList<ArrayList<Permutation>>();
        int[] hydrogens = maygen.setHydrogens(degree);
        int[] newPartition = maygen.getPartition(degree);
        if (maygen.writeSDF)
            maygen.symbolArrayCopy = Arrays.copyOf(maygen.symbolArray, maygen.symbolArray.length);
        final int[] initialPartition;
        if (maygen.writeSDF) {
            initialPartition = maygen.sortWithPartition(newPartition, degree, maygen.symbolArrayCopy, hydrogens);
        } else {
            initialPartition = maygen.sortWithPartition(newPartition, degree, maygen.symbolArray, hydrogens);
        }
        maygen.partSize.set(0);
        int[] connectivityIndices = new int[2];
        learningFromConnectivity[0] = false;
        maygen.learningFromCanonicalTest.set(false);
        int[][] partitionList = new int[maygen.size + 1][1];
        try {
            maygen.partSize.set(maygen.partSize.get() + (maygen.findZeros(initialPartition) - 1));
            maygen.setYZValues(initialPartition);
            partitionList[0] = initialPartition;
            maygen.generate(degree, initialPartition, partitionList, connectivityIndices, learningFromConnectivity,
                    nonCanonicalIndices, formerPermutations, hydrogens);
        } catch (IOException | CloneNotSupportedException | CDKException ignored) {
        }
    }
}
