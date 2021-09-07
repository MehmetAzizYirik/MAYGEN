/*
 MIT License

 Copyright (c) 2021 Mehmet Aziz Yirik

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
 This is the main cass of MAYGEN project for molecular structure generation for a given input
 molecular formula.

 @author Mehmet Aziz Yirik
*/
package maygen;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.interfaces.IAtomContainer;

public class Generation {
    private final MAYGEN maygen;

    public Generation(MAYGEN maygen) {
        this.maygen = maygen;
    }

    public int[] sortWithPartition(
            int[] partitionList, int[] degrees, String[] symbols, int[] hydrogens) {
        int[] partition = maygen.buildArray(partitionList);
        int size = partition.length;
        for (int n = 0; n < size; n++) {
            for (int m = 0; m < (size - 1) - n; m++) {
                if ((partition[m] > partition[m + 1])) {
                    maygen.swap(partition, m, (m + 1));
                    maygen.swap(degrees, m, (m + 1));
                    maygen.swap(hydrogens, m, (m + 1));
                    maygen.swap(symbols, m, (m + 1));
                }
            }
        }
        maygen.reOrder(partition, degrees, symbols, hydrogens);
        return maygen.initialPartition(partition);
    }

    public void run(int[] degree) throws CloneNotSupportedException, CDKException {
        IAtomContainer atomContainer = maygen.builder.newInstance(IAtomContainer.class);
        int[] partSize = new int[] {0};
        int[] r = new int[] {0};
        int[] y = new int[] {0};
        int[] z = new int[] {0};
        String[] symbolArrayCopy = new String[maygen.symbolArray.length];
        int[][] ys = new int[][] {new int[0]};
        int[][] zs = new int[][] {new int[0]};
        boolean[] learningFromCanonicalTest = new boolean[] {false};
        boolean[] learningFromConnectivity = new boolean[] {false};
        int[] nonCanonicalIndices = new int[2];
        ArrayList<ArrayList<Permutation>> formerPermutations = new ArrayList<>();
        int[] hydrogens = maygen.setHydrogens(degree);
        int[] newPartition = maygen.getPartition(degree);
        if (maygen.writeSDF || maygen.printSDF || maygen.writeSMILES || maygen.printSMILES)
            symbolArrayCopy = Arrays.copyOf(maygen.symbolArray, maygen.symbolArray.length);
        final int[] initialPartition;
        if (maygen.writeSDF || maygen.printSDF || maygen.writeSMILES || maygen.printSMILES) {
            initialPartition = sortWithPartition(newPartition, degree, symbolArrayCopy, hydrogens);
        } else {
            initialPartition =
                    sortWithPartition(newPartition, degree, maygen.symbolArray, hydrogens);
        }
        if (maygen.writeSDF || maygen.printSDF || maygen.writeSMILES || maygen.printSMILES)
            atomContainer = maygen.initAC(atomContainer, symbolArrayCopy);
        int[] connectivityIndices = new int[2];
        int[][] partitionList = new int[maygen.size + 1][1];
        try {

            partSize[0] = partSize[0] + (maygen.findZeros(initialPartition) - 1);
            maygen.setYZValues(initialPartition, ys, zs);
            partitionList[0] = initialPartition;
            maygen.generate(
                    atomContainer,
                    degree,
                    initialPartition,
                    partitionList,
                    connectivityIndices,
                    learningFromConnectivity,
                    nonCanonicalIndices,
                    formerPermutations,
                    hydrogens,
                    partSize,
                    r,
                    y,
                    z,
                    ys,
                    zs,
                    learningFromCanonicalTest);
        } catch (IOException ignored) {
        }
    }
}
