/**
 * MIT License
 *
 * <p>Copyright (c) 2021 Mehmet Aziz Yirik
 *
 * <p>Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 * and associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * <p>The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * <p>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * This is the main cass of MAYGEN project for molecular structure generation for a given input
 * molecular formula.
 *
 * @author Mehmet Aziz Yirik
 */
package MAYGEN;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;

public class Generation {
    private final MAYGEN maygen;

    public Generation(MAYGEN maygen) {
        this.maygen = maygen;
    }

    public void run(int[] degree) {
        boolean[] learningFromConnectivity = new boolean[] {false};
        int[] nonCanonicalIndices = new int[2];
        ArrayList<ArrayList<Permutation>> formerPermutations =
                new ArrayList<ArrayList<Permutation>>();
        int[] hydrogens = maygen.setHydrogens(degree);
        int[] newPartition = maygen.getPartition(degree);
        if (maygen.writeSDF)
            maygen.symbolArrayCopy = Arrays.copyOf(maygen.symbolArray, maygen.symbolArray.length);
        final int[] initialPartition;
        if (maygen.writeSDF) {
            initialPartition =
                    maygen.sortWithPartition(
                            newPartition, degree, maygen.symbolArrayCopy, hydrogens);
        } else {
            initialPartition =
                    maygen.sortWithPartition(newPartition, degree, maygen.symbolArray, hydrogens);
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
            maygen.generate(
                    degree,
                    initialPartition,
                    partitionList,
                    connectivityIndices,
                    learningFromConnectivity,
                    nonCanonicalIndices,
                    formerPermutations,
                    hydrogens);
        } catch (IOException | CloneNotSupportedException | CDKException ignored) {
        }
    }
}
