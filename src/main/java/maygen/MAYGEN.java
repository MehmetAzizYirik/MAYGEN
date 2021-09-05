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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class MAYGEN {
    public static final String NUMBERS_FROM_0_TO_9 = "(?=[0-9])";
    public static final String LETTERS_FROM_A_TO_Z = "(?=[A-Z])";
    public int size = 0;
    public int total = 0;
    public boolean tsvoutput = false;
    public boolean writeSDF = false;
    public boolean coordinates = false;
    public SDFWriter sdfOut;
    public Writer smilesOut;
    public boolean writeSMILES = false;
    public boolean printSDF = false;
    public boolean printSMILES = false;
    public boolean multiThread = false;
    public int hIndex = 0;
    public AtomicInteger count = new AtomicInteger();
    public AtomicInteger indexSmiles = new AtomicInteger();
    public AtomicInteger indexSdf = new AtomicInteger();
    public int matrixSize = 0;
    public boolean verbose = false;
    public String formula;
    public String filedir;
    public ArrayList<String> symbols = new ArrayList<>();
    public int[] occurrences;
    public Map<String, Integer> valences;
    public int[] nodeLabels;
    public int graphSize;
    public List<int[]> oxygenSulfur = new ArrayList<>();
    public int[] firstDegrees;
    public int totalHydrogen = 0;
    public ArrayList<String> firstSymbols = new ArrayList<>();
    public int[] firstOccurrences;
    public boolean callHydrogenDistributor = false;
    public boolean justH = false;
    public boolean noHydrogen = false;
    public int sizePart = 0;
    public boolean singleAtom = true;
    public boolean onlyDegree2 = true;
    public boolean OnSm = true;
    public int oxygen = 0;
    public int sulfur = 0;
    public String[] symbolArray;
    public IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
    public SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique);
    public IAtomContainer atomContainer = builder.newInstance(IAtomContainer.class);

    {
        // The atom valences from CDK.
        valences = new HashMap<>();

        valences.put("C", 4);
        valences.put("N", 3);
        valences.put("O", 2);
        valences.put("S", 2);
        valences.put("P", 3);
        valences.put("F", 1);
        valences.put("I", 1);
        valences.put("Cl", 1);
        valences.put("Br", 1);
        valences.put("H", 1);
    }

    /* Basic functions */

    /**
     * Permuting two entries of an Integer array.
     *
     * @param array Integer[] array
     * @param i int first index
     * @param j int second index
     * @return int[]
     */
    public int[] permuteArray(int[] array, int i, int j) {
        int temp = 0;
        temp = array[i];
        array[i] = array[j];
        array[j] = temp;
        return array;
    }

    /**
     * Summing entries of an array.
     *
     * @param array int[]
     * @return int sum
     */
    public int sum(int[] array) {
        int sum = 0;
        for (int i = 0; i < array.length; i++) {
            sum = sum + array[i];
        }
        return sum;
    }

    /**
     * Summing entries of a list until a given index.
     *
     * @param list the int array
     * @param index the index
     * @return int sum
     */
    public int sum(int[] list, int index) {
        int sum = 0;
        for (int i = 0; i <= index; i++) {
            sum += list[i];
        }
        return sum;
    }

    /**
     * Getting the number of atoms' occurrences.
     *
     * @param info String[] atom info
     * @return int
     */
    public int atomOccurrence(String[] info) {
        return info.length > 1 ? Integer.parseInt(info[1]) : 1;
    }

    /**
     * Performing the permutation action on an int array.
     *
     * @param array int[] array
     * @param permutation Permutation permutation
     * @return int[]
     */
    public int[] actArray(int[] array, Permutation permutation) {
        int permLength = permutation.size();
        int newIndex;
        int arrayLength = array.length;
        int[] modified = new int[arrayLength];
        for (int i = 0; i < permLength; i++) {
            newIndex = permutation.get(i);
            modified[newIndex] = array[i];
        }
        return modified;
    }

    /**
     * Values for an id permutation for a given size
     *
     * @param size int permutation size
     * @return int[]
     */
    public int[] idValues(int size) {
        int[] id = new int[size];
        for (int i = 0; i < size; i++) {
            id[i] = i;
        }
        return id;
    }

    /**
     * Builds id permutation.
     *
     * @param size int Permutation size
     * @return Permutation
     */
    public Permutation idPermutation(int size) {
        return new Permutation(size);
    }

    /**
     * The initializer function, reading the formula to set the degrees, partition and file
     * directory variables.
     *
     * @param symbols the symbols
     */
    public void sortAscending(ArrayList<String> symbols) {
        HashMap<String, Integer> inputs = new HashMap<>();
        for (int i = 0; i < symbols.size(); i++) {
            if (inputs.containsKey(symbols.get(i))) {
                Integer count = inputs.get(symbols.get(i)) + 1;
                inputs.put(symbols.get(i), count);
            } else {
                inputs.put(symbols.get(i), 1);
            }
        }

        Set<Entry<String, Integer>> set = inputs.entrySet();
        sort(symbols, set);
    }

    public void sort(ArrayList<String> symbols, Set<Entry<String, Integer>> set) {
        int index = 0;
        int value;
        ArrayList<Entry<String, Integer>> list = new ArrayList<>(set);
        list.sort(
                (value1, value2) -> {
                    int comparison = (value2.getValue()).compareTo(value1.getValue());
                    if (comparison == 0) {
                        return valences.get(value1.getKey())
                                .compareTo(valences.get(value2.getKey()));
                    } else {
                        return value1.getValue().compareTo(value2.getValue());
                    }
                });

        for (Entry<String, Integer> entry : list) {
            value = entry.getValue();
            for (int i = 0; i < value; i++) {
                symbols.set(index + i, entry.getKey());
            }
            index += value;
        }
    }

    public void singleAtomCheck(String[] atoms) {
        String[] info = atoms[0].split(NUMBERS_FROM_0_TO_9, 2);
        String symbol = info[0];
        if (atoms.length == 1) {
            if (symbol.equals("H")) {
                singleAtom = false;
            } else {
                if (atomOccurrence(info) > 1) {
                    singleAtom = false;
                }
            }
        } else {
            for (String atom : atoms) {
                info = atom.split(NUMBERS_FROM_0_TO_9, 2);
                symbol = info[0];
                if (!symbol.equals("H")) {
                    if (atomOccurrence(info) > 1) {
                        singleAtom = false;
                        break;
                    }
                }
            }
        }
    }

    public void checkOxygenSulfur(String[] atoms) {
        String[] info;
        String symbol;
        for (String atom : atoms) {
            info = atom.split(NUMBERS_FROM_0_TO_9, 2);
            symbol = info[0];
            if (valences.get(symbol) != 2) {
                onlyDegree2 = false;
                OnSm = false;
                break;
            } else {
                if (symbol.equals("S")) {
                    sulfur = atomOccurrence(info);
                } else if (symbol.equals("O")) {
                    oxygen = atomOccurrence(info);
                }
            }
        }
        if (onlyDegree2) {
            matrixSize = sulfur + oxygen;
            hIndex = matrixSize;
        }
    }

    public void getSingleAtomVariables() {
        String[] atoms = formula.split(LETTERS_FROM_A_TO_Z);
        ArrayList<String> symbolList = new ArrayList<>();
        String[] info;
        int hydrogens = 0;
        String symbol;
        hIndex = 1;
        for (String atom : atoms) {
            info = atom.split(NUMBERS_FROM_0_TO_9, 2);
            symbol = info[0];
            if (symbol.equals("H")) {
                hydrogens = atomOccurrence(info);
            } else {
                symbolList.add(symbol);
            }
        }
        matrixSize = hydrogens + 1;
        for (int i = 0; i < hydrogens; i++) {
            symbolList.add("H");
        }
        setSymbols(symbolList);
    }

    public void getSymbolOccurrences() {
        String[] atoms = formula.split(LETTERS_FROM_A_TO_Z);
        ArrayList<String> symbolList = new ArrayList<>();
        String[] info;
        int occur;
        int hydrogens = 0;
        String symbol;
        for (String atom : atoms) {
            info = atom.split(NUMBERS_FROM_0_TO_9, 2);
            symbol = info[0];
            if (!symbol.equals("H")) {
                occur = atomOccurrence(info);
                sizePart++;
                for (int i = 0; i < occur; i++) {
                    symbolList.add(symbol);
                    hIndex++;
                }
            } else {
                hydrogens = atomOccurrence(info);
            }
        }
        sortAscending(symbolList);
        for (int i = 0; i < hydrogens; i++) {
            symbolList.add("H");
        }
        firstOccurrences = getPartition(symbolList);
        matrixSize = sum(firstOccurrences);
        setSymbols(symbolList);
        occurrences = getPartition(symbolList);
        if (hydrogens != 0) {
            totalHydrogen += hydrogens;
            if (hIndex == 1) {
                callHydrogenDistributor = false;
            } else if (hIndex == 0) {
                justH = true;
                callHydrogenDistributor = false;
                hIndex = hydrogens;
                matrixSize = hIndex;
            } else {
                callHydrogenDistributor = true;
            }
        } else {
            callHydrogenDistributor = false;
            noHydrogen = true;
        }
    }

    public int[] nextCount(int index, int i, int size, ArrayList<String> symbols, int[] partition) {
        int localCount = 1;
        if (i == (size - 1)) {
            partition[index] = 1;
            index++;
        } else {
            for (int j = i + 1; j < size; j++) {
                if (symbols.get(i).equals(symbols.get(j))) {
                    localCount++;
                    if (j == (size - 1)) {
                        partition[index] = localCount;
                        index++;
                        break;
                    }
                } else {
                    partition[index] = localCount;
                    index++;
                    break;
                }
            }
        }
        return new int[] {localCount, index};
    }

    public int[] getPartition(ArrayList<String> symbols) {
        int i = 0;
        int[] partition = new int[sizePart + 1];
        int localSize = symbols.size();
        int next = 0;
        int index = 0;
        int[] result;
        while (i < localSize) {
            result = nextCount(index, i, localSize, symbols, partition);
            next = (i + result[0]);
            index = result[1];
            if (next == localSize) {
                break;
            } else {
                i = next;
            }
        }
        return partition;
    }
    /**
     * Setting the firstSymbols and symbols global variables for the initial sorted list of symbols.
     *
     * @param symbolList sorted list of atom symbols
     */
    public void setSymbols(ArrayList<String> symbolList) {
        symbolArray = new String[matrixSize];
        int index = 0;
        for (String symbol : symbolList) {
            symbolArray[index] = symbol;
            index++;
            if (!firstSymbols.contains(symbol)) {
                firstSymbols.add(symbol);
            }
        }
    }

    public String normalizeFormula(String formula) {
        String[] from = {"cl", "CL", "c", "n", "o", "s", "p", "f", "i", "br", "BR", "h"};
        String[] to = {"Cl", "Cl", "C", "N", "O", "S", "P", "F", "I", "Br", "Br", "H"};
        return StringUtils.replaceEach(formula, from, to);
    }

    public String[] validateFormula(String formula) {
        String[] from = {"Cl", "C", "N", "O", "S", "P", "F", "I", "Br", "H"};
        String[] to = {"", "", "", "", "", "", "", "", "", ""};
        String result = StringUtils.replaceEach(formula.replaceAll("[0-9]", ""), from, to);
        return result.isEmpty() ? new String[0] : result.split("");
    }

    /**
     * Checking whether a molecular formula can represent a graph or not.
     *
     * <p>For a graph with n nodes, the sum of all its node degrees should be equal or bigger than
     * 2*(n-1). Thus, the minimum number of nodes.
     *
     * @param formula String molecular formula
     * @return boolean
     */
    public boolean canBuildIsomer(String formula) {
        String[] atoms = normalizeFormula(formula).split(LETTERS_FROM_A_TO_Z);
        String[] info;
        String symbol;
        int occur;
        int valence;
        int localSize = 0;
        int sum = 0;
        for (String atom : atoms) {
            info = atom.split(NUMBERS_FROM_0_TO_9, 2);
            symbol = info[0];

            valence = valences.get(symbol);
            occur = atomOccurrence(info);
            localSize += occur;
            sum += (valence * occur);
        }
        total = localSize;
        return sum % 2 == 0 && sum >= 2 * (localSize - 1);
    }

    /** Initial degree arrays are set based on the molecular formula. */
    public void initialDegrees() {
        firstDegrees = new int[matrixSize];
        int index = 0;
        String symbol;
        int length = firstSymbols.size();
        for (int i = 0; i < length; i++) {
            symbol = firstSymbols.get(i);
            for (int j = 0; j < firstOccurrences[i]; j++) {
                firstDegrees[index] = valences.get(symbol);
                index++;
            }
        }
    }

    /**
     * Checks two int[] arrays are equal with respect to an atom partition.
     *
     * @param array1 int[] first array
     * @param array2 int[] second array
     * @param partition int[] atom partition
     * @return boolean
     */
    public boolean equalSetCheck(int[] array1, int[] array2, int[] partition) {
        int[] temp = cloneArray(array2);
        temp = descendingSortWithPartition(temp, partition);
        return equalSetCheck2(partition, array1, temp);
    }

    /**
     * Getting a part of a int array specified by two entry indices
     *
     * @param array int[] array
     * @param begin int beginning index
     * @param end int ending index
     * @return Integer[]
     */
    public int[] getBlocks(int[] array, int begin, int end) {
        return Arrays.copyOfRange(array, begin, end);
    }

    /**
     * Checks two int[] arrays are equal with respect to an atom partition.
     *
     * @param partition int[] atom partition
     * @param array1 int[] array
     * @param array2 int[] array
     * @return boolean
     */
    public boolean equalSetCheck2(int[] partition, int[] array1, int[] array2) {
        boolean check = true;
        int i = 0;
        int limit = findZeros(partition);
        if (partition[size - 1] != 0) {
            for (int d = 0; d < size; d++) {
                if (array1[d] != array2[d]) {
                    check = false;
                    break;
                }
            }
        } else {
            int value;
            for (int s = 0; s < limit; s++) {
                value = partition[s];
                if (compareIndexwise(array1, array2, i, (value + i))) {
                    i = i + value;
                } else {
                    check = false;
                    break;
                }
            }
        }
        return check;
    }

    /**
     * Comparing two int arrays are equal or not for given range of entries.
     *
     * @param array int[] array
     * @param array2 int[] array
     * @param index1 int beginning index
     * @param index2 int last index
     * @return boolean
     */
    public boolean compareIndexwise(int[] array, int[] array2, int index1, int index2) {
        boolean check = true;
        for (int i = index1; i < index2; i++) {
            if (array[i] != array2[i]) {
                check = false;
                break;
            }
        }
        return check;
    }

    /**
     * Comparing two arrays are equal. The second row's index and entries are permuted based on
     * cycle transposition and given permutation.
     *
     * @param index int row index
     * @param A int[][] adjacency matrix
     * @param cycleTransposition Permutation cycle transposition
     * @param permutation Permutation permutation
     * @return boolean
     */
    public boolean equalRowsCheck(
            int index, int[][] A, Permutation cycleTransposition, Permutation permutation) {
        int[] canonical = A[index];
        int[] original;
        int newIndex = findIndex(index, cycleTransposition);
        Permutation pm = permutation.multiply(cycleTransposition);
        original = cloneArray(A[newIndex]);
        original = actArray(original, pm);
        return Arrays.equals(canonical, original);
    }

    /**
     * Sorting entries of a subarray specified by indices.
     *
     * @param array int[] array
     * @param index0 int beginning index
     * @param index1 int last index
     * @return int[]
     */
    public int[] descendingSort(int[] array, int index0, int index1) {
        int temp = 0;
        for (int i = index0; i < index1; i++) {
            for (int j = i + 1; j < index1; j++) {
                if (array[i] < array[j]) {
                    temp = array[i];
                    array[i] = array[j];
                    array[j] = temp;
                }
            }
        }
        return array;
    }

    /**
     * Sorting entries of a int array for a given atom partition.
     *
     * @param array int[] array
     * @param partition int[] atom partition
     * @return int[]
     */
    public int[] descendingSortWithPartition(int[] array, int[] partition) {
        int i = 0;
        int p = 0;
        int limit = findZeros(partition);
        for (int i1 = 0; i1 < limit; i1++) {
            p = partition[i1];
            array = descendingSort(array, i, i + p);
            i = i + p;
        }
        return array;
    }

    /**
     * Checks two arrays in descending order. The second row is sorted in descending order just to
     * check whether there is a possible permutations making the first row non-maximal.
     *
     * @param index int row index
     * @param firstRow int[] array
     * @param check int[] array
     * @param partition int[] atom partition
     * @return boolean
     */
    public boolean biggerCheck(int index, int[] firstRow, int[] check, int[] partition) {
        int[] sorted = cloneArray(check);
        sorted = descendingSortWithPartition(sorted, partition);
        return descendingOrderUpperMatrixCheck(index, partition, firstRow, sorted);
    }

    /**
     * Checks whether there is a permutation making the row bigger in descending order.
     *
     * @param index int row index
     * @param A int[][] adjacency matrix
     * @param permutation Permutation permutation
     * @param partition int[] atom partition
     * @return boolean
     */
    public boolean setBiggest(int index, int[][] A, Permutation permutation, int[] partition) {
        int[] check = row2compare(index, A, permutation);
        return biggerCheck(index, A[index], check, partition);
    }

    /**
     * Get indices from "learning from canonical test" method. Here, the entry makes the row
     * non-canonical is detected. Its indices are set to nonCanonicalIndices global variables.
     *
     * @param index int row index
     * @param A the adjacency matrix
     * @param cycles the list of cycle transpositions
     * @param partition the atom partition
     * @param nonCanonicalIndices the nonCanonicalIndices
     * @param learningFromCanonicalTest the learningFromCanonicalTest
     */
    public void getLernenIndices(
            int index,
            int[][] A,
            ArrayList<Permutation> cycles,
            int[] partition,
            int[] nonCanonicalIndices,
            boolean[] learningFromCanonicalTest) {
        int[] check;
        for (Permutation cycle : cycles) {
            check = row2compare(index, A, cycle);
            if (!biggerCheck(index, A[index], check, partition)) {
                setLernenIndices(
                        index,
                        cycle,
                        A,
                        check,
                        partition,
                        nonCanonicalIndices,
                        learningFromCanonicalTest);
                break;
            }
        }
    }

    /**
     * Setting the nonCanonicalIndices global variable.
     *
     * @param rowIndex1 int first row index
     * @param cycle Permutation cycle transposition
     * @param A int[][] adjacency matrix
     * @param secondRow int[] second row
     * @param partition int[] atom partition
     * @param nonCanonicalIndices the nonCanonicalIndices
     * @param learningFromCanonicalTest the learningFromCanonicalTest
     */
    public void setLernenIndices(
            int rowIndex1,
            Permutation cycle,
            int[][] A,
            int[] secondRow,
            int[] partition,
            int[] nonCanonicalIndices,
            boolean[] learningFromCanonicalTest) {
        System.arraycopy(new int[2], 0, nonCanonicalIndices, 0, 2);
        learningFromCanonicalTest[0] = false;
        int rowIndex2 = cycle.get(rowIndex1);
        Permutation permutation = getNonCanonicalMakerPermutation(secondRow, cycle, partition);
        learningFromCanonicalTest[0] = true;
        System.arraycopy(
                upperIndex(rowIndex1, rowIndex2, A, permutation), 0, nonCanonicalIndices, 0, 2);
    }

    /**
     * Calculating the permutation, permuting the second row and making first row non maximal.
     *
     * @param array the array
     * @param cycle the cycle
     * @param partition the partition
     * @return the Permutation
     */
    public Permutation getNonCanonicalMakerPermutation(
            int[] array, Permutation cycle, int[] partition) {
        int[] sorted = cloneArray(array);
        sorted = descendingSortWithPartition(sorted, partition);
        Permutation permutation = getCanonicalPermutation(sorted, array, partition);
        return permutation.multiply(cycle);
    }

    /**
     * For a row given by index, checking whether it is in maximal form or not. If not, the
     * nonCanonicalIndices is set.
     *
     * @param array int[] array
     * @return boolean
     */
    public boolean zero(int[] array) {
        boolean check = false;
        for (int i = 0; i < size; i++) {
            if (array[i] == 0) {
                check = true;
                break;
            }
        }
        return check;
    }

    public boolean rowDescendingTest(
            int index,
            int[][] A,
            int[] partition,
            int[] nonCanonicalIndices,
            boolean[] learningFromCanonicalTest) {
        boolean check = true;
        if (zero(partition)) {
            if (!descendingOrderCheck(partition, A[index])) {
                check = false;
                int[] array = cloneArray(A[index]);
                array = descendingSortWithPartition(array, partition);
                Permutation canonicalPermutation =
                        getCanonicalPermutation(array, A[index], partition);
                learningFromCanonicalTest[0] = true;
                System.arraycopy(
                        upperIndex(index, index, A, canonicalPermutation),
                        0,
                        nonCanonicalIndices,
                        0,
                        2);
            }
        }
        return check;
    }

    /**
     * By a given permutation, checking which entry is mapped to the index.
     *
     * @param permutation Permutation permutation
     * @param index int entry index in the row
     * @return int
     */
    public int getPermutedIndex(Permutation permutation, int index) {
        int out = 0;
        for (int i = 0; i < permutation.size(); i++) {
            if (permutation.get(i) == index) {
                out += i;
                break;
            }
        }
        return out;
    }

    /**
     * Looking for the upper limit where the original entry is smaller.
     *
     * @param index int row index
     * @param nextRowIndex int index of the row to compare
     * @param A int[][] adjacency matrix
     * @param permutation Permutation permutation from canonical test
     * @return int[]
     */
    public int[] limit(int index, int nextRowIndex, int[][] A, Permutation permutation) {
        int[] original = A[index];
        int[] permuted = A[nextRowIndex];
        int[] limit = new int[2];
        limit[0] = index;
        int newIndex, value, newValue;
        for (int i = index + 1; i < size; i++) {
            newIndex = getPermutedIndex(permutation, i);
            value = original[i];
            newValue = permuted[newIndex];
            if (value != newValue) {
                if (value < newValue) {
                    limit[1] = i;
                }
                break;
            }
        }
        return limit;
    }

    /**
     * Looking for the maximum index where the entry is not zero.
     *
     * @param index int row index
     * @param nextRowIndex int index of the row to compare
     * @param A int[][] adjacency matrix
     * @param permutation Permutation permutation from canonical test
     * @return int[]
     */
    public int[] lowerIndex(int index, int nextRowIndex, int[][] A, Permutation permutation) {
        int max = 0;
        int upperLimit = limit(index, nextRowIndex, A, permutation)[1];
        int[] permuted = A[nextRowIndex];
        int newIndex, newValue;
        for (int i = index + 1; i < upperLimit; i++) {
            newIndex = getPermutedIndex(permutation, i);
            newValue = permuted[newIndex];
            if (newValue > 0) {
                if (max < newIndex) {
                    max = newIndex;
                }
            }
        }
        return new int[] {nextRowIndex, max};
    }

    /**
     * We need to calculate upperIndex. First, we need our j index where the original row become
     * smaller, then calculating the lower index. Based on these two values and the value of j in
     * the permutation, we calculate our upper index.
     *
     * <p>This upper index is used for the 'learning from canonical test' method.
     *
     * @param index int row index
     * @param nextRowIndex int index of the row to compare
     * @param A int[][] adjacency matrix
     * @param permutation Permutation permutation from canonical test
     * @return int[]
     */
    public int[] upperIndex(int index, int nextRowIndex, int[][] A, Permutation permutation) {
        int[] limit = limit(index, nextRowIndex, A, permutation);
        int[] lowerLimit = lowerIndex(index, nextRowIndex, A, permutation);
        int[] upperLimit = new int[2];
        upperLimit[0] = nextRowIndex;
        upperLimit[1] = getPermutedIndex(permutation, limit[1]);
        int[] maximalIndices = getMaximumPair(upperLimit, getMaximumPair(limit, lowerLimit));
        maximalIndices = maximalIndexWithNonZeroEntry(A, maximalIndices);
        return getTranspose(maximalIndices);
    }

    /**
     * In case if the index' entry is zero, updating the index with next index with non-zero entry.
     *
     * @param A int[][] adjacency matrix
     * @param maximalIndices int[] maximal indices for canonical test
     * @return int[]
     */
    public int[] maximalIndexWithNonZeroEntry(int[][] A, int[] maximalIndices) {
        int rowIndex = maximalIndices[0];
        int columnIndex = maximalIndices[1];
        if ((columnIndex > rowIndex) && A[rowIndex][columnIndex] != 0) {
            return maximalIndices;
        } else {
            int[] output = new int[2];
            for (int i = columnIndex; i < size; i++) {
                if (A[rowIndex][i] > 0) {
                    output[0] = rowIndex;
                    output[1] = i;
                    break;
                }
            }
            return output;
        }
    }

    /**
     * For an index pair, getting its transpose.
     *
     * @param indices int[] indices
     * @return int[]
     */
    public int[] getTranspose(int[] indices) {
        int[] out = new int[2];
        if (indices[0] > indices[1]) {
            out[0] = indices[1];
            out[1] = indices[0];
            return out;
        } else {
            return indices;
        }
    }

    /**
     * Between two index pairs, getting the bigger indices.
     *
     * @param a int[] indices
     * @param b int[] indices
     * @return int[]
     */
    public int[] getMaximumPair(int[] a, int[] b) {
        if (a[0] > b[0]) {
            return a;
        } else if (b[0] > a[0]) {
            return b;
        } else {
            if (a[1] > b[1]) {
                return a;
            } else if (b[1] > a[1]) {
                return b;
            } else {
                return a;
            }
        }
    }

    /**
     * Comparing two arrays for specific range of entries, whether the first array is bigger than
     * the second one or not.
     *
     * @param array1 int[] first array
     * @param array2 int[] second array
     * @param index1 int beginning index
     * @param index2 int last index
     * @return boolean
     */
    public boolean compare(int[] array1, int[] array2, int index1, int index2) {
        boolean check = true;
        for (int i = index1; i < index2; i++) {
            if (array1[i] != array2[i]) {
                if (array1[i] < array2[i]) {
                    check = false;
                    break;
                } else if (array1[i] > array2[i]) {
                    check = true;
                    break;
                }
            }
        }
        return check;
    }

    /**
     * Checks the first row is bigger than the second row just in the upper matrix.
     *
     * @param index int row index
     * @param partition int[] atom partition
     * @param firstRow int[] first row
     * @param secondRow int[] second row
     * @return boolean
     */
    public boolean descendingOrderUpperMatrixCheck(
            int index, int[] partition, int[] firstRow, int[] secondRow) {
        boolean check = true;
        int i = index + 1;
        int p;
        int limit = findZeros(partition);
        for (int k = index + 1; k < limit; k++) {
            p = partition[k];
            if (descendingOrderCheck(firstRow, i, (i + p))) {
                if (!compareIndexwise(firstRow, secondRow, i, (i + p))) {
                    check = compare(firstRow, secondRow, i, (i + p));
                    break;
                }
                i = i + p;
            } else {
                check = false;
                break;
            }
        }
        return check;
    }

    /**
     * Checks subarray of specified range of entries, the array is descending order or not.
     *
     * @param array int[] array
     * @param f int first index
     * @param l int last index
     * @return boolean
     */
    public boolean descendingOrderCheck(int[] array, int f, int l) {
        boolean check = true;
        for (int i = f; i < l - 1; i++) {
            if (array[i] < array[i + 1]) {
                check = false;
                break;
            }
        }
        return check;
    }

    /**
     * Checks a int array is in descending order or not with respect to a given atom partition.
     *
     * @param partition int[] atom partition
     * @param array the int array
     * @return boolean
     */
    public boolean descendingOrderCheck(int[] partition, int[] array) {
        boolean check = true;
        int i = 0;
        int value = 0;
        int limit = findZeros(partition);
        for (int s = 0; s < limit; s++) {
            value = partition[s];
            if (!descendingOrderCheck(array, i, (value + i))) {
                check = false;
                break;
            } else {
                i = i + value;
            }
        }
        return check;
    }

    /* ******************************************************************** */

    /* Candidate Matrix Generation Functions */

    /**
     * L; upper triangular matrix like given in 3.2.1. For (i,j), after the index, giving the
     * maximum line capacity.
     *
     * @param degrees the degrees
     * @param max the max
     * @param L the L
     */
    public void upperTriangularL(int[] degrees, int[][][] max, int[][][] L) {
        L[0] = new int[hIndex][hIndex];
        if (hIndex == 2) {
            for (int i = 0; i < hIndex; i++) {
                for (int j = i + 1; j < hIndex; j++) {
                    L[0][i][j] = Math.min(degrees[i], Lsum(i, j, max));
                }
            }
        } else {
            for (int i = 0; i < hIndex; i++) {
                for (int j = i + 1; j < hIndex; j++) {
                    L[0][i][j] = Math.min(degrees[i], Lsum(i, j + 1, max));
                }
            }
        }
    }

    /**
     * C; upper triangular matrix like given in 3.2.1. For (i,j), after the index, giving the
     * maximum column capacity.
     *
     * @param degrees int[] valences
     * @param max the max
     * @param C the C
     */
    public void upperTriangularC(int[] degrees, int[][][] max, int[][][] C) {
        C[0] = new int[hIndex][hIndex];
        if (hIndex == 2) {
            for (int i = 0; i < hIndex; i++) {
                for (int j = i + 1; j < hIndex; j++) {
                    C[0][i][j] = Math.min(degrees[j], Csum(i, j, max));
                }
            }
        } else {
            for (int i = 0; i < hIndex; i++) {
                for (int j = i + 1; j < hIndex; j++) {
                    C[0][i][j] = Math.min(degrees[j], Csum(i + 1, j, max));
                }
            }
        }
    }

    /**
     * Summing ith rows entries starting from the jth column.
     *
     * @param i int row index
     * @param j int column index
     * @param max the max
     * @return the Lsum
     */
    public int Lsum(int i, int j, int[][][] max) {
        int sum = 0;
        for (int k = j; k < hIndex; k++) {
            sum = sum + max[0][i][k];
        }
        return sum;
    }

    /**
     * Summing ith column entries starting from the jth row.
     *
     * @param i int column index
     * @param j int row index
     * @param max the max
     * @return the Csum
     */
    public int Csum(int i, int j, int[][][] max) {
        int sum = 0;
        for (int k = i; k < hIndex; k++) {
            sum = sum + max[0][k][j];
        }
        return sum;
    }

    /**
     * Possible maximal edge multiplicity for the atom pair (i,j).
     *
     * @param degrees the degrees
     * @param max the max
     */
    public void maximalMatrix(int[] degrees, int[][][] max) {
        max[0] = new int[hIndex][hIndex];
        for (int i = 0; i < hIndex; i++) {
            for (int j = 0; j < hIndex; j++) {
                int di = degrees[i];
                int dj = degrees[j];
                if (i == j) {
                    max[0][i][j] = 0;
                } else {
                    if (di != dj) {
                        max[0][i][j] = Math.min(di, dj);
                    } else if (di == dj && i != j) {
                        if (justH) {
                            max[0][i][j] = (di);
                        } else {
                            if (hIndex == 2) {
                                max[0][i][j] = (di);
                            } else {
                                if (di != 1) {
                                    max[0][i][j] = (di - 1);
                                } else {
                                    max[0][i][j] = (di);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Initialization of global variables for the generate of structures for given degree list.
     *
     * @param ac the IAtomContainer
     * @param degreeList int[] valences
     * @param initialPartition the initial partition
     * @param partitionList the partitionList
     * @param connectivityIndices the connectivityIndices
     * @param learningFromConnectivity the learningFromConnectivity
     * @param nonCanonicalIndices the nonCanonicalIndices
     * @param formerPermutations the formerPermutations
     * @param hydrogens the hydrogens
     * @param partSize the partSize
     * @param r the r
     * @param y the y
     * @param z the z
     * @param ys the ys
     * @param zs the zs
     * @param learningFromCanonicalTest the learningFromCanonicalTest
     * @throws IOException in case of IOException
     * @throws CDKException in case of CDKException
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     */
    public void generate(
            IAtomContainer ac,
            int[] degreeList,
            int[] initialPartition,
            int[][] partitionList,
            int[] connectivityIndices,
            boolean[] learningFromConnectivity,
            int[] nonCanonicalIndices,
            ArrayList<ArrayList<Permutation>> formerPermutations,
            int[] hydrogens,
            int[] partSize,
            int[] r,
            int[] y,
            int[] z,
            int[][] ys,
            int[][] zs,
            boolean[] learningFromCanonicalTest)
            throws IOException, CloneNotSupportedException, CDKException {
        int[][] A = new int[matrixSize][matrixSize];
        int[] degrees = degreeList;
        boolean[] flag = new boolean[] {true};
        int[][][] max = new int[][][] {new int[0][0]};
        int[][][] L = new int[][][] {new int[0][0]};
        int[][][] C = new int[][][] {new int[0][0]};
        maximalMatrix(degrees, max);
        upperTriangularL(degrees, max, L);
        upperTriangularC(degrees, max, C);
        int[] indices = new int[2];
        indices[0] = 0;
        indices[1] = 1;
        boolean[] callForward = {true};
        r[0] = 0;
        y[0] = ys[0][r[0]];
        z[0] = zs[0][r[0]];
        while (flag[0]) {
            nextStep(
                    ac,
                    A,
                    indices,
                    degrees,
                    initialPartition,
                    partitionList,
                    callForward,
                    connectivityIndices,
                    learningFromConnectivity,
                    nonCanonicalIndices,
                    formerPermutations,
                    hydrogens,
                    partSize,
                    r,
                    y,
                    z,
                    max,
                    L,
                    C,
                    ys,
                    zs,
                    learningFromCanonicalTest,
                    flag);
            if (!flag[0]) {
                break;
            }
            if (learningFromConnectivity[0]) {
                indices = connectivityIndices;
                findR(indices, initialPartition, r);
                int value = indexYZ(initialPartition, r);
                y[0] = ys[0][value];
                clearFormers(false, y[0], partitionList, formerPermutations);
                learningFromConnectivity[0] = false;
                callForward[0] = false;
            } else {
                if (learningFromCanonicalTest[0]) {
                    indices = successor(nonCanonicalIndices, max[0].length);
                    findR(indices, initialPartition, r);
                    learningFromCanonicalTest[0] = false;
                    callForward[0] = false;
                }
            }
        }
    }

    /**
     * Calculation of the next index pair in a matrix.
     *
     * @param indices int[] index pair.
     * @param size int row length.
     * @return int[]
     */
    public int[] successor(int[] indices, int size) {
        int i0 = indices[0];
        int i1 = indices[1];
        if (i1 < (size - 1)) {
            indices[0] = i0;
            indices[1] = (i1 + 1);
        } else if (i0 < (size - 2) && i1 == (size - 1)) {
            indices[0] = (i0 + 1);
            indices[1] = (i0 + 2);
        }
        return indices;
    }

    /**
     * Calculation of the former index pair in a matrix.
     *
     * @param indices int[] index pair.
     * @param size int row length.
     * @return int[]
     */
    public int[] predecessor(int[] indices, int size) {
        int i0 = indices[0];
        int i1 = indices[1];
        if (i0 == i1 - 1) {
            indices[0] = (i0 - 1);
            indices[1] = (size - 1);
        } else {
            indices[0] = i0;
            indices[1] = (i1 - 1);
        }
        return indices;
    }

    /**
     * Calling
     *
     * @param ac the IAtomContainer
     * @param A the A matrix
     * @param indices the indices
     * @param degrees the degrees
     * @param initialPartition the initial partition
     * @param partitionList the partitionList
     * @param callForward the callForward
     * @param connectivityIndices the connectivityIndices
     * @param learningFromConnectivity the learningFromConnectivity
     * @param nonCanonicalIndices the nonCanonicalIndices
     * @param formerPermutations the formerPermutations
     * @param hydrogens the hydrogens
     * @param partSize the partSize
     * @param r the r
     * @param y the y
     * @param z the z
     * @param max the max
     * @param L the L
     * @param C the C
     * @param ys the ys
     * @param zs the zs
     * @param learningFromCanonicalTest the learningFromCanonicalTest
     * @param flag the flag
     * @throws IOException in case of IOException
     * @throws CDKException in case of CDKException
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     */
    public void nextStep(
            IAtomContainer ac,
            int[][] A,
            int[] indices,
            int[] degrees,
            int[] initialPartition,
            int[][] partitionList,
            boolean[] callForward,
            int[] connectivityIndices,
            boolean[] learningFromConnectivity,
            int[] nonCanonicalIndices,
            ArrayList<ArrayList<Permutation>> formerPermutations,
            int[] hydrogens,
            int[] partSize,
            int[] r,
            int[] y,
            int[] z,
            int[][][] max,
            int[][][] L,
            int[][][] C,
            int[][] ys,
            int[][] zs,
            boolean[] learningFromCanonicalTest,
            boolean[] flag)
            throws IOException, CloneNotSupportedException, CDKException {
        if (callForward[0]) {
            forward(
                    ac,
                    A,
                    indices,
                    degrees,
                    initialPartition,
                    partitionList,
                    callForward,
                    connectivityIndices,
                    learningFromConnectivity,
                    nonCanonicalIndices,
                    formerPermutations,
                    hydrogens,
                    partSize,
                    r,
                    y,
                    z,
                    max,
                    L,
                    C,
                    ys,
                    zs,
                    learningFromCanonicalTest);
        } else {
            backward(A, indices, degrees, initialPartition, callForward, r, max, L, C, flag);
        }
    }

    /**
     * After generating matrices, adding the hydrogen with respect to the pre-hydrogen distribution.
     *
     * @param A the adjacency matrix
     * @param index int beginning index for the hydrogen setting
     * @param hydrogens the hydrogens
     * @return the adjacency matrix
     */
    public int[][] addHydrogens(int[][] A, int index, int[] hydrogens) {
        if (singleAtom) {
            int hIndex = index;
            int hydrogen = valences.get(symbolArray[0]);
            for (int j = hIndex; j < hydrogen + hIndex; j++) {
                A[0][j] = 1;
                A[j][0] = 1;
            }
        } else if (callHydrogenDistributor) {
            int hIndex = index;
            int limit = 0;
            int hydrogen = 0;
            for (int i = 0; i < index; i++) {
                hydrogen = hydrogens[i];
                limit = hIndex + hydrogen;
                for (int j = hIndex; j < limit; j++) {
                    A[i][j] = 1;
                    A[j][i] = 1;
                }
                if (hydrogen != 0) {
                    hIndex = hIndex + hydrogen;
                }
            }
        }
        return A;
    }

    public void findR(int[] indices, int[] initialPartition, int[] r) {
        int block = 0;
        int index = 0;
        int part = 0;
        int rowIndex = indices[0];
        int limit = findZeros(initialPartition);
        for (int i = 0; i < limit; i++) {
            part = initialPartition[i];
            if (index <= rowIndex && rowIndex < (index + part)) {
                break;
            } else {
                block++;
            }
            index = index + part;
        }
        r[0] = block;
    }
    /**
     * The third line of the backward method in Grund 3.2.3. The criteria to decide which function
     * is needed: forward or backward.
     *
     * @param x the value in the adjacency matrix A[i][j]
     * @param lInverse lInverse value of indices {i,j}
     * @param l the l parameter
     * @return the criteria
     */
    public boolean backwardCriteria(int x, int lInverse, int l) {
        int newX = (x - 1);
        return lInverse - newX <= l;
    }

    /**
     * The third step in Grund 3.2.3.
     *
     * <p>Backward step in the algorithm.
     *
     * @param A the adjacency matrix
     * @param indices the indices
     * @param degrees the degrees
     * @param initialPartition the initial partition
     * @param callForward the callForward
     * @param r the r
     * @param max the max
     * @param L the L
     * @param C the C
     * @param flag the flag
     * @return the int[][]
     */
    public int[][] backward(
            int[][] A,
            int[] indices,
            int[] degrees,
            int[] initialPartition,
            boolean[] callForward,
            int[] r,
            int[][][] max,
            int[][][] L,
            int[][][] C,
            boolean[] flag) {
        int i = indices[0];
        int j = indices[1];

        if (i == 0 && j == 1) {
            flag[0] = false;
        } else {
            indices = predecessor(indices, max[0].length);
            // UPDATE
            findR(indices, initialPartition, r);
            i = indices[0];
            j = indices[1];
            int x = A[i][j];
            int l2 = LInverse(i, j, A, degrees);
            int c2 = CInverse(i, j, A, degrees);

            if (x > 0
                    && (backwardCriteria((x), l2, L[0][i][j])
                            && backwardCriteria((x), c2, C[0][i][j]))) {
                A[i][j] = (x - 1);
                A[j][i] = (x - 1);
                indices = successor(indices, max[0].length);
                // UPDATE
                findR(indices, initialPartition, r);
                callForward[0] = true;
            } else {
                callForward[0] = false;
            }
        }
        return A;
    }

    /**
     * Setting successor indices entry if there is a possible filling.
     *
     * @param ac the IAtomContainer
     * @param A the adjacency matrix
     * @param indices the entry indices
     * @param degrees the degrees
     * @param initialPartition the initial partition
     * @param partitionList the partitionList
     * @param callForward the callForward
     * @param connectivityIndices the connectivityIndices
     * @param learningFromConnectivity the learningFromConnectivity
     * @param nonCanonicalIndices the nonCanonicalIndices
     * @param formerPermutations the formerPermutations
     * @param hydrogens the hydrogens
     * @param partSize the partSize
     * @param r the r
     * @param y the y
     * @param z the z
     * @param max the max
     * @param L the L
     * @param C the C
     * @param ys the ys
     * @param zs the zs
     * @param learningFromCanonicalTest the learningFromCanonicalTest
     * @return int[][]
     * @throws IOException in case of IOException
     * @throws CDKException in case of CDKException
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     */
    public int[][] forward(
            IAtomContainer ac,
            int[][] A,
            int[] indices,
            int[] degrees,
            int[] initialPartition,
            int[][] partitionList,
            boolean[] callForward,
            int[] connectivityIndices,
            boolean[] learningFromConnectivity,
            int[] nonCanonicalIndices,
            ArrayList<ArrayList<Permutation>> formerPermutations,
            int[] hydrogens,
            int[] partSize,
            int[] r,
            int[] y,
            int[] z,
            int[][][] max,
            int[][][] L,
            int[][][] C,
            int[][] ys,
            int[][] zs,
            boolean[] learningFromCanonicalTest)
            throws IOException, CloneNotSupportedException, CDKException {
        int i = indices[0];
        int j = indices[1];
        int lInverse = LInverse(i, j, A, degrees);
        int cInverse = CInverse(i, j, A, degrees);
        int minimumValue = Math.min(max[0][i][j], Math.min(lInverse, cInverse));
        int maximumValue = maximalEntry(minimumValue, lInverse, L[0][i][j], cInverse, C[0][i][j]);
        callForward[0] = true;
        return forward(
                ac,
                lInverse,
                cInverse,
                maximumValue,
                i,
                j,
                A,
                indices,
                initialPartition,
                partitionList,
                callForward,
                connectivityIndices,
                learningFromConnectivity,
                nonCanonicalIndices,
                formerPermutations,
                hydrogens,
                partSize,
                r,
                y,
                z,
                max,
                L,
                C,
                ys,
                zs,
                learningFromCanonicalTest);
    }

    public int[][] forward(
            IAtomContainer ac,
            int lInverse,
            int cInverse,
            int maximalX,
            int i,
            int j,
            int[][] A,
            int[] indices,
            int[] initialPartition,
            int[][] partitionList,
            boolean[] callForward,
            int[] connectivityIndices,
            boolean[] learningFromConnectivity,
            int[] nonCanonicalIndices,
            ArrayList<ArrayList<Permutation>> formerPermutations,
            int[] hydrogens,
            int[] partSize,
            int[] r,
            int[] y,
            int[] z,
            int[][][] max,
            int[][][] L,
            int[][][] C,
            int[][] ys,
            int[][] zs,
            boolean[] learningFromCanonicalTest)
            throws IOException, CloneNotSupportedException, CDKException {
        if (((lInverse - maximalX) <= L[0][i][j]) && ((cInverse - maximalX) <= C[0][i][j])) {
            A[i][j] = maximalX;
            A[j][i] = maximalX;
            if (i == (max[0].length - 2) && j == (max[0].length - 1)) {
                if (canonicalTest(
                        A,
                        initialPartition,
                        partitionList,
                        nonCanonicalIndices,
                        formerPermutations,
                        partSize,
                        r,
                        y,
                        z,
                        ys,
                        zs,
                        learningFromCanonicalTest)) {
                    if (connectivityTest(A, connectivityIndices, learningFromConnectivity)) {
                        count.incrementAndGet();
                        if (writeSDF || printSDF) {
                            IAtomContainer ac2 =
                                    buildContainer4SDF(ac, addHydrogens(A, hIndex, hydrogens));
                            if (coordinates)
                                new StructureDiagramGenerator().generateCoordinates(ac2);
                            sdfOut.write(ac2);
                        }
                        if (writeSMILES || printSMILES) {
                            write2smiles(addHydrogens(A, hIndex, hydrogens));
                        }
                        callForward[0] = false;
                    } else {
                        callForward[0] = false;
                        learningFromConnectivity[0] = true;
                    }
                } else {
                    if (!learningFromCanonicalTest[0]) {
                        callForward[0] = false;
                    }
                }
            } else {
                int value = indexYZ(initialPartition, r);
                if (indices[0] == zs[0][value] && indices[1] == (max[0].length - 1)) {
                    callForward[0] =
                            canonicalTest(
                                    A,
                                    initialPartition,
                                    partitionList,
                                    nonCanonicalIndices,
                                    formerPermutations,
                                    partSize,
                                    r,
                                    y,
                                    z,
                                    ys,
                                    zs,
                                    learningFromCanonicalTest);
                    if (callForward[0]) {
                        indices = successor(indices, max[0].length);
                        // update
                        findR(indices, initialPartition, r);
                    } else {
                        callForward[0] = false;
                    }
                } else {
                    indices = successor(indices, max[0].length);
                    // update
                    findR(indices, initialPartition, r);
                    callForward[0] = true;
                }
            }
        } else {
            callForward[0] = false;
        }
        return A;
    }

    /**
     * Calculating the maximal entry for the indices.
     *
     * @param min int minimum of L, C amd maximal matrices for {i,j} indices.
     * @param lInverse int Linverse value of {i,j}
     * @param l int L value of {i,j}
     * @param cInverse int Cinverse value of {i,j}
     * @param c int C value of {i,j}
     * @return int max
     */
    public int maximalEntry(int min, int lInverse, int l, int cInverse, int c) {
        int max = 0;
        for (int v = min; v >= 0; v--) {
            if (((lInverse - v) <= l) && ((cInverse - v) <= c)) {
                max = max + v;
                break;
            }
        }
        return max;
    }

    /**
     * Calculating the sum of the entries in the ith row until the jth column.
     *
     * @param i int row index
     * @param j int column index
     * @param A the adjacency matrix
     * @param degrees the degrees
     * @return int
     */
    public int LInverse(int i, int j, int[][] A, int[] degrees) {
        int sum = 0;
        if (hIndex == 2) {
            for (int s = 0; s <= j; s++) {
                sum = sum + A[i][s];
            }
        } else {
            for (int s = 0; s < j; s++) {
                sum = sum + A[i][s];
            }
        }
        return degrees[i] - sum;
    }

    /**
     * Calculating the sum of the entries in the jth column until the ith row.
     *
     * @param i int row index
     * @param j int column index
     * @param A the adjacency matrix
     * @param degrees the degrees
     * @return int
     */
    public int CInverse(int i, int j, int[][] A, int[] degrees) {
        int sum = 0;
        if (hIndex == 2) {
            for (int s = 0; s <= i; s++) {
                sum = sum + A[s][j];
            }
        } else {
            for (int s = 0; s < i; s++) {
                sum = sum + A[s][j];
            }
        }
        return degrees[j] - sum;
    }

    /**
     * Based on the new degrees and the former partition, getting the new atom partition.
     *
     * @param degrees int[] new atom valences
     * @return int[]
     */
    public int[] getPartition(int[] degrees) {
        int[] newPartition = new int[degrees.length];
        int i = 0;
        int p;
        int length;
        if (justH || noHydrogen) {
            length = firstOccurrences.length;
        } else {
            length = firstOccurrences.length - 1;
        }
        int index = 0;
        for (int part = 0; part < length; part++) {
            p = firstOccurrences[part];
            int[] subArray = getBlocks(degrees, i, p + i);
            for (Integer item : getSubPartition(subArray)) {
                newPartition[index] = item;
                index++;
            }
            i = i + p;
        }
        return newPartition;
    }

    /**
     * Calculating the sub partitions for a given group of degrees.
     *
     * @param degrees int[] valences
     * @return int[]
     */
    public int[] getSubPartition(int[] degrees) {
        int i = 0;
        int size = degrees.length;
        int[] partition = new int[size];
        int next;
        int index = 0;
        int[] result;
        while (i < size) {
            result = nextCount(index, i, size, degrees, partition);
            index = result[1];
            next = (i + result[0]);
            if (next == size) {
                break;
            } else {
                i = next;
            }
        }
        return partition;
    }

    /**
     * Counting the occurrence of a value in a degree.
     *
     * @param index the index
     * @param i the i
     * @param size int number
     * @param degrees int[] valences
     * @param partition int[] partition
     * @return int
     */
    public int[] nextCount(int index, int i, int size, int[] degrees, int[] partition) {
        int count = 1;
        if (i == (size - 1)) {
            partition[index] = 1;
            index++;
        } else {
            for (int j = i + 1; j < size; j++) {
                if (degrees[i] == degrees[j]) {
                    count++;
                    if (j == (size - 1)) {
                        partition[index] = count;
                        index++;
                        break;
                    }
                } else {
                    partition[index] = count;
                    index++;
                    break;
                }
            }
        }
        return new int[] {count, index};
    }

    public boolean checkLengthTwoFormula(String[] atoms) {
        boolean check = true;

        if (atoms.length == 1) {
            String[] info = atoms[0].split(NUMBERS_FROM_0_TO_9, 2);
            if (info[1].equals("2")) {
                if (valences.get(info[0]) > 3) {
                    check = false;
                }
            }
        }
        return check;
    }
    /**
     * Main function to initialize the global variables and calling the generate function.
     *
     * @throws IOException in case of IOException
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     * @throws CDKException in case of CDKException
     */
    public void run() throws IOException, CDKException, CloneNotSupportedException {
        clearGlobals();
        formula = normalizeFormula(formula);
        String[] unsupportedSymbols = validateFormula(formula);
        if (unsupportedSymbols.length > 0) {
            if (verbose)
                System.out.println(
                        "The input formula consists user defined element types: "
                                + String.join(", ", unsupportedSymbols));
        } else {
            long startTime = System.nanoTime();
            if (verbose) System.out.println("MAYGEN is generating isomers of " + formula + "...");
            if (writeSDF || printSDF) {
                if (writeSDF) {
                    new File(filedir).mkdirs();
                    sdfOut =
                            new SDFWriter(
                                    new FileWriter(
                                            filedir + "/" + normalizeFormula(formula) + ".sdf"));
                } else {
                    sdfOut = new SDFWriter(new PrintWriter(System.out));
                }
            }
            if (writeSMILES || printSMILES) {
                if (writeSMILES) {
                    new File(filedir).mkdirs();
                    smilesOut = new FileWriter(filedir + "/" + normalizeFormula(formula) + ".smi");
                } else {
                    smilesOut = new PrintWriter(System.out);
                }
            }
            String[] atoms = formula.split(LETTERS_FROM_A_TO_Z);
            if (checkLengthTwoFormula(atoms)) {
                singleAtomCheck(atoms);
                if (singleAtom) {
                    getSingleAtomVariables();
                    initSingleAC();
                    writeSingleAtom(new int[] {});
                    displayStatistic(startTime);
                } else {
                    checkOxygenSulfur(atoms);
                    if (onlyDegree2) {
                        if (oxygen == 0 || sulfur == 0) {
                            degree2graph();
                        } else {
                            if (printSMILES) System.out.println("Formula is not supported");
                            distributeSulfurOxygen();
                        }
                        displayStatistic(startTime);
                    } else {
                        if (canBuildIsomer(formula)) {
                            getSymbolOccurrences();
                            initialDegrees();
                            structureGenerator();
                            displayStatistic(startTime);
                        } else {
                            if (verbose)
                                System.out.println(
                                        "The input formula, "
                                                + formula
                                                + ", does not represent any molecule.");
                        }
                    }
                }
            } else {
                if (verbose)
                    System.out.println(
                            "The input formula, " + formula + ", does not represent any molecule.");
            }
        }
    }

    public void displayStatistic(long startTime) throws IOException {
        if (writeSDF || printSDF) {
            sdfOut.close();
        }
        if (writeSMILES || printSMILES) {
            smilesOut.close();
        }
        long endTime = System.nanoTime() - startTime;
        double seconds = (double) endTime / 1000000000.0;
        DecimalFormat d = new DecimalFormat(".###");
        d.setDecimalFormatSymbols(DecimalFormatSymbols.getInstance(Locale.ENGLISH));

        if (verbose) {
            if (printSMILES || printSDF) System.out.println();
            System.out.println("The number of structures is: " + count);
            System.out.println("Time: " + d.format(seconds) + " seconds");
        }

        if (tsvoutput) {
            System.out.println(
                    formula
                            + "\t"
                            + count
                            + "\t"
                            + d.format(seconds)
                            + "\t"
                            + (multiThread ? size : 1));
        }
    }

    /**
     * If there are hydrogens in the formula, calling the hydrogenDistributor. This is the
     * pre-hydrogen distribution. Then, the new list of degrees is defined for each hydrogen
     * distribution.
     *
     * @return the list of int arrays
     */
    public ArrayList<int[]> distributeHydrogens() {
        ArrayList<int[]> degreeList = new ArrayList<>();
        if (!callHydrogenDistributor) {
            degreeList.add(firstDegrees);
        } else {
            List<int[]> distributions =
                    new HydrogenDistributor().run(firstOccurrences, firstDegrees);
            if (hIndex == 2) {
                for (int[] dist : distributions) {
                    int[] newDegree = new int[size];
                    for (int i = 0; i < size; i++) {
                        newDegree[i] = (firstDegrees[i] - dist[i]);
                    }
                    if (newDegree[0] == newDegree[1]) degreeList.add(newDegree);
                }
            } else {
                for (int[] dist : distributions) {
                    int[] newDegree = new int[size];
                    for (int i = 0; i < size; i++) {
                        newDegree[i] = (firstDegrees[i] - dist[i]);
                    }
                    degreeList.add(newDegree);
                }
            }
        }
        return degreeList;
    }

    /**
     * Setting the y and z values for each block. y is the beginning index and z is the last index
     * of a block in the adjacency matrix.
     *
     * @param initialPartition the initial partition
     * @param ys the ys
     * @param zs the zs
     */
    public void setYZValues(int[] initialPartition, int[][] ys, int[][] zs) {
        ys[0] = new int[size];
        zs[0] = new int[size];
        int limit = findZeros(initialPartition);
        int value = 0;
        int index = 0;
        int y = 0;
        int z = 0;
        for (int i = 0; i < limit; i++) {
            value = initialPartition[i];
            y = findY(i, initialPartition);
            z = findZ(i, initialPartition);
            for (int j = 0; j < value; j++) {
                ys[0][index] = y;
                zs[0][index] = z;
                index++;
            }
        }
    }

    /**
     * For a block index r, calculating its first row index.
     *
     * @param r int block index
     * @param initialPartition the initial partition
     * @return int
     */
    public int findY(int r, int[] initialPartition) {
        return sum(initialPartition, r - 1);
    }

    /**
     * For a block index r, calculating its last row index.
     *
     * @param r int block index
     * @param initialPartition the initial partition
     * @return int
     */
    public int findZ(int r, int[] initialPartition) {
        return (sum(initialPartition, r) - 1);
    }

    public void writeSingleAtom(int[] hydrogens)
            throws IOException, CDKException, CloneNotSupportedException {
        int[][] A = new int[matrixSize][matrixSize];
        count.incrementAndGet();
        if (writeSDF || printSDF) {
            IAtomContainer ac = buildContainer4SDF(addHydrogens(A, hIndex, hydrogens));
            if (coordinates) new StructureDiagramGenerator().generateCoordinates(ac);
            sdfOut.write(ac);
        }
        if (writeSMILES || printSMILES) {
            write2smiles(addHydrogens(A, hIndex, hydrogens));
        }
    }

    /**
     * Calling the generate function for each degree values after the hydrogen distribution.
     *
     * @param degree the degree
     * @return the hydrogens
     */
    public int[] setHydrogens(int[] degree) {
        int[] hydrogens = new int[size];
        for (int i = 0; i < size; i++) {
            hydrogens[i] = firstDegrees[i] - degree[i];
        }
        return hydrogens;
    }

    public void structureGenerator() {
        if (noHydrogen) {
            size = sum(firstOccurrences, firstOccurrences.length - 1);
        } else if (justH) {
            size = hIndex;
        } else {
            size = sum(firstOccurrences, firstOccurrences.length - 2);
        }
        ArrayList<int[]> newDegrees = distributeHydrogens();

        if (multiThread) {
            try {
                new ForkJoinPool(size)
                        .submit(
                                () ->
                                        newDegrees
                                                .parallelStream()
                                                .forEach(new Generation(this)::run))
                        .get();
            } catch (InterruptedException | ExecutionException ex) {
                if (verbose) {
                    Logger.getLogger(MAYGEN.class.getName()).log(Level.SEVERE, null, ex);
                }
                Thread.currentThread().interrupt();
            }
        } else {
            newDegrees.forEach(new Generation(this)::run);
        }
    }

    /** For several calls of the run function, setting the global variables. */
    public void clearGlobals() {
        singleAtom = true;
        onlyDegree2 = true;
        OnSm = true;
        oxygen = 0;
        sulfur = 0;
        graphSize = 0;
        callHydrogenDistributor = false;
        total = 0;
        totalHydrogen = 0;
        size = 0;
        sizePart = 0;
        hIndex = 0;
        count.set(0);
        indexSmiles.set(0);
        indexSdf.set(0);
        matrixSize = 0;
        noHydrogen = false;
        justH = false;
        noHydrogen = false;
        oxygenSulfur = new ArrayList<>();
        symbols = new ArrayList<>();
        occurrences = null;
        symbolArray = null;
        firstSymbols = new ArrayList<>();
        symbols = new ArrayList<>();
        firstOccurrences = null;
    }

    /**
     * Finding the neighbors of a given index.
     *
     * @param index int row (atom) index
     * @param total int number of atoms.
     * @param mat the adjacency matrix
     * @return the nValues
     */
    public Set<Integer> nValues(int index, int total, int[][] mat) {
        Set<Integer> nValues = new HashSet<>();
        nValues.add(index);
        int[] theRow = mat[index];
        for (int i = (index + 1); i < total; i++) {
            if (theRow[i] > 0) {
                nValues.add(i);
            }
        }
        return nValues;
    }

    /**
     * Finding the W values of neighbors in the former connectivity partition.
     *
     * @param nValues the N values
     * @param Kformer the K values of the former step
     * @return the wValues
     */
    public Set<Integer> wValues(Set<Integer> nValues, int[] Kformer) {
        Set<Integer> wValues = new HashSet<>();
        for (Integer i : nValues) {
            wValues.add(Kformer[i]);
        }
        return wValues;
    }

    /**
     * Finding the connectivity partition, so the smallest index in the neighborhood.
     *
     * @param total the total
     * @param wValues the wValues
     * @param kFormer the K values of the former step
     * @return int[]
     */
    public int[] kValues(int total, Set<Integer> wValues, int[] kFormer) {
        int[] kValues = new int[total];
        int min = Collections.min(wValues);
        for (int i = 0; i < total; i++) {
            if (wValues.contains(kFormer[i])) {
                kValues[i] = min;
            } else {
                kValues[i] = kFormer[i];
            }
        }
        return kValues;
    }

    /**
     * Initializing the first connectivity partition.
     *
     * @param total int number of atoms.
     * @return int[]
     */
    public int[] initialKList(int total) {
        int[] k = new int[total];
        for (int i = 0; i < total; i++) {
            k[i] = i;
        }
        return k;
    }

    /* 3.6.2. Connectivity Test */

    /**
     * Test whether an adjacency matrix is connected or disconnected.
     *
     * @param mat int[][] adjacency matrix
     * @param connectivityIndices the connectivityIndices
     * @param learningFromConnectivity the learningFromConnectivity
     * @return boolean
     */
    public boolean connectivityTest(
            int[][] mat, int[] connectivityIndices, boolean[] learningFromConnectivity) {
        learningFromConnectivity[0] = false;
        boolean check = false;
        int[] kValues = initialKList(hIndex);
        Set<Integer> nValues;
        Set<Integer> wValues;
        Set<Integer> zValues = new HashSet<>();
        int zValue = 0;
        for (int i = 0; i < hIndex; i++) {
            nValues = nValues(i, hIndex, mat);
            wValues = wValues(nValues, kValues);
            zValue = Collections.min(wValues);
            zValues.add(zValue);
            kValues = kValues(hIndex, wValues, kValues);
        }
        if (zValue == 0 && allIs0(kValues)) {
            check = true;
        } else {
            setLearningFromConnectivity(
                    zValues, kValues, connectivityIndices, learningFromConnectivity);
        }
        return check;
    }

    /**
     * If matrix is not connected, setting learninfFromConnectivity global variables.
     *
     * @param zValues the minimum index values of each atom's neighborhoods.
     * @param kValues the connectivity partition
     * @param connectivityIndices the connectivityIndices
     * @param learningFromConnectivity the learningFromConnectivity
     */
    public void setLearningFromConnectivity(
            Set<Integer> zValues,
            int[] kValues,
            int[] connectivityIndices,
            boolean[] learningFromConnectivity) {
        learningFromConnectivity[0] = true;
        connectivityIndices[0] = minComponentIndex(zValues, kValues);
        connectivityIndices[1] = hIndex - 1;
    }

    /**
     * Getting the minimum component index. Here, components are compared based on their last
     * indices and sizes.
     *
     * @param zValues the minimum index values of each atom's neighborhoods.
     * @param kValues the connectivity partition
     * @return int
     */
    public int minComponentIndex(Set<Integer> zValues, int[] kValues) {
        int index = findMaximalIndexInComponent(kValues, 0);
        int value;
        for (Integer i : zValues) {
            value = findMaximalIndexInComponent(kValues, i);
            if (value < index) {
                index = value;
            }
        }
        return index;
    }

    /**
     * Finding the maximal index in a component to compare with other components.
     *
     * @param kValues int[] connectivity partition
     * @param value int minimum neighborhood index
     * @return int
     */
    public int findMaximalIndexInComponent(int[] kValues, int value) {
        int maxIndex = hIndex;
        for (int i = hIndex - 1; i > 0; i--) {
            if (kValues[i] == value) {
                maxIndex = i;
                break;
            }
        }
        return maxIndex;
    }

    /**
     * Checks whether all the entries are equal to 0 or not.
     *
     * @param list int[]
     * @return boolean
     */
    public boolean allIs0(int[] list) {
        boolean check = true;
        for (int i = 0; i < list.length; i++) {
            if (list[i] != 0) {
                check = false;
                break;
            }
        }
        return check;
    }

    /**
     * Based on the molecules automorphisms, testing an adjacency matrix is canonical or not.
     *
     * @param initialPartition the initial partition
     * @param r the r
     * @return int
     */
    public int indexYZ(int[] initialPartition, int[] r) {
        int index = 0;
        for (int i = 0; i <= r[0]; i++) {
            index = index + initialPartition[i];
        }
        return index - 1;
    }

    public boolean canonicalTest(
            int[][] A,
            int[] initialPartition,
            int[][] partitionList,
            int[] nonCanonicalIndices,
            ArrayList<ArrayList<Permutation>> formerPermutations,
            int[] partSize,
            int[] r,
            int[] y,
            int[] z,
            int[][] ys,
            int[][] zs,
            boolean[] learningFromCanonicalTest) {
        boolean check = true;
        learningFromCanonicalTest[0] = false;
        int value = indexYZ(initialPartition, r);
        y[0] = ys[0][value];
        z[0] = zs[0][value];
        if (partSize[0] == r[0] && z[0] != 1) {
            z[0] = z[0] - 1;
        }
        clearFormers(false, y[0], partitionList, formerPermutations);
        boolean test = true;
        for (int i = y[0]; i <= z[0]; i++) {
            test =
                    rowCanonicalTest(
                            i,
                            r,
                            A,
                            partitionList[i],
                            canonicalPartition(i, partitionList[i]),
                            initialPartition,
                            partitionList,
                            nonCanonicalIndices,
                            formerPermutations,
                            y,
                            ys,
                            learningFromCanonicalTest);
            if (!test) {
                check = false;
                break;
            }
        }
        clearFormers(check, y[0], partitionList, formerPermutations);
        return check;
    }

    /**
     * When an adjacency matrix is non-canonical, cleaning the formerPermutations and partitionList
     * from the first row of the tested block.
     *
     * @param check boolean canonical test result
     * @param y int first row of the tested block
     * @param partitionList the partitionList
     * @param formerPermutations the formerPermutations
     */
    public void clearFormers(
            boolean check,
            int y,
            int[][] partitionList,
            ArrayList<ArrayList<Permutation>> formerPermutations) {
        if (!check) {
            int formerSize = formerPermutations.size() - 1;
            if (formerSize >= y) {
                formerPermutations.subList(y, formerSize + 1).clear();
            }

            int partitionSize = partitionList.length - 1;
            for (int i = partitionSize; i > y; i--) {
                partitionList[i] = null;
            }
        }
    }

    /**
     * Calculating all candidate permutations for row canonical test.
     *
     * <p>The DFS multiplication of former automorphisms list with the list of cycle transpositions
     * of the row.
     *
     * @param index int row index
     * @param cycles the cycle transpositions
     * @param formerPermutations the formerPermutations
     */
    public void candidatePermutations(
            int index,
            ArrayList<Permutation> cycles,
            ArrayList<ArrayList<Permutation>> formerPermutations) {
        ArrayList<Permutation> newList = new ArrayList<>(cycles);
        if (index != 0) {
            ArrayList<Permutation> formers = formerPermutations.get(index - 1);
            for (Permutation form : formers) {
                if (!form.isIdentity()) {
                    newList.add(form);
                }
            }
            ArrayList<Permutation> newForm = new ArrayList<>();
            for (Permutation frm : formers) {
                if (!frm.isIdentity()) {
                    newForm.add(frm);
                }
            }
            ArrayList<Permutation> newCycles = new ArrayList<>();
            if (cycles.size() != 1) {
                for (Permutation cyc : cycles) {
                    if (!cyc.isIdentity()) {
                        newCycles.add(cyc);
                    }
                }
            }
            for (Permutation perm : newForm) {
                for (Permutation cycle : newCycles) {
                    Permutation newPermutation = cycle.multiply(perm);
                    if (!newPermutation.isIdentity()) {
                        newList.add(newPermutation);
                    }
                }
            }
        }
        formerPermutations.add(index, newList);
    }

    /**
     * Canonical test for a row in the tested block.
     *
     * @param index int row index
     * @param r int block index
     * @param A int[][] adjacency matrix
     * @param partition int[] former partition
     * @param newPartition int[] canonical partition
     * @param initialPartition the initial partition
     * @param partitionList the partitionList
     * @param nonCanonicalIndices the nonCanonicalIndices
     * @param formerPermutations the formerPermutations
     * @param y the y
     * @param ys the ys
     * @param learningFromCanonicalTest the learningFromCanonicalTest
     * @return boolean
     */
    public boolean rowCanonicalTest(
            int index,
            int[] r,
            int[][] A,
            int[] partition,
            int[] newPartition,
            int[] initialPartition,
            int[][] partitionList,
            int[] nonCanonicalIndices,
            ArrayList<ArrayList<Permutation>> formerPermutations,
            int[] y,
            int[][] ys,
            boolean[] learningFromCanonicalTest) {
        boolean check;
        if (!rowDescendingTest(
                index, A, newPartition, nonCanonicalIndices, learningFromCanonicalTest)) {
            check = false;
        } else {
            int value = indexYZ(initialPartition, r);
            y[0] = ys[0][value];
            ArrayList<Permutation> cycles = new ArrayList<>();
            if (partition[size - 1] != 0) {
                Permutation id = new Permutation(size);
                cycles.add(id);
            } else {
                cycles = cycleTranspositions(index, partition);
            }
            candidatePermutations(index, cycles, formerPermutations);
            check = check(index, y[0], size, A, newPartition, formerPermutations);
            if (!check) {
                if (cycles.size() != 1) {
                    getLernenIndices(
                            index,
                            A,
                            cycles,
                            newPartition,
                            nonCanonicalIndices,
                            learningFromCanonicalTest);
                }
            } else {
                addPartition(index, newPartition, A, partitionList);
            }
        }
        return check;
    }

    /**
     * Updating canonical partition list.
     *
     * @param index row index
     * @param newPartition atom partition
     * @param A int[][] adjacency matrix
     * @param partitionList the partitionList
     */
    public void addPartition(int index, int[] newPartition, int[][] A, int[][] partitionList) {
        if (newPartition[size - 1] != 0) {
            partitionList[index + 1] = newPartition;
        } else {
            partitionList[index + 1] = refinedPartitioning(newPartition, A[index]);
        }
    }

    /**
     * Refining the input partition based on the row entries.
     *
     * @param partition int[] atom partition
     * @param row int[] row
     * @return int[]
     */
    public int[] refinedPartitioning(int[] partition, int[] row) {
        int[] refined = new int[size];
        int index = 0;
        int count = 1;
        int refinedIndex = 0;
        int limit = findZeros(partition);
        for (int s = 0; s < limit; s++) {
            if (partition[s] != 1) {
                for (int i = index; i < partition[s] + index - 1; i++) {
                    if (i + 1 < partition[s] + index - 1) {
                        if (row[i] == row[i + 1]) {
                            count++;
                        } else {
                            refined[refinedIndex] = count;
                            refinedIndex++;
                            count = 1;
                        }
                    } else {
                        if (row[i] == row[i + 1]) {
                            count++;
                            refined[refinedIndex] = count;
                            refinedIndex++;
                            count = 1;
                        } else {
                            refined[refinedIndex] = count;
                            refinedIndex++;
                            refined[refinedIndex] = 1;
                            refinedIndex++;
                            count = 1;
                        }
                    }
                }
                index = index + partition[s];
            } else {
                index++;
                refined[refinedIndex] = 1;
                refinedIndex++;
                count = 1;
            }
        }
        return refined;
    }

    /**
     * For a row given by index, detecting the other row to compare in the block. For the detection
     * of the next row index, cycle transposition is used.
     *
     * @param index int row index
     * @param A int[][] adjacency matrix
     * @param cycleTransposition Permutation cycle transposition
     * @return int[]
     */
    public int[] row2compare(int index, int[][] A, Permutation cycleTransposition) {
        int[] array = cloneArray(A[findIndex(index, cycleTransposition)]);
        return actArray(array, cycleTransposition);
    }

    /**
     * With the cycle permutation, mapping the row index to another row in the block.
     *
     * @param index int row index
     * @param cycle Permutation cycle transposition
     * @return int
     */
    public int findIndex(int index, Permutation cycle) {
        int size = cycle.size();
        int output = 0;
        for (int i = 0; i < size; i++) {
            if (cycle.get(i) == index) {
                output = i;
                break;
            }
        }
        return output;
    }

    /**
     * Cloning int array
     *
     * @param array int[] array
     * @return int[]
     */
    public int[] cloneArray(int[] array) {
        return array.clone();
    }

    /**
     * Calculating the canonical permutation of a row.
     *
     * <p>In a block, the original and the other rows are compared; if there is a permutation
     * mapping rows to each other, canonical permutation, else id permutation is returned.
     *
     * @param originalRow int[] original row
     * @param rowToCheck int[] row to compare with
     * @param partition int[] partition
     * @return Permutation
     */
    public Permutation getCanonicalPermutation(
            int[] originalRow, int[] rowToCheck, int[] partition) {
        int[] cycles = getCanonicalPermutation2(partition, originalRow, rowToCheck);
        int[] perm = new int[size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == cycles[j]) {
                    perm[i] = j;
                }
            }
        }
        return new Permutation(perm);
    }

    /**
     * Calculating the canonical permutation of a row.
     *
     * <p>In a block, the original and the other rows are compared; if there is a permutation
     * mapping rows to each other, canonical permutation, else id permutation is returned.
     *
     * @param partition int[] partition
     * @param max int[] max
     * @param check int[] check
     * @return int[]
     */
    public int[] getCanonicalPermutation2(int[] partition, int[] max, int[] check) {
        int[] values = idValues(sum(partition));
        int i = 0;

        if (!equalSetCheck(max, check, partition)) {
            return values;
        } else {
            int limit = findZeros(partition);
            for (int s = 0; s < limit; s++) {
                int[] can = getBlocks(max, i, partition[s] + i);
                int[] non = getBlocks(check, i, partition[s] + i);
                values = getCyclesList(can, non, i, values);
                i = i + partition[s];
            }
            return values;
        }
    }

    public int[] getCyclesList(int[] max, int[] non, int index, int[] values) {
        int i = 0;
        int permutationIndex = 0;
        while (i < max.length && max[i] != 0) {
            if (max[i] != non[i]) {
                permutationIndex = findMatch(max, non, max[i], i);
                if (i != permutationIndex) {
                    non = permuteArray(non, i, permutationIndex);
                }
                int temp = values[i + index];
                values[i + index] = values[permutationIndex + index];
                values[permutationIndex + index] = temp;
            }
            i++;
        }
        return values;
    }

    /**
     * @param max the max
     * @param non the non
     * @param value the vvalue
     * @param start the start
     * @return the index
     */
    public int findMatch(int[] max, int[] non, int value, int start) {
        int size = non.length;
        int index = start;
        for (int i = start; i < size; i++) {
            if (non[i] == value) {
                if (max[i] != non[i]) {
                    index = i;
                    break;
                }
            }
        }
        return index;
    }

    public Permutation getEqualPerm(
            Permutation cycleTransposition, int index, int[][] A, int[] newPartition) {
        int[] check = row2compare(index, A, cycleTransposition);
        return getCanonicalPermutation(A[index], check, newPartition);
    }

    public Permutation getCanonicalCycle(
            int index, int total, int[][] A, int[] newPartition, Permutation cycleTransposition) {
        Permutation canonicalPermutation = idPermutation(total);
        if (!equalRowsCheck(index, A, cycleTransposition, canonicalPermutation)) {
            canonicalPermutation = getEqualPerm(cycleTransposition, index, A, newPartition);
        }
        return canonicalPermutation;
    }

    public boolean check(
            int index,
            int y,
            int total,
            int[][] A,
            int[] newPartition,
            ArrayList<ArrayList<Permutation>> formerPermutations) {
        boolean check = true;
        ArrayList<Permutation> formerList = new ArrayList<>();
        ArrayList<Permutation> form = formerPermutations.get(index);
        for (Permutation permutation : form) {
            boolean biggest = setBiggest(index, A, permutation, newPartition);
            if (biggest) {
                Permutation canonicalPermutation =
                        getCanonicalCycle(index, total, A, newPartition, permutation);
                int[] test = row2compare(index, A, permutation);
                test = actArray(test, canonicalPermutation);
                if (descendingOrderUpperMatrixCheck(index, newPartition, A[index], test)) {
                    if (canonicalPermutation.isIdentity()) {
                        if (equalSetCheck2(newPartition, A[index], test)) {
                            formerList.add(permutation);
                        }
                    } else {
                        Permutation newPermutation = canonicalPermutation.multiply(permutation);
                        formerList.add(newPermutation);
                    }
                } else {
                    formerList.clear();
                    check = false;
                    break;
                }
            } else {
                formerList.clear();
                check = false;
                break;
            }
        }
        if (check) {
            formerPermutations.get(index).clear();
            formerPermutations.set(index, formerList);
        }
        return check;
    }

    public ArrayList<Permutation> cycleTranspositions(int index, int[] partition) {
        ArrayList<Permutation> perms = new ArrayList<>();
        int lValue = LValue(partition, index);
        int[] values;
        int former;
        for (int i = 0; i < lValue; i++) {
            values = idValues(size);
            former = values[index];
            values[index] = values[index + i];
            values[index + i] = former;
            Permutation p = new Permutation(values);
            perms.add(p);
        }
        return perms;
    }

    /**
     * Grund Thesis 3.3.3. To calculate the number of conjugacy classes, used in cycle transposition
     * calculation.
     *
     * @param partEx int[] former atom partition
     * @param degree the degree
     * @return the LValue
     */
    public int LValue(int[] partEx, int degree) {
        return (sum(partEx, (degree)) - (degree));
    }

    /**
     * To get the canonical partition like in Grund Thesis 3.3.11
     *
     * @param i int row index
     * @param partition int[] partition
     * @return int[]
     */
    public int[] canonicalPartition(int i, int[] partition) {
        return partitionCriteria(partition, i + 1);
    }

    /**
     * Add number of 1s into an ArrayList
     *
     * @param list the integer array
     * @param number the number
     */
    public void addOnes(int[] list, int number) {
        for (int i = 0; i < number; i++) {
            list[i] = 1;
        }
    }

    public int findZeros(int[] array) {
        int index = size;
        for (int i = 0; i < size; i++) {
            if (array[i] == 0) {
                index = i;
                break;
            }
        }
        return index;
    }

    /**
     * Grund Thesis 3.3.2 Partitioning criteria (DONE)
     *
     * @param partEx the former partition
     * @param degree degree of the partitioning.
     * @return the partition criteria
     */
    public int[] partitionCriteria(int[] partEx, int degree) {
        int[] partNew = new int[size];
        int limit = findZeros(partEx);
        if (zero(partEx)) {
            addOnes(partNew, degree);
            int index = degree;
            int oldValue = partEx[degree - 1];
            if (oldValue > 1) {
                partNew[index] = oldValue - 1;
                index++;
                for (int k = degree; k < limit; k++) {
                    partNew[index] = partEx[k];
                    index++;
                }
            } else if (oldValue == 1) {
                for (int k = degree; k < limit; k++) {
                    partNew[index] = partEx[k];
                    index++;
                }
            }
            return partNew;
        } else {
            return partEx;
        }
    }

    public void orderDegreeSymbols(
            int[] degree, String[] symbol, int index0, int index1, int[] hydrogens) {
        int temp = 0;
        int temp2 = 0;
        for (int i = index0; i < index1; i++) {
            for (int j = i + 1; j < index1; j++) {
                if (degree[i] > degree[j]) {
                    swap(symbol, i, j);
                    temp = degree[i];
                    degree[i] = degree[j];
                    degree[j] = temp;
                    temp2 = hydrogens[i];
                    hydrogens[i] = hydrogens[j];
                    hydrogens[j] = temp2;
                }
            }
        }
    }

    public void swap(String[] array, int i, int j) {
        String temp;
        temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    public void swap(int[] array, int i, int j) {
        int swapString = array[i];
        array[i] = array[j];
        array[j] = swapString;
    }

    public int[] sortWithPartition(
            int[] partitionList, int[] degrees, String[] symbols, int[] hydrogens) {
        int[] partition = buildArray(partitionList);
        int size = partition.length;
        for (int n = 0; n < size; n++) {
            for (int m = 0; m < (size - 1) - n; m++) {
                if ((partition[m] > partition[m + 1])) {
                    swap(partition, m, (m + 1));
                    swap(degrees, m, (m + 1));
                    swap(hydrogens, m, (m + 1));
                    swap(symbols, m, (m + 1));
                }
            }
        }
        reOrder(partition, degrees, symbols, hydrogens);
        return initialPartition(partition);
    }

    public int[] initialPartition(int[] partition) {
        int index = 0;
        int index2 = 0;
        int part;
        int[] init = new int[size];
        while (index != hIndex) {
            part = partition[index];
            init[index2++] = part;
            index += part;
        }
        return init;
    }

    public int[] buildArray(int[] partition) {
        int[] partitionArray = new int[sum(partition)];
        int index = 0;
        for (int p : partition) {
            for (int i = 0; i < p; i++) {
                partitionArray[index] = p;
                index++;
            }
        }
        return partitionArray;
    }

    public void reOrder(int[] partition, int[] degrees, String[] symbols, int[] hydrogens) {
        int index = 0;
        int part;
        while (index != (hIndex)) {
            part = partition[index];
            orderDegreeSymbols(degrees, symbols, index, (index + part), hydrogens);
            index = index + part;
        }
    }

    public void parseArgs(String[] args) throws ParseException {
        Options options = setupOptions();
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);
            this.formula = cmd.getOptionValue("formula");
            if (cmd.hasOption("outputFile")) {
                String filedir = cmd.getOptionValue("outputFile");
                this.filedir = Objects.isNull(filedir) ? "." : filedir;
                if (cmd.hasOption("smi")) {
                    this.writeSMILES = true;
                }
                if (cmd.hasOption("sdf")) {
                    this.writeSDF = true;
                } else if (cmd.hasOption("sdfCoord")) {
                    this.writeSDF = true;
                    this.coordinates = true;
                }
            } else {
                if (cmd.hasOption("smi") && !cmd.hasOption("sdf")) {
                    this.printSMILES = true;
                }
                if (cmd.hasOption("sdf")) {
                    this.printSDF = true;
                } else if (cmd.hasOption("sdfCoord")) {
                    this.printSDF = true;
                    this.coordinates = true;
                }
            }
            if (cmd.hasOption("verbose")) this.verbose = true;
            if (cmd.hasOption("tsvoutput")) this.tsvoutput = true;
            if (cmd.hasOption("multithread")) this.multiThread = true;
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(null);
            String header =
                    "\nGenerates molecular structures for a given molecular formula."
                            + "\nThe input is a molecular formula string."
                            + "\n\nFor example 'C2OH4'."
                            + "\n\nIf user wants to store output file in a specific directory, that is needed to be specified."
                            + "\nIt is also possible to generate SMILES instead of an SDF file, but it slows down"
                            + "the generation time. For this, use the '-smi' option."
                            + "\n\n";
            String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/MAYGEN";
            formatter.printHelp("java -jar MAYGEN.jar", header, options, footer, true);
            throw new ParseException("Problem parsing command line");
        }
    }

    public Options setupOptions() {
        Options options = new Options();
        Option formula =
                Option.builder("f")
                        .required(true)
                        .hasArg()
                        .longOpt("formula")
                        .desc("formula (required)")
                        .build();
        options.addOption(formula);
        Option verbose =
                Option.builder("v")
                        .required(false)
                        .longOpt("verbose")
                        .desc("print message")
                        .build();
        options.addOption(verbose);
        Option tvsoutput =
                Option.builder("t")
                        .required(false)
                        .longOpt("tsvoutput")
                        .desc(
                                "Output formula, number of structures and execution time in CSV format."
                                        + "In multithread, the 4th column in the output is the number of threads.")
                        .build();
        options.addOption(tvsoutput);
        Option filedir =
                Option.builder("o")
                        .required(false)
                        .hasArg()
                        .optionalArg(true)
                        .longOpt("outputFile")
                        .desc("Store output file")
                        .build();
        options.addOption(filedir);
        Option multithread =
                Option.builder("m")
                        .required(false)
                        .longOpt("multithread")
                        .desc("Use multi thread")
                        .build();
        options.addOption(multithread);
        Option smiles =
                Option.builder("smi")
                        .required(false)
                        .longOpt("SMILES")
                        .desc("Output in SMILES format")
                        .build();
        options.addOption(smiles);
        Option sdf =
                Option.builder("sdf")
                        .required(false)
                        .longOpt("SDF")
                        .desc("Output in SDF format")
                        .build();
        options.addOption(sdf);
        Option coordinates =
                Option.builder("sdfCoord")
                        .required(false)
                        .longOpt("coordinates")
                        .desc("Output in SDF format with atom coordinates")
                        .build();
        options.addOption(coordinates);
        return options;
    }

    public int[][] buildOnSm(int[] ar) {
        int[][] mat = new int[graphSize][graphSize];
        int oxygen = 0;
        int sulfur = graphSize / 2;
        int x = 0;
        int y = 0;

        if (ar[graphSize - 1] == 0) {
            x = sulfur - 1;
        } else {
            x = graphSize - 1;
        }

        if (ar[0] == 0) {
            y = 0;
        } else {
            y = sulfur;
        }
        int temp = 0;
        if (x > y) {
            temp = x;
            x = y;
            y = temp;
        }
        mat[x][y] = 1;
        for (int i = 0; i < graphSize - 1; i++) {
            if (ar[i] == 0) {
                x = oxygen;
                oxygen++;
            } else {
                x = sulfur;
                sulfur++;
            }
            if (ar[i + 1] == 0) {
                y = oxygen;
            } else {
                y = sulfur;
            }
            if (x > y) {
                temp = x;
                x = y;
                y = temp;
            }
            mat[x][y] = 1;
        }
        return mat;
    }

    public void write2smiles(int[][] mat) throws IOException, CloneNotSupportedException {

        IAtomContainer atomContainer = buildAtomContainer(mat);
        try {
            String smilesString = smilesGenerator.create(atomContainer);
            smilesOut.write(smilesString + "\n");
        } catch (CDKException ex) {
            if (verbose) {
                Logger.getLogger(MAYGEN.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    /**
     * Building an atom container from a string of atom-implicit hydrogen information.
     *
     * @param ac the IAtomContainer
     * @param symbolArrayCopy the symbolArrayCopy
     * @return the IAtomContainer
     */
    public IAtomContainer initAC(IAtomContainer ac, String[] symbolArrayCopy) {
        for (int i = 0; i < symbolArrayCopy.length; i++) {
            ac.addAtom(new Atom(symbolArrayCopy[i]));
        }

        for (IAtom atom : ac.atoms()) {
            atom.setImplicitHydrogenCount(0);
        }
        return ac;
    }

    public void initSingleAC() {
        atomContainer = builder.newInstance(IAtomContainer.class);
        for (int i = 0; i < symbolArray.length; i++) {
            atomContainer.addAtom(new Atom(symbolArray[i]));
        }

        for (IAtom atom : atomContainer.atoms()) {
            atom.setImplicitHydrogenCount(0);
        }
    }

    public void intAC(String formula) {
        String[] atoms = formula.split(LETTERS_FROM_A_TO_Z);
        ArrayList<String> symbolList = new ArrayList<>();
        String[] info;
        int occur;
        String symbol;
        for (String atom : atoms) {
            info = atom.split(NUMBERS_FROM_0_TO_9, 2);
            symbol = info[0];
            occur = atomOccurrence(info);
            for (int i = 0; i < occur; i++) {
                symbolList.add(symbol);
            }
        }
        atomContainer = builder.newInstance(IAtomContainer.class);
        for (String s : symbolList) {
            atomContainer.addAtom(new Atom(s));
        }
        for (IAtom atom : atomContainer.atoms()) {
            atom.setImplicitHydrogenCount(0);
        }
    }
    /**
     * Building an atom container from a string of atom-implicit hydrogen information.
     *
     * @param symbol the symbol
     */
    public void initAC(String symbol) {
        atomContainer = builder.newInstance(IAtomContainer.class);
        for (int i = 0; i < matrixSize; i++) {
            atomContainer.addAtom(new Atom(symbol));
        }

        for (IAtom atom : atomContainer.atoms()) {
            atom.setImplicitHydrogenCount(0);
        }
    }

    /**
     * Building an atom container for an adjacency matrix.
     *
     * @param ac the IAtomContainer
     * @param mat int[][] adjacency matrix
     * @return IAtomContainer
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     */
    public IAtomContainer buildContainer4SDF(IAtomContainer ac, int[][] mat)
            throws CloneNotSupportedException {
        IAtomContainer ac2 = ac.clone();
        for (int i = 0; i < mat.length; i++) {
            for (int j = i + 1; j < mat.length; j++) {
                if (mat[i][j] == 1) {
                    ac2.addBond(i, j, Order.SINGLE);
                } else if (mat[i][j] == 2) {
                    ac2.addBond(i, j, Order.DOUBLE);
                } else if (mat[i][j] == 3) {
                    ac2.addBond(i, j, Order.TRIPLE);
                }
            }
        }

        ac2 = AtomContainerManipulator.removeHydrogens(ac2);
        return ac2;
    }

    /**
     * Building an atom container for an adjacency matrix.
     *
     * @param mat int[][] adjacency matrix
     * @return IAtomContainer
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     */
    public IAtomContainer buildContainer4SDF(int[][] mat) throws CloneNotSupportedException {
        IAtomContainer ac2 = atomContainer.clone();
        for (int i = 0; i < mat.length; i++) {
            for (int j = i + 1; j < mat.length; j++) {
                if (mat[i][j] == 1) {
                    ac2.addBond(i, j, Order.SINGLE);
                } else if (mat[i][j] == 2) {
                    ac2.addBond(i, j, Order.DOUBLE);
                } else if (mat[i][j] == 3) {
                    ac2.addBond(i, j, Order.TRIPLE);
                }
            }
        }
        return AtomContainerManipulator.removeHydrogens(ac2);
    }

    /**
     * Building an atom container from a string of atom-implicit hydrogen information.
     *
     * @param mat the mat
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     * @return the IAtomContainer
     */
    public IAtomContainer buildAtomContainer(int[][] mat) throws CloneNotSupportedException {
        IAtomContainer atomContainer = this.builder.newAtomContainer();
        for (String s : this.symbolArray) {
            atomContainer.addAtom(new Atom(s));
        }

        for (IAtom atom : atomContainer.atoms()) {
            atom.setImplicitHydrogenCount(0);
        }

        atomContainer = buildContainer4SDF(atomContainer, mat);

        return atomContainer;
    }

    public void degree2graph() throws IOException, CDKException, CloneNotSupportedException {
        int[][] mat = new int[matrixSize][matrixSize];
        mat[0][1] = 1;
        mat[0][2] = 1;
        for (int i = 1; i < matrixSize - 2; i++) {
            mat[i][i + 2] = 1;
        }
        mat[matrixSize - 2][matrixSize - 1] = 1;
        count.incrementAndGet();
        String symbol = "";
        if (oxygen == 0) {
            symbol = "S";
        } else {
            symbol = "O";
        }
        initAC(symbol);
        if (writeSDF || printSDF) {
            IAtomContainer ac = buildContainer4SDF(mat);
            if (coordinates) new StructureDiagramGenerator().generateCoordinates(ac);
            sdfOut.write(ac);
        }
        if (writeSMILES || printSMILES) {
            smilesOut.write("Formula is not supported" + "\n");
        }
    }

    /*
     The following functions are for the distribution of given number of sulfurs and oxygens in
     regular graphs as node labelling. These functions are the Java implementation of Sawada's
     method[1][2] in mathematical chemistry.

     <p>References:

     <p>[1] Sawada J. Generating bracelets in constant amortized time. SIAM Journal on Computing.
     2001;31(1):259-68. [2] http://www.cis.uoguelph.ca/~sawada/prog/necklaces.c
    */
    /**
     * The function to compare a node labelling array when the number of left consecutive entries is
     * equal to the right consecutive entries.
     *
     * <p>From this comparison, function returns :
     *
     * <p>1: if node labelling is same as its reversal.
     *
     * <p>0: if node labelling is smaller than its reversal.
     *
     * <p>-1: if node labelling is bigger than its reversal.
     *
     * <p>Reversal check helps to avoid duplicates, easy way of rotational symmetry check.
     *
     * @param length int node labelling length
     * @param index int starting index in the array
     * @return the result
     */
    public int reverseComparison(int length, int index) {
        for (int i = index + 1; i <= (length + 1) / 2; i++) {
            if (nodeLabels[i] < nodeLabels[length - i + 1]) {
                return 0;
            } else if (nodeLabels[i] > nodeLabels[length - i + 1]) {
                return -1;
            }
        }
        return 1;
    }

    public int[] build() {
        int[] arr = new int[graphSize];
        if (graphSize >= 0) {
            System.arraycopy(nodeLabels, 1, arr, 0, graphSize);
        }
        return arr;
    }

    /**
     * Main function for the distribution of atom symbols: O and S for OnSm form formulae.
     *
     * @param oxy int number of oxygens to distribute
     * @param sul int number of sulfur to distribute
     * @param nextSize int length of the next labelling
     * @param currentSize int length of the current labelling
     * @param reversedLength int longest node labelling, equal to its reversal
     * @param leftEquivalents int the number of consequtively equivalent values at the left side of
     *     the array
     * @param rightEquivalents int the number of consequtively equivalent values at the left side of
     *     the array
     * @param reversalIsSmaller boolean from the reversal comparison, using the boolean variable to
     *     know reveral is smaller than the bode labelling or not.
     * @throws CloneNotSupportedException in case of CloneNotSupportedException
     * @throws CDKException in case of CDKException
     * @throws IOException in case of IOException
     */
    public void distributeSymbols(
            int oxy,
            int sul,
            int nextSize,
            int currentSize,
            int reversedLength,
            int leftEquivalents,
            int rightEquivalents,
            boolean reversalIsSmaller)
            throws CDKException, CloneNotSupportedException, IOException {

        if (2 * (nextSize - 1) > (graphSize + reversedLength)) {
            if (nodeLabels[nextSize - 1] > nodeLabels[graphSize - nextSize + 2 + reversedLength])
                reversalIsSmaller = false;
            else if (nodeLabels[nextSize - 1]
                    < nodeLabels[graphSize - nextSize + 2 + reversedLength])
                reversalIsSmaller = true;
        }
        if (nextSize > graphSize) {
            if (!reversalIsSmaller && (graphSize % currentSize) == 0) {
                count.incrementAndGet();
                if (writeSDF || printSDF) {
                    IAtomContainer ac = buildContainer4SDF(buildOnSm(build()));
                    if (coordinates) new StructureDiagramGenerator().generateCoordinates(ac);
                    sdfOut.write(ac);
                }
                if (writeSMILES || printSMILES) {
                    write2smiles(buildOnSm(build()));
                }
            }
        } else {
            int oxy2 = oxy;
            int sul2 = sul;

            nodeLabels[nextSize] = nodeLabels[nextSize - currentSize];

            if (nodeLabels[nextSize] == 0) {
                oxy2--;
            } else {
                sul2--;
            }
            if (nodeLabels[nextSize] == nodeLabels[1]) {
                rightEquivalents++;
            } else {
                rightEquivalents = 0;
            }
            if ((leftEquivalents == (nextSize - 1))
                    && (nodeLabels[nextSize - 1] == nodeLabels[1])) {
                // left consecutive element number incremention.
                leftEquivalents++;
            }

            if (oxy2 >= 0
                    && sul2 >= 0
                    && !(nextSize == graphSize
                            && leftEquivalents != graphSize
                            && nodeLabels[graphSize] == nodeLabels[1])) {
                if (leftEquivalents == rightEquivalents) {
                    int reverse = reverseComparison(nextSize, leftEquivalents);
                    if (reverse == 0) {
                        distributeSymbols(
                                oxy2,
                                sul2,
                                nextSize + 1,
                                currentSize,
                                reversedLength,
                                leftEquivalents,
                                rightEquivalents,
                                reversalIsSmaller);

                    } else if (reverse == 1) {
                        distributeSymbols(
                                oxy2,
                                sul2,
                                nextSize + 1,
                                currentSize,
                                nextSize,
                                leftEquivalents,
                                rightEquivalents,
                                false);
                    }
                } else {
                    distributeSymbols(
                            oxy2,
                            sul2,
                            nextSize + 1,
                            currentSize,
                            reversedLength,
                            leftEquivalents,
                            rightEquivalents,
                            reversalIsSmaller);
                }
            }

            if (leftEquivalents == nextSize) {
                leftEquivalents--;
            }

            if (nodeLabels[nextSize - currentSize] == 0 && sul > 0) {
                nodeLabels[nextSize] = 1;

                if (nextSize == 1) {
                    distributeSymbols(
                            oxy, sul - 1, nextSize + 1, nextSize, 1, 1, 1, reversalIsSmaller);
                } else {
                    distributeSymbols(
                            oxy,
                            sul - 1,
                            nextSize + 1,
                            nextSize,
                            reversedLength,
                            leftEquivalents,
                            0,
                            reversalIsSmaller);
                }
            }
        }
    }

    public void distributeSulfurOxygen()
            throws CDKException, CloneNotSupportedException, IOException {
        graphSize = oxygen + sulfur;
        nodeLabels = new int[graphSize + 1];
        nodeLabels[0] = 0;
        intAC(formula);
        distributeSymbols(oxygen, sulfur, 1, 1, 0, 0, 0, false);
    }

    public static void main(String[] args) {
        MAYGEN gen = new MAYGEN();
        try {
            gen.parseArgs(args);
            gen.run();
        } catch (Exception ex) {
            if (gen.verbose) {
                Logger.getLogger(MAYGEN.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
}
