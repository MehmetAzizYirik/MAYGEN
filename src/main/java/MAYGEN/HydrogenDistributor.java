package MAYGEN;
/**
 * MIT License
 *
 * Copyright (c) 2018 Mehmet Aziz Yirik
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


/**
 * This class is for the distribution of hydrogens into chemical structures.
 * For a given chemical formula, all the possible distributions of the hydrogen
 * atoms to the hetero atoms are generated. 
 * 
 * @author Mehmet Aziz Yirik
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

public class HydrogenDistributor {
	public static Map<Integer, Integer> capacities;
	public static int isotopes;
	public static int[] capacity;
	public static int[] valences;
	public static int totalHydrogen; // Total number of hydrogens.
	public static int[] totalAtom; // Total number of atoms.
	public static int hydrogens2distribute;
	
	static {
		//The atom capacities from MOLGEN book. Capacity of an atom equals to 
		capacities = new HashMap<Integer, Integer>();
		capacities.put(6, 3);
		capacities.put(7, 2);
		capacities.put(8, 1);
		capacities.put(16, 1);
		capacities.put(15, 2);
		capacities.put(9, 0);
		capacities.put(53, 0);
		capacities.put(1, 0);
		capacities.put(17, 0);
		capacities.put(35, 0);
		capacities.put(53, 0);
		
	}
	
	/**
	 * The basic functions used in the hydrogen distributor.
	 */
	
	public static int[] addElement(int[] a, int e) {
        a  = Arrays.copyOf(a, a.length + 1);
        a[a.length - 1] = e;
        return a;
    }
	
	public static int sum(int[] partition, int index) {
		int sum=0;
		for(int i=0;i<=index;i++) {
			sum=sum+partition[i];
		}
	    return sum;
	}
	
	public static int sum(ArrayList<Integer> list, int index) {
		int sum=0;
		for(int i=0;i<=index;i++) {
			sum=sum+list.get(i);
		}
		return sum;
	}
	
	public static int[] setValues(ArrayList<Integer> partition, int[] degrees) {
		int partitionSize= partition.size();
		int[] capacity = new int[partitionSize];
		int[] valences = new int[partitionSize];
		int[] totalAtom = new int[partitionSize];
		int i=0;
		int sum=0;
		for(int j=0;j<partitionSize;j++) {
			totalAtom[i]=partition.get(i);
			sum=sum(partition,i);
			valences[i]=degrees[sum-1]-1;
			capacity[i]=(degrees[sum-1]-1)*partition.get(i);
			i++;
		}

		HydrogenDistributor.capacity=capacity;
		HydrogenDistributor.valences=valences;
		HydrogenDistributor.totalAtom=totalAtom;
		return capacity;
	}
	
	public static int sum(int[] array) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+array[i];
		}
		return sum;
	}
	
	public static int[] mergeArrays(ArrayList<int[]> arrays) {
		int size = 0;
		for (int[] array : arrays) {
			size += array.length;
		}
		int[] mergedArray = new int[size];
		int index = 0;
		for (int[] array : arrays) {
			for (int i : array) {
				mergedArray[index++] = i;
		    }
		}
		return mergedArray;
	}
	
	public static int[] arraySum(int[] a, int[] b) {
		ArrayList<int[]> arrays= new ArrayList<int[]>();
		arrays.add(a);
		arrays.add(b);
		return mergeArrays(arrays);
	}

	public static ArrayList<int[]> combineArrays(LinkedList<ArrayList <int[]>> lists) {
		ArrayList<int[]> comb = new ArrayList<int[]>();
	    for (int[] s: lists.removeFirst()) {
	    	comb.add(s);
	    }
	    while (!lists.isEmpty()) {
	        List<int[]> list = lists.removeFirst();
			ArrayList<int[]> newComb =  new ArrayList<int[]>();
	        for (int[] arr1: comb) { 
	            for (int[] arr2 : list) { 
	            	newComb.add(arraySum(arr1,arr2));
	            }
	        }
	        comb = newComb;
	    }
	    return comb;
	}
	
	/**
	 * To initialise the inputs and run the functions while recording the duration time.
	 */
	
	public static ArrayList<int[]> run(ArrayList<Integer> partition, int[] degrees) throws CloneNotSupportedException {
		int partitionSize= partition.size();
		int hydrogen= partition.get(partitionSize-1);
		HydrogenDistributor.isotopes=partitionSize-1;
		setValues(partition,degrees);
		HydrogenDistributor.totalHydrogen=hydrogen;
		ArrayList<int[]> result= new ArrayList<int[]>();
		if(isotopes==1) {
			ArrayList<int[]> iarrays= new ArrayList<int[]>();
			int[] array = new int[0];
			HydrogenDistributor.hydrogens2distribute=totalHydrogen;
			distribute(iarrays,totalHydrogen,array,valences[0],totalAtom[0]);
			result= iarrays;
		}else {
			ArrayList<int[]> distributions= new ArrayList<int[]>();
			for(int[] dene:partition(totalHydrogen,isotopes,0)){
				LinkedList<ArrayList<int[]>> lists = new LinkedList<ArrayList <int[]>>();
				for(int i=0;i<dene.length;i++) {
					HydrogenDistributor.hydrogens2distribute=dene[i];
					ArrayList<int[]> iarrays= new ArrayList<int[]>();
					int[] array = new int[0];
					distribute(iarrays,dene[i],array,valences[i],totalAtom[i]);
					lists.add(iarrays);
				}
				ArrayList<int[]> combined=combineArrays(lists);
				distributions.addAll(combined);
			}
			result=distributions;
		}
		return result;
	}
	
	/**
	 * These functions are built for the integer partitioning problem.
	 */
	
	public static List<int[]> partition(int n, int d,int depth) {
		if(d==depth) {
			List<int[]> array= new ArrayList<int[]>();
			int[] take=new int[0];
			array.add(take);
			return array;
		}
		return buildArray(n,d,depth);
		
	}
	
	public static ArrayList<int[]> buildArray(int n,int d, int depth){
		ArrayList<int[]> array= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,n);
		for(int i:range.toArray()) {
			for(int[] item: partition(n-i,d,depth+1)) {
				if(i<=capacity[item.length]) {
					item=addElement(item,i);
			        if(item.length==d) {
			        	if(sum(item)==totalHydrogen) {
			        		array.add(item);
			        	}
			        }else {
			        	array.add(item);
			        }
				}
			}
		}
		return array;
	}
	
	public static int[] addZeros(int[] array, int zeros) {
		for(int i=0;i<zeros;i++) {
			array=addElement(array,0);
		}
		return array;
	}
	
	public static int[] descendingOrderArray(int[] arr) {
		return Arrays.stream(arr).boxed().sorted().mapToInt(Integer::intValue).toArray();
	}
	
	public static void distribute(ArrayList<int[]> arrays,int hydrogen,int[]arr,int valence, int numAtom) {
		if(hydrogen==0 && sum(arr)==hydrogens2distribute){
			if(arr.length!=numAtom) {
				arr=addZeros(arr,(numAtom-arr.length));
			}
			arr=descendingOrderArray(arr);
			arrays.add(arr);
		}else if((numAtom-arr.length)==1) {
			int add=Math.min(hydrogen,valence);
			hydrogen=hydrogen-add;
			if(arr.length==0) {
				distribute(arrays,0,addElement(arr,add),valence,numAtom); 
			}
			if((arr.length)>0) {
				if(arr[arr.length-1]<=add){ 
					distribute(arrays,0,addElement(arr,add),valence,numAtom);
				}
			}
		}else { 
			for(int i = Math.min(valence,hydrogen); i > 0; i--) {
				if(arr.length==0) {
					distribute(arrays,(hydrogen-i),addElement(arr,i),valence,numAtom); 
				}
				if((arr.length)>0) {
					if(arr[arr.length-1]<=i){ 
						distribute(arrays,(hydrogen-i),addElement(arr,i),valence,numAtom);
					}
				}
			}
		}
	}
}

