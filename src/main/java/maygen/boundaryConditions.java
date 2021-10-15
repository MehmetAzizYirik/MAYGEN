package maygen;

/**
 * The class includes the early boundary conditions for the chemical graph generation. 
 */

public class boundaryConditions {
	
    /**
     * No triple bonds
     * 
     * @param mat 	int[][] adjacency matrix
     * @return 		boolean 
     */
    
    public boolean detectTripleBonds(int[][] mat) {
    	boolean check=false;
    	int length=mat.length;
    	outer: for(int i=0;i<length;i++) {
    		for(int j=0;j<length;j++) {
    			if(mat[i][j]==3) {
        			check=true;
        			break outer;
        		}
    		}
    	}
    	return check;
    }
    
    /**
     * No adjacent double bonds
     * 
     * @param mat 	int[][] adjacency matrix
     * @return 		boolean 
     */
    
    public boolean detectAdjacentDoubleBonds(int[][] mat) {
    	boolean check =false;
    	int count = 0;
    	outer : for(int i=0;i<mat.length;i++) {
    		count=0;
    		for(int j=0; j<mat.length; j++) {
    			if(mat[i][j]==2) count++;
    		}
    		if(count>=2){
    			check=true;
    			break outer;
    		}
    	}
    	return check;
    }
    
    /**
     * No allenes
     * 
     * @param mat 	int[][] adjacency matrix
     * @return 		boolean 
     */
    
    public boolean detectAllenes(int[][] mat, String[] symbols) {
    	boolean check =false;
    	int count = 0;
    	outer : for(int i=0;i<mat.length;i++) {
    		count=0;
    		if(symbols[i].equals("C")) {
    			for(int j=0; j<mat.length; j++) {
        			if(mat[i][j]==2 && symbols[j].equals("C")) {
        				count++;
        			}
        		}
    			if(count>=2){
    				check=true;
    				break outer;
    			}
    		}
    	}
    	return check;
    }
    
    /**
     * After defining the above early boundary conditions, they need to be added to the boundaryConditionCheck function. 
     * 
     * This class will help users to easily define the badlist or the any sort of filtering in the generation process. 
     * This will help to avoid post processing filtering.
     * 
     * @param A		int[][] adjacency matrix
     * @return		boolean 
     */
    
    public boolean boundaryConditionCheck(int[][] A, String[] symbolArray) {
    	boolean check = true;
    	
    	/**
    	 * Here, users can define the functions as the boundary conditions. Example conditions
    	 * are given below. 
    	 */
    	
    	//if(detectAllenes(A, symbolArray)) check = false;
    	//if(detectAdjacentDoubleBonds(A)) check = false;
    	//if(detectTripleBonds(A)) check = false;
    	
    	return check;
    }
    
}
