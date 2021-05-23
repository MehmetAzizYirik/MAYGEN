package MAYGEN;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.net.URISyntaxException;

import org.junit.Test;
import org.openscience.cdk.exception.CDKException;

public class MaygenTest {
	
	/**
	 * This is the junit test class for MAYGEN. Randomly selected 40 molecular
	 * formulas are used. The number of generated structures are checked. The
	 * number of isomers are also tested with MOLGEN algorithm. MAYGEN generates 
	 * same number of isomers like MOLGEN.
	 * @throws URISyntaxException 
	 */
	
       
	@Test

	public void test_C3Cl2H4() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C3Cl2H4";
		MAYGEN.run();
		assertEquals(7,MAYGEN.count);
	}

	@Test

	public void test_C2NO2H5() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C2NO2H5";
		MAYGEN.run();
		assertEquals(84,MAYGEN.count);
	}

	@Test

	public void test_C6H6() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C6H6";
		MAYGEN.run();
		assertEquals(217,MAYGEN.count);
	}

	@Test

	public void test_C3O3H4() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C3O3H4";
		MAYGEN.run();
		assertEquals(152,MAYGEN.count);
	}

	@Test

	public void test_Cl2C5H4() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="Cl2C5H4";
		MAYGEN.run();
		assertEquals(217,MAYGEN.count);
	}

	@Test

	public void test_C5H9ClO() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5H9ClO";
		MAYGEN.run();
		assertEquals(334,MAYGEN.count);
	}

	@Test

	public void test_C6OF2H12() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C6OF2H12";
		MAYGEN.run();
		assertEquals(536,MAYGEN.count);
	}

	@Test

	public void test_C7H10() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C7H10";
		MAYGEN.run();
		assertEquals(575, MAYGEN.count);
	}

	@Test

	public void test_C6O2H12() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C6O2H12";
		MAYGEN.run();
		assertEquals(1313, MAYGEN.count);
	}

	@Test

	public void test_F2P3BrNO2H() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="F2P3BrNO2H";
		MAYGEN.run();
		assertEquals(1958,MAYGEN.count);
	}

	@Test

	public void test_C6OH6() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C6OH6";
		MAYGEN.run();
		assertEquals(2237,MAYGEN.count);
	}

	@Test

	public void test_C5H6BrN() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5H6BrN";
		MAYGEN.run();
		assertEquals(2325, MAYGEN.count);
	}

	@Test

	public void test_C6H7F2I() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C6H7F2I";
		MAYGEN.run();
		assertEquals(3523,MAYGEN.count);
	}

	@Test

	public void test_C5F2O2H2() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5F2O2H2";
		MAYGEN.run();
		assertEquals(7094, MAYGEN.count);
	}

	@Test

	public void test_C7OH10() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C7OH10";
		MAYGEN.run();
		assertEquals(7166, MAYGEN.count);
	}

	@Test

	public void test_C4ClHF2O3() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C4ClHF2O3";
		MAYGEN.run();
		assertEquals(7346,MAYGEN.count);
	}

	@Test

	public void test_C4O5H6() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C4O5H6";
		MAYGEN.run();
		assertEquals(8070, MAYGEN.count);
	}

	@Test

	public void test_C5ClHF2O2() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5ClHF2O2";
		MAYGEN.run();
		assertEquals(12400, MAYGEN.count);
	}

	@Test

	public void test_C5H10BrF2OP() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5H10BrF2OP";
		MAYGEN.run();
		assertEquals(15009, MAYGEN.count);
	}

	@Test

	public void test_C9H12() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C9H12";
		MAYGEN.run();
		assertEquals(19983, MAYGEN.count);
	}

	@Test

	public void test_C6H10O2Br2() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C6H10O2Br2";
		MAYGEN.run();
		assertEquals(24201, MAYGEN.count);
	}
	
	@Test
	
	public void test_C10H16() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C10H16";
		MAYGEN.run();
		assertEquals(24938, MAYGEN.count);
	}

	@Test

	public void test_C6H6ClOI() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C6H6ClOI";
		MAYGEN.run();
		assertEquals(30728, MAYGEN.count);
	}

	@Test

	public void test_C4H5O2Br2N() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C4H5O2Br2N";
		MAYGEN.run();
		assertEquals(41067, MAYGEN.count);
	}

	@Test

	public void test_C4H10NOSP() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C4H10NOSP";
		MAYGEN.run();
		assertEquals(52151, MAYGEN.count);
	}

	@Test

	public void test_C7O2H10() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C7O2H10";
		MAYGEN.run();
		assertEquals(54641, MAYGEN.count);
	}

	@Test

	public void test_P3O3NCl2() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="P3O3NCl2";
		MAYGEN.run();
		assertEquals(665, MAYGEN.count);
	}

	@Test

	public void test_C5H5SI5() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5H5SI5";
		MAYGEN.run();
		assertEquals(2619, MAYGEN.count);
	}

	@Test

	public void test_C3O3NH5() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C3O3NH5";
		MAYGEN.run();
		assertEquals(2644,MAYGEN.count);
	}

	@Test

	public void test_C5H9ClOS() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5H9ClOS";
		MAYGEN.run();
		assertEquals(3763, MAYGEN.count);
	}

	@Test

	public void test_C3NO2SH7() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C3NO2SH7";
		MAYGEN.run();
		assertEquals(3838, MAYGEN.count);
	}

	@Test

	public void test_C4H8Cl3O2P() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C4H8Cl3O2P";
		MAYGEN.run();
		assertEquals(9313, MAYGEN.count);
	}

	@Test

	public void test_C5H2F2SO() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5H2F2SO";
		MAYGEN.run();
		assertEquals(13446, MAYGEN.count);
	}

	@Test

	public void test_C7H11ClS() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C7H11ClS";
		MAYGEN.run();
		assertEquals(15093, MAYGEN.count);
	}

	@Test

	public void test_C4NO3H7() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C4NO3H7";
		MAYGEN.run();
		assertEquals(18469,MAYGEN.count);
	}

	@Test

	public void test_C4H5O2F2P() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C4H5O2F2P";
		MAYGEN.run();
		assertEquals(41067,MAYGEN.count);
	}

	@Test

	public void test_C3N3O2H7() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C3N3O2H7";
		MAYGEN.run();
		assertEquals(45626,MAYGEN.count);
	}

	@Test

	public void test_C5N3H9() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5N3H9";
		MAYGEN.run();
		assertEquals(46125, MAYGEN.count);
	}

	@Test

	public void test_C3O6PH5() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C3O6PH5";
		MAYGEN.run();
		assertEquals(51323,MAYGEN.count);
	}

	@Test

	public void test_C5H5POBr2() throws IOException, CDKException, CloneNotSupportedException, URISyntaxException {
		MAYGEN.formula="C5H5POBr2";
		MAYGEN.run();
		assertEquals(62886, MAYGEN.count);
	}
}

