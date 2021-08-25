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
 * This is the junit test class for MAYGEN. Randomly selected 40 molecular formulas are used.
 * The number of generated structures are checked. The number of isomers are also tested with
 * MOLGEN algorithm. MAYGEN generates same number of isomers like MOLGEN.
 *
 * @author Mehmet Aziz Yirik
 */
package maygen;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import org.junit.Test;

public class MaygenTest {

    @Test
    public void test_C3Cl2H4() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C3Cl2H4";
        maygen.run();
        assertEquals(7, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(7, maygen.count.get());
    }

    @Test
    public void test_O13S7() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "O13S7";
        maygen.run();
        assertEquals(1980, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(1980, maygen.count.get());
    }

    @Test
    public void test_O10S10() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "O10S10";
        maygen.run();
        assertEquals(4752, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(4752, maygen.count.get());
    }

    @Test
    public void test_S27() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "S27";
        maygen.run();
        assertEquals(1, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(1, maygen.count.get());
    }

    @Test
    public void test_O18() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "O18";
        maygen.run();
        assertEquals(1, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(1, maygen.count.get());
    }

    @Test
    public void test_C2NO2H5() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C2NO2H5";
        maygen.run();
        assertEquals(84, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(84, maygen.count.get());
    }

    @Test
    public void test_C6H6() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C6H6";
        maygen.run();
        assertEquals(217, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(217, maygen.count.get());
    }

    @Test
    public void test_C3O3H4() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C3O3H4";
        maygen.run();
        assertEquals(152, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(152, maygen.count.get());
    }

    @Test
    public void test_Cl2C5H4() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "Cl2C5H4";
        maygen.run();
        assertEquals(217, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(217, maygen.count.get());
    }

    @Test
    public void test_C5H9ClO() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5H9ClO";
        maygen.run();
        assertEquals(334, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(334, maygen.count.get());
    }

    @Test
    public void test_C6OF2H12() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C6OF2H12";
        maygen.run();
        assertEquals(536, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(536, maygen.count.get());
    }

    @Test
    public void test_C7H10() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C7H10";
        maygen.run();
        assertEquals(575, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(575, maygen.count.get());
    }

    @Test
    public void test_C6O2H12() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C6O2H12";
        maygen.run();
        assertEquals(1313, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(1313, maygen.count.get());
    }

    @Test
    public void test_F2P3BrNO2H() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "F2P3BrNO2H";
        maygen.run();
        assertEquals(1958, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(1958, maygen.count.get());
    }

    @Test
    public void test_C6OH6() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C6OH6";
        maygen.run();
        assertEquals(2237, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(2237, maygen.count.get());
    }

    @Test
    public void test_C5H6BrN() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5H6BrN";
        maygen.run();
        assertEquals(2325, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(2325, maygen.count.get());
    }

    @Test
    public void test_C6H7F2I() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C6H7F2I";
        maygen.run();
        assertEquals(3523, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(3523, maygen.count.get());
    }

    @Test
    public void test_C5F2O2H2() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5F2O2H2";
        maygen.run();
        assertEquals(7094, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(7094, maygen.count.get());
    }

    @Test
    public void test_C7OH10() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C7OH10";
        maygen.run();
        assertEquals(7166, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(7166, maygen.count.get());
    }

    @Test
    public void test_C4ClHF2O3() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C4ClHF2O3";
        maygen.run();
        assertEquals(7346, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(7346, maygen.count.get());
    }

    @Test
    public void test_C4O5H6() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C4O5H6";
        maygen.run();
        assertEquals(8070, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(8070, maygen.count.get());
    }

    @Test
    public void test_C5ClHF2O2() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5ClHF2O2";
        maygen.run();
        assertEquals(12400, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(12400, maygen.count.get());
    }

    @Test
    public void test_C5H10BrF2OP() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5H10BrF2OP";
        maygen.run();
        assertEquals(15009, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(15009, maygen.count.get());
    }

    @Test
    public void test_C9H12() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C9H12";
        maygen.run();
        assertEquals(19983, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(19983, maygen.count.get());
    }

    @Test
    public void test_C6H10O2Br2() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C6H10O2Br2";
        maygen.run();
        assertEquals(24201, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(24201, maygen.count.get());
    }

    @Test
    public void test_C10H16() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C10H16";
        maygen.run();
        assertEquals(24938, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(24938, maygen.count.get());
    }

    @Test
    public void test_C6H6ClOI() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C6H6ClOI";
        maygen.run();
        assertEquals(30728, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(30728, maygen.count.get());
    }

    @Test
    public void test_C4H5O2Br2N() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C4H5O2Br2N";
        maygen.run();
        assertEquals(41067, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(41067, maygen.count.get());
    }

    @Test
    public void test_C4H10NOSP() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C4H10NOSP";
        maygen.run();
        assertEquals(52151, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(52151, maygen.count.get());
    }

    @Test
    public void test_C7O2H10() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C7O2H10";
        maygen.run();
        assertEquals(54641, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(54641, maygen.count.get());
    }

    @Test
    public void test_P3O3NCl2() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "P3O3NCl2";
        maygen.run();
        assertEquals(665, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(665, maygen.count.get());
    }

    @Test
    public void test_C5H5SI5() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5H5SI5";
        maygen.run();
        assertEquals(2619, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(2619, maygen.count.get());
    }

    @Test
    public void test_C3O3NH5() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C3O3NH5";
        maygen.run();
        assertEquals(2644, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(2644, maygen.count.get());
    }

    @Test
    public void test_C5H9ClOS() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5H9ClOS";
        maygen.run();
        assertEquals(3763, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(3763, maygen.count.get());
    }

    @Test
    public void test_C3NO2SH7() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C3NO2SH7";
        maygen.run();
        assertEquals(3838, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(3838, maygen.count.get());
    }

    @Test
    public void test_C4H8Cl3O2P() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C4H8Cl3O2P";
        maygen.run();
        assertEquals(9313, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(9313, maygen.count.get());
    }

    @Test
    public void test_C5H2F2SO() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5H2F2SO";
        maygen.run();
        assertEquals(13446, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(13446, maygen.count.get());
    }

    @Test
    public void test_C7H11ClS() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C7H11ClS";
        maygen.run();
        assertEquals(15093, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(15093, maygen.count.get());
    }

    @Test
    public void test_C4NO3H7() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C4NO3H7";
        maygen.run();
        assertEquals(18469, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(18469, maygen.count.get());
    }

    @Test
    public void test_C4H5O2F2P() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C4H5O2F2P";
        maygen.run();
        assertEquals(41067, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(41067, maygen.count.get());
    }

    @Test
    public void test_C3N3O2H7() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C3N3O2H7";
        maygen.run();
        assertEquals(45626, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(45626, maygen.count.get());
    }

    @Test
    public void test_C5N3H9() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5N3H9";
        maygen.run();
        assertEquals(46125, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(46125, maygen.count.get());
    }

    @Test
    public void test_C3O6PH5() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C3O6PH5";
        maygen.run();
        assertEquals(51323, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(51323, maygen.count.get());
    }

    @Test
    public void test_C5H5POBr2() throws IOException {
        MAYGEN maygen = new MAYGEN();
        maygen.formula = "C5H5POBr2";
        maygen.run();
        assertEquals(62886, maygen.count.get());
        maygen.multiThread = true;
        maygen.run();
        assertEquals(62886, maygen.count.get());
    }
}
