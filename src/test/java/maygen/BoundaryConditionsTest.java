/*
 * Copyright (c) 2021 Mehmet Aziz Yirik <mehmetazizyirik@outlook.com> <0000-0001-7520-7215@orcid.org>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package maygen;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * Unit test class for the BundaryConditions class of MAYGEN.
 *
 * @author MehmetAzizYirik mehmetazizyirik@outlook.com 0000-0001-7520-7215@orcid.org
 */
public class BoundaryConditionsTest {

    @Test
    public void detectTripleBonds() {
        assertFalse(BoundaryConditions.detectTripleBonds(new int[][] {}));
    }

    @Test
    public void detectTripleBonds2() {
        assertTrue(
                BoundaryConditions.detectTripleBonds(
                        new int[][] {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}}));
    }

    @Test
    public void detectAdjacentDoubleBonds() {
        assertFalse(BoundaryConditions.detectAdjacentDoubleBonds(new int[][] {}));
    }

    @Test
    public void detectAdjacentDoubleBonds2() {
        assertTrue(
                BoundaryConditions.detectAdjacentDoubleBonds(
                        new int[][] {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}}));
    }

    @Test
    public void detectAllenes() {
        assertFalse(BoundaryConditions.detectAllenes(new int[][] {}, new String[] {}));
    }

    @Test
    public void detectAllenes2() {
        assertTrue(
                BoundaryConditions.detectAllenes(
                        new int[][] {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}},
                        new String[] {"C", "C", "C"}));
    }

    @Test
    public void boundaryConditionCheck() {
        assertTrue(BoundaryConditions.boundaryConditionCheck(new int[][] {}, new String[] {}));
    }

    @Test
    public void boundaryConditionCheck2() {
        assertFalse(
                BoundaryConditions.boundaryConditionCheck(
                        new int[][] {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}},
                        new String[] {"C", "C", "C"}));
    }
}
