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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.junit.Test;

/**
 * Unit test class for the MaygenCLI class.
 *
 * @author MehmetAzizYirik mehmetazizyirik@outlook.com 0000-0001-7520-7215@orcid.org
 */
public class MaygenCLITest {

    @Test
    public void parseArgs() throws ParseException {
        assertFalse(
                new MaygenCLI()
                        .parseArgs(new String[] {"-f", "C5N3H9", "-v", "-o", "-smi"}));
        assertFalse(
                new MaygenCLI()
                        .parseArgs(new String[] {"-f", "C5N3H9", "-v", "-o", "-smi", "-sdf"}));
        assertTrue(
                new MaygenCLI()
                        .parseArgs(new String[] {"-f", "C5N3H9", "-v", "-o", "-smi", "-h"}));
    }

    @Test
    public void displayHelpMessage() {
        MaygenCLI maygenConsole = new MaygenCLI();
        Options options = maygenConsole.setupOptions();
        maygenConsole.displayHelpMessage(options);
        assertTrue(true);
    }

    @Test
    public void main() {
        MaygenCLI.main(new String[] {"-f", "C5N3H9", "-v", "-o", "-smi", "-h"});
        assertTrue(true);
    }
}
