/*
 MIT License

 Copyright (c) 2021-2022 Mehmet Aziz Yirik <mehmetazizyirik@outlook.com> <0000-0001-7520-7215@orcid.org>

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

package maygen;

import static org.junit.Assert.*;

import org.apache.commons.cli.*;
import org.junit.Test;

/**
 * Unit test class for the MaygenConsole class.
 *
 * @author MehmetAzizYirik mehmetazizyirik@outlook.com 0000-0001-7520-7215@orcid.org
 */
public class MaygenConsoleTest {

    @Test
    public void parseArgs() throws ParseException {
        assertFalse(
                new MaygenConsole(new Maygen())
                        .parseArgs(new String[] {"-f", "C5N3H9", "-v", "-o", "-smi"}));
        assertTrue(
                new MaygenConsole(new Maygen())
                        .parseArgs(new String[] {"-f", "C5N3H9", "-v", "-o", "-smi", "-h"}));
    }

    @Test
    public void checkSmiAndSdf() throws ParseException {
        Options options = new MaygenConsole(new Maygen()).setupOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine commandLine =
                parser.parse(options, new String[] {"-f", "C5N3H9", "-v", "-o", "-smi"});
        new MaygenConsole(new Maygen()).checkSmiAndSdf(commandLine);
        assertTrue(commandLine.hasOption("smi"));
    }

    @Test
    public void displayHelpMessage() {
        MaygenConsole maygenConsole = new MaygenConsole(new Maygen());
        Options options = maygenConsole.setupOptions();
        maygenConsole.displayHelpMessage(options);
        assertTrue(true);
    }
}
