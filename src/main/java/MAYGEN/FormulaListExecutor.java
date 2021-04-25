package MAYGEN;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.Instant;
import java.util.List;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.joining;
import static java.util.stream.Collectors.toList;

public class FormulaListExecutor {
    private String formulaList;

    private void parseArgs(String[] args) throws ParseException {
        Options options = setupOptions(args);
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);
            this.formulaList = cmd.getOptionValue("formula-list");
        } catch (ParseException e) {
            displayHelp(options);
            System.out.println("Problem parsing command line");
        }
    }

    private void displayHelp(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.setOptionComparator(null);
        String header = "\nExecutes generator for the formulas list."
                + " The input is component strings."
                + "For example '-fl formula-list.txt'.\n\n";
        String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/MAYGEN";
        formatter.printHelp("java -cp MAYGEN.jar MAYGEN.FormulaListExecutor", header, options, footer, true);
    }

    private Options setupOptions(String[] args) {
        Options options = new Options();
        Option fl = Option.builder("fl")
                .required(true)
                .hasArg()
                .longOpt("formula-list")
                .desc("Formula list")
                .build();
        options.addOption(fl);
        return options;
    }

    public static void main(String[] args) throws IOException, ParseException {
        FormulaListExecutor gen = new FormulaListExecutor();
        gen.parseArgs(args);
        gen.executeFormulas();
    }

    private void executeFormulas() throws IOException {
        class Line {
            String formula;
            double averageTime;
            long numberOfStructure;

            public Line(String formula, double averageTime, long numberOfStructure) {
                this.formula = formula;
                this.averageTime = averageTime;
                this.numberOfStructure = numberOfStructure;
            }

            @Override
            public String toString() {
                return "Line{" +
                        "formula='" + formula + '\'' +
                        ", averageTime=" + averageTime +
                        ", numberOfStructure=" + numberOfStructure +
                        '}';
            }
        }
        List<Line> lines = Files.lines(Paths.get(this.formulaList)).map(formula -> {
            double average = IntStream.rangeClosed(1, 5).asLongStream().map(
                    i -> execMorgen(formula)).average().getAsDouble();
            return new Line(formula, average, MAYGEN.count);
        }).collect(toList());
        Files.write(Paths.get(this.formulaList + ".csv"),
                lines.stream().map(line ->
                        String.join("|", MAYGEN.normalizeFormula(line.formula),
                                String.valueOf(line.averageTime), String.valueOf(line.numberOfStructure)))
                        .collect(joining("\n")).getBytes(StandardCharsets.UTF_8)
        );
        System.out.println(lines.size() + " formulas have been executed.");
    }

    private static long execMorgen(String formula) {
        Instant start = Instant.now();
        MAYGEN.count = 0;
        try {
            MAYGEN.formula = formula;
            MAYGEN.run();
        } catch (Exception e) {
            e.printStackTrace();
        }
        Instant finish = Instant.now();
        return Duration.between(start, finish).toMillis();
    }
}
