import java.math.BigInteger;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.BitSet;
import java.util.List;
import java.util.LinkedList;

public class Factorization {

    public static void main(String[] args) {
        String numberToFactor = "130314078057199508079913"; //change here to factor other numbers
        System.out.println("Factoring " + numberToFactor);
        long startTime = System.currentTimeMillis();

        BigInteger bigIntegerToFactor = new BigInteger(numberToFactor);

        int extraEquations = 10;
        int l = 1024+extraEquations;

        List<BigInteger> factorbase = factorbaseGenerator(l - extraEquations);
        System.out.println("Factorbase generated in " +
                (System.currentTimeMillis()-startTime) + " ms");

        List<BigInteger> rSquareValues = new LinkedList();
        List<BigInteger> rValues = new LinkedList();
        List<BitSet> binaryMatrix = new LinkedList();

        rValueGenerator(rValues, rSquareValues, l, bigIntegerToFactor, factorbase, binaryMatrix);
        System.out.println("R:s generated " +
                (System.currentTimeMillis()-startTime) + " ms");
        try {
            createInputFile(binaryMatrix, binaryMatrix.size(), factorbase.size());
            System.out.println("Input file creation complete " +
                    (System.currentTimeMillis()-startTime) + " ms");
            findSolutions();
            System.out.println("Find solution complete " +
                    (System.currentTimeMillis()-startTime) + " ms");
        }
        catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }

        factorFinder(binaryMatrix, rValues, rSquareValues, bigIntegerToFactor);

        long endTime = System.currentTimeMillis();
        System.out.println("It took " + (endTime-startTime) + " ms");
        System.out.println("It took3 " + (endTime-startTime) + " ms");
    }

    private static BigInteger squareRoot(BigInteger x) {
        BigInteger right = x, left = BigInteger.ZERO, mid;
        while (right.subtract(left).compareTo(BigInteger.ONE) > 0) {
            mid = (right.add(left)).shiftRight(1);
            if (mid.multiply(mid).compareTo(x) > 0)
                right = mid;
            else
                left = mid;
        }
        return left;
    }


   
    private static void createInputFile(List<BitSet> bitSets, int M, int N)
            throws IOException {
        System.out.println("Building Binary Matrix");
        try {
            File file = new File("./input_file.txt");
            BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file));
            bufferedWriter.write(M + " " + N + "\n");
            for (int i = 0; i < M; i++) writeBitSet(bitSets.get(i), bufferedWriter, N);
            bufferedWriter.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Binary Matrix generated");
    }

    private static void writeBitSet (BitSet bitSet, BufferedWriter bufferedWriter, int length) throws Exception{
        for (int i = 0; i < length; i++) {
            if(bitSet.get(i)) bufferedWriter.write("1 ");
            else bufferedWriter.write("0 ");
        }
        bufferedWriter.write("\n");
    }

    private static void factorFinder(List<BitSet> binaryMatrix, List<BigInteger> rValues,
                                     List<BigInteger> rSquareValues, BigInteger N) {
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(
                    "output_file.txt"));
            int numberOfSolutions = Integer.parseInt(bufferedReader.readLine());
            System.out.println("Looking for factors");
            for (int i = 0; i < numberOfSolutions; i++) {
                String solution = bufferedReader.readLine().replaceAll(" ", "");
                BitSet currentSolution = toBitset(solution);

                if (checkSolution(currentSolution, N, rValues, rSquareValues, binaryMatrix.size())) break;
            }
            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static BitSet toBitset(String s) {
        BitSet t = new BitSet(s.length());
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == '1') {
                t.set(i);
            }
        }
        return t;
    }
    private static boolean checkSolution(BitSet solution, BigInteger N,
                                       List<BigInteger> rs, List<BigInteger> r2s, int rows) {
        BigInteger x2 = BigInteger.ONE;
        BigInteger y2 = BigInteger.ONE;
        for (int j = 0; j < rows; j++) {
            if (solution.get(j)) {
                x2 = x2.multiply(rs.get(j));
                y2 = y2.multiply(r2s.get(j));
            }
        }
        x2 = x2.mod(N);
        y2 = squareRoot(y2).mod(N);
        BigInteger gcd = N.gcd(y2.subtract(x2));

        if (!gcd.equals(BigInteger.ONE) && !gcd.equals(N)) {
            BigInteger p = gcd;
            BigInteger q = N.divide(p);
            
            System.out.println("\np is: " + p);
            System.out.println("q is: " + q);
            return true;
        }
        return false;
    }

    private static void findSolutions() throws IOException, InterruptedException {
        System.out.println("Preparing to run GaussBin");
        Process compileGaussBin = Runtime.
                getRuntime().
                exec("g++ -o gauss GaussBin.cpp", null,new File("./src"));
        compileGaussBin.waitFor();

        String gauss = "./src/gauss";
        ProcessBuilder processBuilder =
                new ProcessBuilder(gauss, "./input_file.txt", "./output_file.txt");

        System.out.println("Starting Gauss Solver");
        Process gaussSolver = processBuilder.start();
        gaussSolver.waitFor();
        System.out.println("Gauss Solver finished");
    }

    private static List<BigInteger> factorbaseGenerator(int wantedFactorbaseSize) {
        List<BigInteger> factorbase = new ArrayList();
        BigInteger prime = new BigInteger("2"); //ignore the number 1 in the factorbase
     boolean noRoot = true;
        while(wantedFactorbaseSize > factorbase.size()) {
            for (BigInteger integer : factorbase) {
                if (prime.mod(integer).equals(BigInteger.ZERO)) {
                    noRoot = false;
                    break;
                }
            }
            if (noRoot) {
                factorbase.add(prime);
            }
            else {
                noRoot = true; //reset noRoot var
            }
            prime = prime.add(BigInteger.ONE); //enumerate
        }
        return factorbase;
    }





    private static void rValueGenerator(List<BigInteger> rValues, List<BigInteger> r2Values, int l,
                                        BigInteger n, List<BigInteger> factorbase, List<BitSet> m){
            System.out.println("Generating R values");
            for (int k = 1; rValues.size() < l; k++) {
                BigInteger bigK = new BigInteger(Integer.toString(k));
                int jMax = 2*k;
                for (int j = 1 ;j <= jMax && rValues.size() < l; j++) {
                    BigInteger bigJ = new BigInteger(Integer.toString(j));
                    BigInteger r = squareRoot(bigK.multiply(n)).add(bigJ);
                    BigInteger rSquare = r.pow(2).mod(n);
                    BitSet binaryRow = generateBinaryRow(rSquare, factorbase);
                    if (binaryRow != null && !m.contains(binaryRow)) {
                        m.add(binaryRow);
                        rValues.add(r);
                        r2Values.add(rSquare);
                    }
                }
            }
            System.out.println("R values generated");
    }

    private static BitSet generateBinaryRow(BigInteger number, List<BigInteger> factorbase) {

        BitSet row = new BitSet(factorbase.size());
        for (int j = 0; j < factorbase.size(); j++) {
            BigInteger prime = factorbase.get(j);
            int factorFrequency = 0;
            while (number.mod(prime).equals(BigInteger.ZERO)) {
                number = number.divide(prime);
                factorFrequency++;
            }
            if (factorFrequency % 2 != 0) {
                row.set(j);
            }
        }
        if (number.compareTo(BigInteger.ONE) == 0){ //all factors in Factorbase
            return  row;
        }
        return null; //Some factor not in factorbase (too large factors)
    }


}






