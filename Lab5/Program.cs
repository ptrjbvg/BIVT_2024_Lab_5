using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }
    #region Level 1
    public long Task_1_1(int n, int k)
    {
        long answer = 0;
        static long Combinations(int n, int k)
        {
            if (n < 0 || k < 0 || k > n)
            {
                return 0;
            }
            return Factorial(n)/(Factorial(n- k) * Factorial(k));
        }
        static long Factorial(int n)
        {
            long a = 0;
            for (int i = 1; i <= n; i++)
            {
                a *= i;
            }
            return a;
        }
        answer = Combinations(n, k);
        return answer;
    }

    public int Task_1_2(double[] triangle1, double[] triangle2)
    {
        if (triangle1 == null || triangle2 == null) return -1;

        if (triangle1.Length < 3 || triangle2.Length < 3) return -1;
        double area1 = GeronArea(triangle1[0], triangle1[1], triangle1[2]);
        double area2 = GeronArea(triangle2[0], triangle2[1], triangle2[2]);
        if (area1 < 0 || area2 < 0) return -1;

        if (area1 > area2) return 1;
        if (area1 < area2) return 2;
        return 0; 
    }

    public static double GeronArea(double sideA, double sideB, double sideC)
    {
        if (sideA <= 0 || sideB <= 0 || sideC <= 0) return -1;
        if (sideA + sideB <= sideC || sideA + sideC <= sideB || sideB + sideC <= sideA) return -1;

        double semiPerimeter = (sideA + sideB + sideC) / 2;
        return Math.Sqrt(semiPerimeter * (semiPerimeter - sideA) * (semiPerimeter - sideB) * (semiPerimeter - sideC));
    }

    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        double s_first, s_second;

        s_first = GetDistance(v1, a1, time);
        s_second = GetDistance(v2, a2, time);

        if (s_first < 0 || s_second < 0)
        {
            answer = -1; 
        }
        else if (s_first > s_second)
        {
            answer = 1; 
        }
        else if (s_first < s_second)
        {
            answer = 2; 
        }
        else
        {
            answer = 0;
        }

        return answer;
    }

    public static double GetDistance(double v, double a, int t)
    {
        if (v < 0 || a < 0 || t < 0) return -1; 
        return v * t + (a * t * t) / 2; 
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int time = 1;
        double distanceFirst, distanceSecond;
        int result = 0;

        while (true)
        {
            distanceFirst = GetDistance(v1, a1, time);
            distanceSecond = GetDistance(v2, a2, time);

            if (distanceSecond >= distanceFirst)
            {
                result = time;
                break;
            }
            else
            {
                time++;
            }
        }

        return result;
    }

    #endregion

    #region Level 2
    public void Task_2_1(int[,] A, int[,] B)
    {
        int maxRowA, maxColA;
        int maxValueA = FindMax(A, out maxRowA, out maxColA);

        int maxRowB, maxColB;
        int maxValueB = FindMax(B, out maxRowB, out maxColB);

        (A[maxRowA, maxColA], B[maxRowB, maxColB]) = (B[maxRowB, maxColB], A[maxRowA, maxColA]);
    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }

    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        int diagonalMaxIndexB = FindDiagonalMaxIndex(B);
        B = DeleteRow(B, diagonalMaxIndexB);

        int diagonalMaxIndexC = FindDiagonalMaxIndex(C);
        C = DeleteRow(C, diagonalMaxIndexC);
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        // create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }

    public void Task_2_5(int[,] A, int[,] B)
    {
        // Получаем индексы строк с максимальными значениями в колонках
        int rowAIndex, rowBIndex, columnIndex = 0;
        int maxAValue = FindMaxInColumn(A, columnIndex, out rowAIndex);
        int maxBValue = FindMaxInColumn(B, columnIndex, out rowBIndex);

        // Перевернем значения между строками A и B
        int columnCount = A.GetLength(1);
        for (int col = 0; col < columnCount; col++)
        {
            // Временное значение для обмена
            int tempValue = A[rowAIndex, col];
            A[rowAIndex, col] = B[rowBIndex, col];
            B[rowBIndex, col] = tempValue;
        }
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
    }

    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here

        int highestRowPositiveCount = -1;
        int rowWithMaxPositiveCount = -1;

        for (int row = 0; row < 4; row++)
        {
            if (CountRowPositive(B, row) > highestRowPositiveCount)
            {
                highestRowPositiveCount = CountRowPositive(B, row);
                rowWithMaxPositiveCount = row;
            }
        }

        int highestColumnPositiveCount = -1;
        int columnWithMaxPositiveCount = -1;

        for (int column = 0; column < 6; column++)
        {
            if (CountColumnPositive(C, column) > highestColumnPositiveCount)
            {
                highestColumnPositiveCount = CountColumnPositive(C, column);
                columnWithMaxPositiveCount = column;
            }
        }

        int[,] updatedMatrixB = new int[5, 5];

        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (i <= rowWithMaxPositiveCount)
                    updatedMatrixB[i, j] = B[i, j];
                else if (i == rowWithMaxPositiveCount + 1)
                    updatedMatrixB[i, j] = C[j, columnWithMaxPositiveCount];
                else
                    updatedMatrixB[i, j] = B[i - 1, j];
            }
        }

        B = updatedMatrixB;
    }
    int CountRowPositive(int[,] matrix, int rowIndex)
    {
        int result = 0;
        int m = matrix.GetLength(1);
        for (int j = 0; j < m; j++)
        {
            if (matrix[rowIndex, j] > 0) result++;
        }
        return result;
    }
    int CountColumnPositive(int[,] matrix, int colIndex)
    {
        int result = 0;
        int n = matrix.GetLength(0);
        for (int i = 0; i < n; i++)
        {
            if (matrix[i, colIndex] > 0) result++;
        }
        return result;
    }

    public static int FindMax(int[,] matrix, out int row, out int column)
    {
        row = 0; column = 0;
        int maxValue = matrix[0, 0];

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > maxValue)
                {
                    maxValue = matrix[i, j];
                    row = i;
                    column = j;
                }
            }
        }
        return maxValue;
    }

    public static int FindDiagonalMaxIndex(int[,] matrix)
    {
        int maxIndex = 0;
        for (int i = 1; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, i] > matrix[maxIndex, maxIndex])
                maxIndex = i;

        }
        return maxIndex;
    }

    public static int FindMaxInColumn(int[,] matrix, int columnIndex, out int rowIndex)
    {
        rowIndex = 0;
        int maxValue = matrix[0, columnIndex];

        for (int i = 1; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, columnIndex] > maxValue)
            {
                maxValue = matrix[i, columnIndex];
                rowIndex = i;
            }
        }
        return maxValue;
    }


    public static int CountPositiveElements(int[,] matrix, int rowIndex)
    {
        int count = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[rowIndex, j] > 0)
                count++;
        }
        return count;
    }

    public static int FindColumnWithMaxPositives(int[,] matrix)
    {
        int maxPositivesColumn = 0;
        int maxCount = CountPositiveElementsInColumn(matrix, 0);

        for (int j = 1; j < matrix.GetLength(1); j++)
        {
            int count = CountPositiveElementsInColumn(matrix, j);
            if (count > maxCount)
            {
                maxCount = count;
                maxPositivesColumn = j;
            }
        }
        return maxPositivesColumn;
    }

    public static int CountPositiveElementsInColumn(int[,] matrix, int columnIndex)
    {
        int count = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, columnIndex] > 0)
                count++;
        }
        return count;
    }

    public static int[,] DeleteRow(int[,] matrix, int rowToDelete)
    {
        int totalRows = matrix.GetLength(0);
        int totalCols = matrix.GetLength(1);
        int[,] newMatrix = new int[totalRows - 1, totalCols];

        for (int i = 0, newRow = 0; i < totalRows; i++)
        {
            if (i != rowToDelete)
            {
                for (int j = 0; j < totalCols; j++)
                {
                    newMatrix[newRow, j] = matrix[i, j];
                }
                newRow++;
            }
        }
        return newMatrix;
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {

        // code here

        int[] positiveColumnSumsA = CalculatePositiveElementSums(ref A);
        int[] positiveColumnSumsC = CalculatePositiveElementSums(ref C);
        int[] answer = new int[positiveColumnSumsA.Length + positiveColumnSumsC.Length];

        for (int index = 0; index < answer.Length; index++)
        {
            if (index < A.GetLength(1)) answer[index] = positiveColumnSumsA[index];
            else answer[index] = positiveColumnSumsC[index - positiveColumnSumsA.Length];
        }

        // create and use CalculatePositiveElementSums(matrix);

        int[] CalculatePositiveElementSums(ref int[,] inputMatrix)
        {
            int[] sums = new int[inputMatrix.GetLength(1)];
            for (int column = 0; column < inputMatrix.GetLength(1); column++)
            {
                int columnSum = 0;
                for (int row = 0; row < inputMatrix.GetLength(0); row++)
                {
                    if (inputMatrix[row, column] > 0) columnSum += inputMatrix[row, column];
                }
                sums[column] = columnSum;
            }
            return sums;
        }
        // end

        return answer;
    }
    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here

        // use FindMax(matrix); from Task_2_1

        int maxRowA, maxColA, maxRowB, maxColB;
        FindMax(A, out maxRowA, out maxColA);
        FindMax(B, out maxRowB, out maxColB);
        int maxElementA, maxElementB;
        maxElementA = A[maxRowA, maxColA];
        maxElementB = B[maxRowB, maxColB];
        A[maxRowA, maxColA] = maxElementB;
        B[maxRowB, maxColB] = maxElementA;
        
    }

    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }

    public void Task_2_13(ref int[,] matrix)
    {
        // code here

        // create and use RemoveRow(matrix, rowIndex);

        // end
        int minRowIndex = -1;
        int minValue = int.MaxValue;
        int maxRowIndex = -1;
        int maxValue = int.MinValue;

        for (int rowIndex = 0; rowIndex < matrix.GetLength(0); rowIndex++)
        {
            for (int colIndex = 0; colIndex < matrix.GetLength(1); colIndex++)
            {
                if (matrix[rowIndex, colIndex] > maxValue)
                {
                    maxValue = matrix[rowIndex, colIndex];
                    maxRowIndex = rowIndex;
                }
            }
        }

        for (int rowIndex = 0; rowIndex < matrix.GetLength(0); rowIndex++)
        {
            for (int colIndex = 0; colIndex < matrix.GetLength(1); colIndex++)
            {
                if (matrix[rowIndex, colIndex] < minValue)
                {
                    minValue = matrix[rowIndex, colIndex];
                    minRowIndex = rowIndex;
                }
            }
        }

        RemoveRow(ref matrix, minRowIndex);
        if (minRowIndex < maxRowIndex)
        {
            RemoveRow(ref matrix, maxRowIndex - 1);
        }
        else if (minRowIndex > maxRowIndex)
        {
            RemoveRow(ref matrix, maxRowIndex);
        }
    }
    public void RemoveRow(ref int[,] matrix, int rowToRemove)
    {
        int[,] updatedMatrix = new int[matrix.GetLength(0) - 1, matrix.GetLength(1)];
        int currentRow = 0;
        for (int rowIndex = 0; rowIndex < matrix.GetLength(0); rowIndex++)
        {
            if (rowIndex != rowToRemove)
            {
                for (int colIndex = 0; colIndex < matrix.GetLength(1); colIndex++)
                {
                    updatedMatrix[currentRow, colIndex] = matrix[rowIndex, colIndex];
                }
                currentRow++;
            }
        }
        matrix = updatedMatrix;
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }

    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        double averageA = GetAverageWithoutMinMax(A);
        double averageB = GetAverageWithoutMinMax(B);
        double averageC = GetAverageWithoutMinMax(C);


        if (averageA < averageB && averageB < averageC)
        {
            answer = 1;
        }
        else if (averageA > averageB && averageB > averageC)
        {
            answer = -1;
        }
        else
        {
            answer = 0;
        }

        return answer;
    }
    double GetAverageWithoutMinMax(int[,] matrix)
    {
        int rowCount = matrix.GetLength(0);
        int columnCount = matrix.GetLength(1);
        int elementCount = rowCount * columnCount;

        if (elementCount <= 2) return 0;

        int minimum = int.MaxValue, maximum = int.MinValue, totalSum = 0;

        for (int row = 0; row < rowCount; row++)
        {
            for (int col = 0; col < columnCount; col++)
            {
                int currentValue = matrix[row, col];
                totalSum += currentValue;

                if (currentValue < minimum) minimum = currentValue;
                if (currentValue > maximum) maximum = currentValue;
            }
        }

        totalSum -= (minimum + maximum);
        double average = (double)totalSum / (elementCount - 2);
        return average;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        // create and use SortNegative(array);

        // end
    }

    public void Task_2_17(int[,] A, int[,] B)
    {

        A = SortRowsByMaxElements(A);
        B = SortRowsByMaxElements(B);
 
    }
    int[,] SortRowsByMaxElements(int[,] inputMatrix)
    {
        int totalRows = inputMatrix.GetLength(0);
        int totalCols = inputMatrix.GetLength(1);
        int[] maxElements = new int[totalRows];
        int[] indices = new int[totalRows];
        int[,] resultMatrix = new int[totalRows, totalCols];

        for (int i = 0; i < totalRows; i++)
        {
            int maxElement = int.MinValue;
            for (int j = 0; j < totalCols; j++)
            {
                if (inputMatrix[i, j] > maxElement)
                {
                    maxElement = inputMatrix[i, j];
                }
            }
            maxElements[i] = maxElement;
            indices[i] = i;
        }

        for (int i = 0; i < maxElements.Length; i++)
        {
            for (int j = 0; j < maxElements.Length - i - 1; j++)
            {
                if (maxElements[j] < maxElements[j + 1])
                {
                    int tempMax = maxElements[j + 1];
                    maxElements[j + 1] = maxElements[j];
                    maxElements[j] = tempMax;

                    int tempIndex = indices[j + 1];
                    indices[j + 1] = indices[j];
                    indices[j] = tempIndex;
                }
            }
        }

        for (int i = 0; i < totalRows; i++)
        {
            for (int j = 0; j < totalCols; j++)
            {
                resultMatrix[i, j] = inputMatrix[indices[i], j];
            }
        }

        for (int i = 0; i < totalRows; i++)
        {
            for (int j = 0; j < totalCols; j++)
            {
                inputMatrix[i, j] = resultMatrix[i, j];
            }
        }

        return inputMatrix;
    }


    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        int rowCount = matrix.GetLength(0);
        int colCount = matrix.GetLength(1);
        for (int i = 0; i < rowCount; i++)
        {
            bool hasZero = false;
            for (int j = 0; j < colCount; j++)
            {
                if (matrix[i, j] == 0)
                {
                    hasZero = true;
                    break;
                }
            }
            if (!hasZero) continue;
            RemoveRow(ref matrix, i);
            i--;
            rowCount--; 
        }
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);
    }
    public int[] CreateArrayFromMins(int[,] matrix)
    {
        int[] resultArray = new int[matrix.GetLength(0)];

        for (int rowIndex = 0; rowIndex < matrix.GetLength(0); rowIndex++)
        {
            resultArray[rowIndex] = matrix[rowIndex, rowIndex];
            for (int colIndex = rowIndex + 1; colIndex < matrix.GetLength(1); colIndex++)
            {
                if (matrix[rowIndex, colIndex] < resultArray[rowIndex])
                {
                    resultArray[rowIndex] = matrix[rowIndex, colIndex];
                }
            }
        }

        return resultArray;
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }

    public void Task_2_23(double[,] A, double[,] B)
    {
        MatrixValuesChange(A);
        MatrixValuesChange(B);
    }
    public void MatrixValuesChange(double[,] matrix)
    {
        int totalElements = matrix.GetLength(0) * matrix.GetLength(1);
        double[] tempArray = new double[totalElements];

        for (int row = 0; row < matrix.GetLength(0); row++)
        {
            for (int col = 0; col < matrix.GetLength(1); col++)
            {
                tempArray[row * matrix.GetLength(1) + col] = matrix[row, col];
            }
        }

        for (int i = 1; i < tempArray.Length; i++)
        {
            int currentIndex = i - 1;

            while (currentIndex >= 0 && tempArray[currentIndex] < tempArray[currentIndex + 1])
            {
                (tempArray[currentIndex], tempArray[currentIndex + 1]) = (tempArray[currentIndex + 1], tempArray[currentIndex]);
                currentIndex--;
            }
        }

        for (int row = 0; row < matrix.GetLength(0); row++)
        {
            for (int col = 0; col < matrix.GetLength(1); col++)
            {
                bool existsInTopFive = false;

                for (int topFiveIndex = 0; topFiveIndex < 5; topFiveIndex++)
                {
                    if (matrix[row, col] == tempArray[topFiveIndex])
                    {
                        existsInTopFive = true;
                        matrix[row, col] = (matrix[row, col] > 0) ? matrix[row, col] * 2 : matrix[row, col] / 2;
                        break;
                    }
                }

                if (!existsInTopFive)
                {
                    matrix[row, col] = (matrix[row, col] > 0) ? matrix[row, col] / 2 : matrix[row, col] * 2;
                }
            }
        }
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }

    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = FindRowWithMaxNegativeCount(A);
        maxB = FindRowWithMaxNegativeCount(B);
    }
    public int FindRowWithMaxNegativeCount(int[,] matrix)
    {
        int maximumCount = int.MinValue, rowIndexWithMaxCount = -1;

        for (int currentRow = 0; currentRow < matrix.GetLength(0); currentRow++)
        {
            int negativeCount = CountNegativeInRow(matrix, currentRow);

            if (negativeCount > maximumCount)
            {
                maximumCount = negativeCount;
                rowIndexWithMaxCount = currentRow;
            }
        }

        return rowIndexWithMaxCount;
    }

    public int CountNegativeInRow(int[,] matrix, int rowIndex)
    {
        int negativeCount = 0;

        for (int column = 0; column < matrix.GetLength(1); column++)
        {
            if (matrix[rowIndex, column] < 0) negativeCount++;
        }

        return negativeCount;
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }

    public void Task_2_27(int[,] A, int[,] B)
    {
        ReplaceMatrixElements(A);
        ReplaceMatrixElements(B);
    }
    public void FindRowMaxIndex(int[,] matrix, int rowIndex, out int columnIndex)
    {
        int maximumValue = matrix[rowIndex, 0];
        columnIndex = 0;

        for (int col = 0; col < matrix.GetLength(1); col++)
        {
            if (matrix[rowIndex, col] > maximumValue)
            {
                maximumValue = matrix[rowIndex, col];
                columnIndex = col;
            }
        }
    }

    public void ReplaceMaxElementOdd(int[,] matrix, int row, int column)
    {
        matrix[row, column] *= (column + 1);
    }

    public void ReplaceMaxElementEven(int[,] matrix, int row, int column)
    {
        matrix[row, column] = 0;
    }

    public void ReplaceMatrixElements(int[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int maxIndex;
            FindRowMaxIndex(matrix, i, out maxIndex);

            if ((i + 1) % 2 != 0)
            {
                ReplaceMaxElementOdd(matrix, i, maxIndex);
            }
            else
            {
                ReplaceMaxElementEven(matrix, i, maxIndex);
            }
        }
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion

    #region Level 3
    public delegate double SumFunction(double x, int i, ref int iCalculated);
    public delegate double YFunction(double x);

    public double SumFunction1(double x, int i, ref int iFactorial)
    {
        if (i != 0) iFactorial *= i;

        return Math.Cos(i * x) / (iFactorial * (i == 0 ? 1 : 1));
    }

    public double SumFunction2(double x, int i, ref int iSign)
    {
        iSign = -iSign;

        return iSign * Math.Cos(i * x) / Math.Pow(i, 2);
    }

    public double YFunction1(double x)
    {
        return Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
    }

    public double YFunction2(double x)
    {
        return (x * x - Math.PI * Math.PI / 3) / 4;
    }

    public double CalculateSum(SumFunction sumFunction, double x, int i)
    {
        int iCalculated = 1;
        double totalSum = 0;
        double term = sumFunction(x, i++, ref iCalculated);

        while (Math.Abs(term) > 1e-4)
        {
            totalSum += term;
            term = sumFunction(x, i++, ref iCalculated);
        }

        return totalSum;
    }

    public void GetSumAndY(double[,] sumAndY, SumFunction sumFunction, YFunction yFunction, double a, double h, int firstI)
    {
        for (int index = 0; index < sumAndY.GetLength(0); index++)
        {
            double x = a + h * index;

            sumAndY[index, 0] = CalculateSum(sumFunction, x, firstI);
            sumAndY[index, 1] = yFunction(x);
        }
    }
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        double start1 = 0.1, end1 = 1, step1 = 0.1;
        firstSumAndY = new double[(int)((end1 - start1) / step1) + 1, 2];
        GetSumAndY(firstSumAndY, SumFunction1, YFunction1, start1, step1, 0);

        double start2 = Math.PI / 5, end2 = Math.PI, step2 = Math.PI / 25;
        secondSumAndY = new double[(int)((end2 - start2) / step2) + 1, 2];
        GetSumAndY(secondSumAndY, SumFunction2, YFunction2, start2, step2, 1);
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }

    public double Task_3_3(double[] array)
    {
        SwapDirector swapper;

        double average = 0;

        foreach (double element in array)
        {
            average += element;
        }

        average /= array.Length;

        if (array[0] > average) swapper = SwapLeft;
        else swapper = SwapRight;

        swapper(array);

        return GetSum(array, 1, 2);
    }
    public delegate void SwapDirector(double[] array);

    public void SwapRight(double[] array)  // start swapping from the last element
    {
        for (int index = array.Length - 1; index > 0; index -= 2)
        {
            (array[index], array[index - 1]) = (array[index - 1], array[index]);
        }
    }

    public void SwapLeft(double[] array)  // start swapping from the first element
    {
        for (int index = 0; index < array.Length - 1; index += 2)
        {
            (array[index], array[index + 1]) = (array[index + 1], array[index]);
        }
    }

    public double GetSum(double[] array, int start, int h)
    {
        double total = 0;

        for (int i = start; i < array.Length; i += h)
        {
            total += array[i];
        }

        return total;
    }
    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }

    public void Task_3_5(out int func1, out int func2)
    {
        func1 = CountSignFlips(FirstFunction, 0, 2, 0.1);
        func2 = CountSignFlips(SecondFunction, -1, 1, 0.2);

        double FirstFunction(double x)
        {
            return Math.Pow(x, 2) - Math.Sin(x);
        }

        double SecondFunction(double x)
        {
            return Math.Exp(x) - 1;
        }
    }

    public int CountSignFlips(YFunction yFunction, double a, double b, double h)
    {
        int size = (int)Math.Ceiling((b - a) / h) + 1;
        int index = 0, flipCount = 0;
        double[] results = new double[size];

        for (double x = a; x <= b; x += h)
        {
            results[index++] = yFunction(x);
        }

        double prevValue = results[0];
        Console.WriteLine(string.Join(", ", results));

        for (int i = 1; i < size; i++)
        {
            if (results[i] == 0) continue;

            if ((prevValue > 0 && results[i] < 0) || (prevValue < 0 && results[i] > 0))
            {
                flipCount++;
            }

            prevValue = results[i];
        }

        return flipCount;
    }


    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }

    public void Task_3_7(ref int[,] B, int[,] C)
    {
        CountPositive operation = CountRowPositive;

        int positiveCount;
        int highestRowPositives = int.MinValue, highestRowIndex = -1;
        int highestColPositives = int.MinValue, highestColIndex = -1;

        for (int row = 0; row < B.GetLength(0); row++)
        {
            positiveCount = operation(B, row);
            if (positiveCount > highestRowPositives)
            {
                highestRowPositives = positiveCount;
                highestRowIndex = row;
            }
        }

        operation = CountColumnPositive;
        for (int col = 0; col < C.GetLength(1); col++)
        {
            positiveCount = operation(C, col);
            if (positiveCount > highestColPositives)
            {
                highestColPositives = positiveCount;
                highestColIndex = col;
            }
        }

        B = InsertColumn(B, highestRowIndex, C, highestColIndex);
    }
    public delegate int CountPositive(int[,] matrix, int index);

    public int[,] InsertColumn(int[,] matrixB, int insertAfterRow, int[,] matrixC, int columnIndexC)
    {
        int rowCountB = matrixB.GetLength(0), colCountB = matrixB.GetLength(1), rowCountC = matrixC.GetLength(0);
        int[] columnToInsert = new int[rowCountC];

        for (int i = 0; i < rowCountC; i++)
        {
            columnToInsert[i] = matrixC[i, columnIndexC];
        }

        int[,] augmentedB = new int[rowCountB + 1, colCountB];

        for (int i = 0; i < (rowCountB + 1); i++)
        {
            if (i <= insertAfterRow)
            {
                for (int j = 0; j < colCountB; j++)
                {
                    augmentedB[i, j] = matrixB[i, j];
                }
            }
            else if (insertAfterRow + 1 == i)
            {
                for (int j = 0; j < colCountB; j++)
                {
                    augmentedB[i, j] = columnToInsert[j];
                }
            }
            else
            {
                for (int j = 0; j < colCountB; j++)
                {
                    augmentedB[i, j] = matrixB[i - 1, j];
                }
            }
        }

        return augmentedB;
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here

        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)

        // end
    }

    public void Task_3_13(ref int[,] matrix)
    {
        RemoveRows(ref matrix, FindMaxIndex, FindMinIndex);
    }
    public void FindMaxIndex(int[,] matrix, out int maxI, out int maxJ)
    {
        int maximumValue = matrix[0, 0];
        maxI = 0;
        maxJ = 0;

        for (int row = 0; row < matrix.GetLength(0); row++)
        {
            for (int col = 0; col < matrix.GetLength(1); col++)
            {
                if (matrix[row, col] > maximumValue)
                {
                    maximumValue = matrix[row, col];
                    maxI = row;
                    maxJ = col;
                }
            }
        }
    }

    public void FindMinIndex(int[,] matrix, out int indexI, out int indexJ)
    {
        int minimumValue = matrix[0, 0];
        indexI = 0;
        indexJ = 0;

        for (int row = 0; row < matrix.GetLength(0); row++)
        {
            for (int col = 0; col < matrix.GetLength(1); col++)
            {
                if (matrix[row, col] < minimumValue)
                {
                    minimumValue = matrix[row, col];
                    indexI = row;
                    indexJ = col;
                }
            }
        }
    }

    public delegate void FindElementDelegate(int[,] matrix, out int indexI, out int indexJ);

    public void RemoveRows(ref int[,] matrix, FindElementDelegate findMax, FindElementDelegate findMin)
    {
        int maxRowIndex, maxColIndex, minRowIndex, minColIndex;
        findMax(matrix, out maxRowIndex, out maxColIndex);
        findMin(matrix, out minRowIndex, out minColIndex);

        if (minRowIndex < maxRowIndex)
        {
            RemoveRow(ref matrix, maxRowIndex);
            RemoveRow(ref matrix, minRowIndex);
        }
        else if (minRowIndex > maxRowIndex)
        {
            RemoveRow(ref matrix, minRowIndex);
            RemoveRow(ref matrix, maxRowIndex);
        }
        else
        {
            RemoveRow(ref matrix, minRowIndex);
        }
    }

    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }

    public void Task_3_27(int[,] A, int[,] B)
    {
        EvenOddRowsTransform(A, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        EvenOddRowsTransform(B, ReplaceMaxElementOdd, ReplaceMaxElementEven);
    }
    public delegate void ReplaceMaxElement(int[,] matrix, int rowIndex, int max);

    public void EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement replaceMaxElementOdd,
        ReplaceMaxElement replaceMaxElementEven)
    {
        for (int row = 0; row < matrix.GetLength(0); row++)
        {
            int maxIndex;
            if ((row + 1) % 2 != 0)
            {
                FindRowMaxIndex(matrix, row, out maxIndex);
                replaceMaxElementOdd(matrix, row, maxIndex);
            }
            else
            {
                FindRowMaxIndex(matrix, row, out maxIndex);
                replaceMaxElementEven(matrix, row, maxIndex);
            }
        }
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion
    #region bonus part
    public double[,] Task_4(double[,] matrix, int index)
    {
        // MatrixConverter[] mc = new MatrixConverter[4]; - uncomment me

        // code here

        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end

        return matrix;
    }
    #endregion
}
