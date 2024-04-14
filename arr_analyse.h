#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define N 3 // Size of the square matrix defined for finding eigen values and gaussian elimination

#define ALLOC(p, n)                                         \
    do                                                      \
    {                                                       \
        if (!((p) = calloc((n), sizeof(*(p)))))             \
        {                                                   \
            fprintf(stderr, "Memory allocation failure\n"); \
            exit(1);                                        \
        }                                                   \
    } while (0)

#define PI 3.14159265
#define float double

int size(int *arr)
{
    int out = sizeof(arr) / sizeof(arr[0]);
}

void gaussian_elimination(double A[N][N+1]) {
    for (int i = 0; i < N; i++) {
        // Find pivot for this column
        double max_val = 0.0;
        int max_row = i;
        for (int k = i; k < N; k++) {
            if (abs(A[k][i]) > max_val) {
                max_val = abs(A[k][i]);
                max_row = k;
            }
        }

        // Swap maximum row with current row (if needed)
        if (max_row != i) {
            for (int k = i; k < N+1; k++) {
                double temp = A[i][k];
                A[i][k] = A[max_row][k];
                A[max_row][k] = temp;
            }
        }

        // Make all rows below this one 0 in current column
        for (int k = i+1; k < N; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < N+1; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }
}

// Function to back-substitute to find the solution
void back_substitute(double A[N][N+1], double x[N]) {
    for (int i = N-1; i >= 0; i--) {
        x[i] = A[i][N];
        for (int j = i+1; j < N; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

// Function to solve a system of linear equations
void solve(double A[N][N], double b[N], double x[N]) {
    // Augment matrix A with vector b
    double augmented[N][N+1];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][N] = b[i];
    }

    // Perform Gaussian elimination
    gaussian_elimination(augmented);

    // Back-substitute to find the solution
    back_substitute(augmented, x);
}

double dot_product(double *a, double *b, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += a[i] * b[i];
    }
    return result;
}

// Function to normalize a vector
void normalize(double *v, int n) {
    double norm = sqrt(dot_product(v, v, n));
    for (int i = 0; i < n; i++) {
        v[i] /= norm;
    }
}

// Function to calculate the eigenvalues and eigenvectors of a matrix
void eig(double A[N][N], double eigenvalues[N], double eigenvectors[N][N]) {
    // Initialize eigenvectors matrix as identity matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Power iteration method
    int max_iterations = 1000;
    double epsilon = 1e-6;
    for (int k = 0; k < N; k++) {
        double x[N] = {1.0}; // Initial guess for eigenvector
        double lambda_old = 0.0;

        for (int iter = 0; iter < max_iterations; iter++) {
            double y[N];
            for (int i = 0; i < N; i++) {
                y[i] = 0.0;
                for (int j = 0; j < N; j++) {
                    y[i] += A[i][j] * x[j];
                }
            }

            // Calculate eigenvalue
            double lambda = dot_product(y, x, N);
            if (fabs(lambda - lambda_old) < epsilon) {
                eigenvalues[k] = lambda;
                break;
            }
            lambda_old = lambda;

            // Normalize eigenvector
            normalize(y, N);

            // Update eigenvector
            for (int i = 0; i < N; i++) {
                x[i] = y[i];
                eigenvectors[i][k] = x[i];
            }
        }
    }
}

void matrixmultiply(int n,int c,int m,int a[n][c],int b[c][m])
{
  int ans[n][m];
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<m;j++)
    {
      int sum = 0;
      for(int k=0;k<c;k++)
      {
        sum += a[i][k]*b[k][j];
      }
      ans[i][j] = sum;
    }
  }
  printf("the multiplication matrix is : \n");
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<m;j++)
    {
      printf("%d ",ans[i][j]);
    }
    printf("\n");
  }

}
void diagflat(int n,int a[],int offset)
{
  int x = offset;
  if(x<0)
  {
    x = -x;
  }
  int m = n+x;
  int ans[m][m];
  for(int i=0;i<m;i++)
  {
    for(int j=0;j<m;j++)
    {
      ans[i][j] = 0;
    }
  }
  if(offset>=0)
  {
    for(int i=0;i<m-offset;i++)
    {
      ans[i][i+offset] = a[i];
    }
  }
  else
  {
    int pnt = 0;
    for(int i=x;i<m;i++)
    {
      ans[i][i+offset] = a[pnt];
      pnt++;
    }
  }
  printf("the diagflat matrix is : \n");
  for(int i=0;i<m;i++)
  {
    for(int j=0;j<m;j++)
    {
      printf("%d ",ans[i][j]);
    }
    printf("\n");
  }

}

int value_of_poly(int *arr, int degree, int x)
{
    int n = degree;
    int value = 0;

    int i = 0, j = 0, temp;
    while (i < (n + 1))
    {
        temp = arr[i];
        while (j <= n - i)
        {
            temp *= x;
            j++;
        }
        value += temp;
        i++;
    }
    return value;
}


void trace(long long n,long long a[][n])
{
    long long trace = 0;
    for(long long i=0;i<n;i++)
    {
        trace += a[i][i];
    }
    printf("the trace of the martix is %ld\n",trace);
}

void triu(long long n,long long a[][n])
{
    printf("the upper triangular matrix is :\n");
    
    long long triu[n][n];
    for(long long i=0;i<n;i++)
    {
        for(long long j=0;j<n;j++)
        {
          if(i>=j)
          {
            triu[i][j] = a[i][j];
          }
          else
          {
            triu[i][j] = 0;
          }
        }
    }

    for(long long i=0;i<n;i++)
    {
        for(long long j=0;j<n;j++)
        {
            printf("%ld ",triu[i][j]);
        }
        printf("\n");
    }
}

void tril(long long n,long long a[][n])
{
    printf("the lower triangular matrix is :\n");
    
    long long tril[n][n];
    int i =0, j =0;
    while(i<n)
    {
        while(j<n)
        {
          if(i>=j)
          {
            tril[i][j] = 0;
          }
          else
          {
            tril[i][j] = a[i][j];
          }
            j++;
        }
        i++;
    }

    i = 0;
    j = 0;
    while(i<n)
    {
        while(j<n)
        {
            printf("%ld ",tril[i][j]);
            j++;
        }
        printf("\n");
        i++;
    }
  
}

// To create a function to calculate the length of an array
int length(int arr[])
{
    int length = 0;
    int i = 0;
    while (arr[i] != '\0')
    {
        length++;
        i++;
    }
    return length;
}

// To create a function that tells the dimension of an array
int dimension(void *arr)
{
    int dim = 1;
    int size = sizeof(int);
    while (1)
    {
        if ((void *)arr != NULL)
        {
            dim++;
            arr = ((void **)arr)[1]; // increasing the dimension of the array
        }
        else
            break;
    }
    if (dim >= 1)
        return dim;
    else
        return -1;
}

// Bubble sort function
void bubbleSort(int arr[], int n)
{
    int i = 0;
    while (i < n - 1)
    {
        int j = 0;
        while (j < n - i - 1)
        {
            if (arr[j] > arr[j + 1])
            {
                // swap arr[j] and arr[j+1]
                int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
            j++;
        }
        i++;
    }
}

// Insertion sort function
void insertionSort(int arr[], int n)
{
    int i = 1;
    while (i < n)
    {
        int key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
        i++;
    }
}

// Selection sort function
void selectionSort(int arr[], int n)
{
    int i = 0;
    while (i < n - 1)
    {
        int minIndex = i;
        int j = i + 1;
        while (j < n)
        {
            if (arr[j] < arr[minIndex])
            {
                minIndex = j;
            }
            j++;
        }
        // swap arr[i] and arr[minIndex]
        int temp = arr[i];
        arr[i] = arr[minIndex];
        arr[minIndex] = temp;
        i++;
    }
}

// Function to input arrays
void exampleArrays(int numArrays)
{
    int i, j, n;

    for (i = 0; i < numArrays; i++)
    {
        printf("Enter the size of array %d: ", i + 1);
        scanf("%d", &n);

        int arr[n];
        printf("Enter %d elements for array %d: ", n, i + 1);
        for (j = 0; j < n; j++)
        {
            scanf("%d", &arr[j]);
        }

        printf("Array %d: ", i + 1);
        for (j = 0; j < n; j++)
        {
            printf("%d ", arr[j]);
        }
        printf("\n");
    }
}

// To create a function to determine the shape of an array
int shape(int *arr)
{
    int dim_size = size(arr);
    int shape = dim_size / sizeof(int);
    return shape;
}

// To return an equally spaced array
double *linspace(double start, double end, int n)
{
    double *x;
    int i = 0;
    double step = (end - start) / ((double)n - 1); // the space between to consecutive elements
    x = (double *)calloc(n, sizeof(double));       // allocation of memory and initialisation woth 0
    x[i] = start;

    i = 1;
    while (i < n)
    {
        x[i] = x[i - 1] + step;
        i++;
    }
    x[n - 1] = end;
    return x;
}

// to return an array equally spaced on the logarithmic scale
double *logspace(double start, double end, int n, double base)
{
    double *result;
    result = (double *)calloc(n, sizeof(double));
    double log_start = log10(start) / log10(base);
    double log_end = log10(end) / log10(base);
    double step = (log_end - log_start) / (n - 1); // calculating the difference between 2 numbers ont he logarithmic scale

    int i = 0;
    while (i < n)
    {
        result[i] = pow(base, log_start + i * step);
        i++;
    }
    return result;
}

// Values are generated within the half-open interval [start, stop), with spacing between values given by step.
double *arange(double start, double end, double step)
{
    int n = (int)((end - start) / step) + 1;
    double *result = (double *)calloc(n, sizeof(double));

    for (int i = 0; i < n; i++)
    {
        result[i] = start + i * step;
    }

    return result;
}

// To calculate mean of an array
int mean(int *arr)
{
    int n = length(arr); // Number of data points
    double sum = 0.0;
    double mean;

    int i = 0;
    while (i < n)
    {
        sum += arr[i];
        i++;
    }
    mean = sum / n;
    return mean;
}

typedef struct
{
    int row;
    int col;
} Index;

Index *argwhere(int **indicesArray, int rows, int cols, int reqCondition, int *length)
{

    Index *indices = (Index *)malloc(rows * cols * sizeof(Index));
    if (indices == NULL)
    {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    int count = 0;

    // Loop through the array to find indices where the value is greater than the threshold
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (indicesArray[i][j] > reqCondition)
            {
                indices[count].row = i;
                indices[count].col = j;
                count++;
            }
        }
    }

    *length = count;

    return indices;
}

void pad_array(int *array, int rows, int cols, int pad_top, int pad_bottom, int pad_left, int pad_right)
{
    int new_rows = rows + pad_top + pad_bottom;
    int new_cols = cols + pad_left + pad_right;

    int *padded_array = (int *)malloc(new_rows * new_cols * sizeof(int));

    // Fill padded_array with original array
    for (int i = 0; i < new_rows; i++)
    {
        for (int j = 0; j < new_cols; j++)
        {
            if (i >= pad_top && i < rows + pad_top && j >= pad_left && j < cols + pad_left)
            {
                padded_array[i * new_cols + j] = array[(i - pad_top) * cols + (j - pad_left)];
            }
            else
            {
                padded_array[i * new_cols + j] = 0; // Assuming padding with zeros
            }
        }
    }

    // Copy padded_array back to array
    for (int i = 0; i < new_rows; i++)
    {
        for (int j = 0; j < new_cols; j++)
        {
            array[i * new_cols + j] = padded_array[i * new_cols + j];
        }
    }

    free(padded_array);
}

// Function to convert degrees to radians
double toRadians(double degree)
{
    return degree * (PI / 180);
}

// Function to calculate factorial
double factorial(int n)
{
    double result = 1;
    for (int i = 2; i <= n; ++i)
    {
        result *= i;
    }
    return result;
}

// Function to calculate sine using Taylor series
double sine(double angle)
{
    angle = toRadians(angle);
    double result = 0;
    double power = angle;
    double sign = 1;

    for (int i = 1; i <= 10; ++i)
    {
        result += sign * power / factorial(2 * i - 1);
        power *= angle * angle;
        sign = -sign;
    }

    return result;
}

// Function to calculate cosine using Taylor series
double cosine(double angle)
{
    angle = toRadians(angle);
    double result = 0;
    double power = 1;
    double sign = 1;

    for (int i = 0; i <= 10; ++i)
    {
        result += sign * power / factorial(2 * i);
        power *= angle * angle;
        sign = -sign;
    }

    return result;
}

// Function to calculate tangent using Taylor series
double tangent(double angle)
{
    angle = toRadians(angle);
    double result1 = 0;
    double power1 = angle;
    double sign1 = 1;

    for (int i = 1; i <= 10; ++i)
    {
        result1 += sign1 * power1 / factorial(2 * i - 1);
        power1 *= angle * angle;
        sign1 = -sign1;
    }
    angle = toRadians(angle);
    double result2 = 0;
    double power2 = 1;
    double sign2 = 1;

    for (int i = 0; i <= 10; ++i)
    {
        result2 += sign2 * power2 / factorial(2 * i);
        power2 *= angle * angle;
        sign2 = -sign2;
    }
    double result = 0;
    result = result1 / result2;
    return result;
}

// Function to calculate cotangent using Taylor series
double cotangent(double angle)
{
    angle = toRadians(angle);
    double result1 = 0;
    double power1 = angle;
    double sign1 = 1;

    for (int i = 1; i <= 10; ++i)
    {
        result1 += sign1 * power1 / factorial(2 * i - 1);
        power1 *= angle * angle;
        sign1 = -sign1;
    }
    angle = toRadians(angle);
    double result2 = 0;
    double power2 = 1;
    double sign2 = 1;

    for (int i = 0; i <= 10; ++i)
    {
        result2 += sign2 * power2 / factorial(2 * i);
        power2 *= angle * angle;
        sign2 = -sign2;
    }
    double result = 0;
    result = result2 / result1;
    return result;
}

// Function to calculate cosecant using Taylor series
double cosecant(double angle)
{
    angle = toRadians(angle);
    double result = 0;
    double power = angle;
    double sign = 1;

    for (int i = 1; i <= 10; ++i)
    {
        result += sign * power / factorial(2 * i - 1);
        power *= angle * angle;
        sign = -sign;
    }
    result = 1 / result;
    return result;
}

// Function to calculate secant using Taylor series
double secant(double angle)
{
    angle = toRadians(angle);
    double result = 0;
    double power = 1;
    double sign = 1;

    for (int i = 0; i <= 10; ++i)
    {
        result += sign * power / factorial(2 * i);
        power *= angle * angle;
        sign = -sign;
    }
    result = 1 / result;
    return result;
}

void polynomial(int *arr, int length)
{
    int deg = length - 1;
    int n = 0;
    int j = 0;
    while (j < length)
    {
        if (arr[j] == 0)
            break;

        if (arr[j] != 0)
        {
            n++;
        }
        j++;
    }

    int terms = n; // No. of terms
    int *exp = (int *)malloc(terms * sizeof(int));
    int exponent[terms];

    int k = 0;
    while (k < terms)
    {
        exponent[k] = deg - k;
        k++;
    }

    int i = 0;
    while (i < length)
    {
        if (arr[i] != 0)
        {
            printf("%dx^%d", arr[i], exponent[i]);

            if (i < n - 1)
            {
                printf(" + ");
            }
        }
        i++;
    }
    free(exp);
}

void integrate(int *arr, int length)
{
    int deg = length - 1;
    int n = 0;
    int j = 0;
    while (j < length)
    {
        if (arr[j] == 0)
            break;

        if (arr[j] != 0)
        {
            n++;
        }
        j++;
    }

    int terms = n; // No. of terms
    int *exp = (int *)malloc(terms * sizeof(int));
    int exponent[terms];

    int k = 0;
    while (k < terms)
    {
        exponent[k] = deg - k;
        k++;
    }

    int i = 0;
    int expt, coeff;
    while (i < length)
    {
        if (arr[i] != 0)
        {
            coeff = arr[i] / (exponent[i] + 1);
            expt = exponent[i] + 1;

            if (coeff != 0)
            {
                printf("%dx^%d", coeff, expt);
            }

            if (i < n - 1)
            {
                printf(" + ");
            }
        }
        i++;
    }
    free(exp);
}

void differentiate(int *arr, int length)
{
    int deg = length - 1;
    int n = 0;
    int j = 0;
    while (j < length)
    {
        if (arr[j] == 0)
            break;

        if (arr[j] != 0)
        {
            n++;
        }
        j++;
    }

    int terms = n; // No. of terms
    int *exp = (int *)malloc(terms * sizeof(int));
    int exponent[terms];

    int k = 0;
    while (k < terms)
    {
        exponent[k] = deg - k;
        k++;
    }

    int i = 0;
    int expt, coeff;
    while (i < length)
    {
        if (arr[i] != 0)
        {
            coeff = arr[i] * exponent[i];
            expt = exponent[i] - 1;

            if (coeff != 0)
            {
                printf("%dx^%d", coeff, expt);
            }

            if (i < n - 1)
            {
                printf(" + ");
            }
        }
        i++;
    }
    free(exp);
}

void *reshape_2d_3d(size_t id1, size_t id2, int iar[][id2],
                    size_t od1, size_t od2, size_t od3)
{
    // oar is a pointer to a multidimensional array; in this case, it will point to the first element of an array of arrays (of arrays).
    int(*oar)[od2][od3];
    size_t size1 = id1 * id2;
    size_t size2 = od1 * od2 * od3;
    size_t min_size = (size1 <= size2) ? size1 : size2;

    ALLOC(oar, od1);

    for (size_t i = 0; i < min_size; i++)
    {
        oar[i / (od2 * od3)][(i / od3) % od2][i % od3] = iar[i / id2][i % id2];
    }
    return oar;
}

void *reshape_1d_2d(size_t id1, int *iar,
                    size_t od1, size_t od2)
{
    // oar is a pointer to a multidimensional array; in this case, it will point to the first element of an array of arrays (of arrays).
    int(*oar)[od2];
    size_t size1 = id1;
    size_t size2 = od1 * od2;
    size_t min_size = (size1 <= size2) ? size1 : size2;

    ALLOC(oar, od1);

    for (size_t i = 0; i < min_size; i++)
    {
        oar[i / (od2)][i % od2] = iar[i];
    }
    return oar;
}

void *reshape_1d_3d(size_t id1, int *iar,
                    size_t od1, size_t od2, size_t od3)
{
    // oar is a pointer to a multidimensional array; in this case, it will point to the first element of an array of arrays (of arrays).
    int(*oar)[od2][od3];
    size_t size1 = id1;
    size_t size2 = od1 * od2 * od3;
    size_t min_size = (size1 <= size2) ? size1 : size2;

    ALLOC(oar, od1);

    for (size_t i = 0; i < min_size; i++)
    {
        oar[i / (od2 * od3)][(i / od3) % od2][i % od3] = iar[i];
    }
    return oar;
}

void regression(int n, int ax[], int ay[])
{
    int i;
    float x, y, m, c, d;
    float sumx = 0, sumxsq = 0, sumy = 0, sumxy = 0;
    float sumysq = 0;
    float mae = 0;
    float mse = 0;
    float p = n;
    for (int i = 0; i < n; i++)
    {
        sumx = sumx + x;
        sumxsq = sumxsq + (x * x);
        sumy = sumy + y;
        sumxy = sumxy + (x * y);
        sumysq = sumysq + (y * y);
        ax[i] = x;
        ay[i] = y;
    }
    d = p * sumxsq - sumx * sumx;
    m = (p * sumxy - sumx * sumy) / d;
    c = (sumy * sumxsq - sumx * sumxy) / d;
    while (1)
    {
        printf("1.compute sum\n");
        printf("2.find best fitting line\n");
        printf("3.end\n");
        int num;
        scanf("%d", &num);
        if (num == 1)
        {
            printf("sum x = %lf\n", sumx);
            printf("sum y = %lf\n", sumy);
            printf("sum xy = %lf\n", sumxy);
            printf("sum xsq = %lf\n", sumxsq);
            printf("sum ysq = %lf\n", sumysq);
        }
        else if (num == 2)
        {
            float p = n;
            printf("the best fitting line is y = %lfx+%lf\n", m, c);
            for (int i = 0; i < n; i++)
            {
                mae += abs(ay[i] - (m * ax[i] + c));
                mse += (ay[i] - (m * ax[i] + c)) * (ay[i] - (m * ax[i] + c));
            }
            printf("mean absolute error = %lf\n", mae / p);
            printf("mean squared error = %lf\n", mse / p);
        }
        else if (num == 3)
        {
            printf("END");
        }
    }
}

// Function to calculate the square root of a number
double sqrt(double x)
{
    double guess = x / 2.0;
    double prevGuess;
    do
    {
        prevGuess = guess;
        guess = (guess + x / guess) / 2.0;
    } while (prevGuess - guess > 0.00001); // Change the threshold for desired precision
    return guess;
}

// Function to calculate the inverse trigonometric arcsine (asin)
double arcsin(double x)
{
    if (x < -1.0 || x > 1.0)
        return NAN; // Not a Number for invalid input

    double result = 0.0;
    double term = x;
    double xSquared = x * x;
    double coefficient = 1.0;
    for (int n = 1; n <= 100; ++n)
    { // Change the number of terms for desired precision
        result += term * coefficient;
        term *= xSquared * (2 * n - 1) / (2 * n);
        coefficient *= (2 * n - 1) / (2 * n);
    }

    return PI / 2 - result;
}

// Function to calculate the inverse trigonometric arccosine (acos)
double arccos(double x)
{
    if (x < -1.0 || x > 1.0)
        return NAN; // Not a Number for invalid input

    double result = 0.0;
    double term = x;
    double xSquared = x * x;
    double coefficient = 1.0;
    for (int n = 1; n <= 100; ++n)
    { // Change the number of terms for desired precision
        result += term * coefficient;
        term *= xSquared * (2 * n - 1) / (2 * n);
        coefficient *= (2 * n - 1) / (2 * n);
    }
    return result;
}

// Function to calculate the inverse trigonometric arctangent (atan)
double arctan(double x)
{
    if (x == 0.0)
        return 0.0;
    if (x < -1.0 || x > 1.0)
        return NAN; // Not a Number for invalid input

    double result = 0.0;
    double term = x;
    double xSquared = x * x;
    for (int n = 1; n <= 100; ++n)
    { // Change the number of terms for desired precision
        result += term / n;
        term *= -xSquared;
    }
    return result;
}

// Function to calculate the inverse trigonometric arccotangent
double arccot(double x)
{
    if (x == 0.0)
        return PI / 2;
    if (x < -1.0 || x > 1.0)
        return NAN; // Not a Number for invalid input

    double result = 0.0;
    double term = x;
    double xSquared = x * x;
    for (int n = 1; n <= 100; ++n)
    { // Change the number of terms for desired precision
        result += term / n;
        term *= -xSquared;
    }
    return PI / 2 - result;
}

// Function to calculate the inverse trigonometric arcsecant
double arcsec(double x)
{
    if (x <= -1.0 || x >= 1.0)
        return NAN; // Not a Number for invalid input

    double result = 0.0;
    double term = 1 / x;
    double xSquared = 1 / x * 1 / x;
    double coefficient = 1.0;
    for (int n = 1; n <= 100; ++n)
    { // Change the number of terms for desired precision
        result += term * coefficient;
        term *= xSquared * (2 * n - 1) / (2 * n);
        coefficient *= (2 * n - 1) / (2 * n);
    }
    return result;
}

// Function to calculate the inverse trigonometric arccosecant
double arccsc(double x)
{
    if (x <= -1.0 || x >= 1.0)
        return NAN; // Not a Number for invalid input

    double result = 0.0;
    double term = 1 / x;
    double xSquared = 1 / x * 1 / x;
    double coefficient = 1.0;
    for (int n = 1; n <= 100; ++n)
    { // Change the number of terms for desired precision
        result += term * coefficient;
        term *= xSquared * (2 * n - 1) / (2 * n);
        coefficient *= (2 * n - 1) / (2 * n);
    }
    return result;
}

void *flat(int *arr)
{
    int n = shape(arr);
    int *flatArray = (int *)malloc(n * sizeof(int));
    if (flatArray == NULL)
    {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    int *ptr = flatArray; // Pointer to the beginning of flatArray
    while (*arr != '\0')
    {
        *ptr = *arr; // Assign the value pointed to by arr to the value pointed to by ptr
        ptr++;       // Move ptr to the next element in flatArray
        arr++;       // Move arr to the next element in arr
    }

    return flatArray;
}

// Function to print a matrix
void printMatrix(int rows, int cols, int *matrix[])
{
    printf("Matrix:\n");
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function to add two matrices
void matrixAddition(int rows, int cols, int *matrix1[], int *matrix2[], int *result[])
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
}

// Function to subtract two matrices
void matrixSubtraction(int rows, int cols, int *matrix1[], int *matrix2[], int *result[])
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
}

typedef struct
{
    double real;
    double complex;
} Complex;

double *real(double *arr)
{
    double a = arr[0];
    double b = arr[1];
    double c = arr[2];
    double D = pow(b, 2) - 4 * a * c;

    if (D != 0)
    {
        double r1 = (-b - sqrt(D)) / (2 * a);
        double r2 = (-b + sqrt(D)) / (2 * a);

        double root[2] = {r1, r2};
        return root;
    }

    else
    {
        double r = -b / (2 * a);
        double root[1] = {r};
        return root;
    }
}

Complex *com(double *arr)
{
    Complex r1, r2;

    double a = arr[0];
    double b = arr[1];
    double c = arr[2];
    double D = pow(b, 2) - 4 * a * c;

    r1.real = (-b) / (2 * a);
    r2.real = (-b) / (2 * a);

    r1.complex = sqrt(-D) / (2 * a);
    r2.complex = -sqrt(-D) / (2 * a);

    Complex root[2] = {r1, r2};
    return root;
}

void root(double *arr)
{
    double a = arr[0];
    double b = arr[1];
    double c = arr[2];

    double D = pow(b, 2) - 4 * a * c;
    if (D >= 0)
    {
        real(arr);
    }
    else
    {
        com(arr);
    }
}
