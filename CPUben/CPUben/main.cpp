#include <iostream>
#include <chrono>
#include <cstdlib>
using namespace std;
using namespace std::chrono;

//About 1.03 seconds
#define N 9
//About 1.02 seconds
#define M 23000

double matrix1[N][N] = {
	{2,1,1,1,1,1,1,1,1},
	{1,2,1,1,1,1,1,1,1},
	{1,1,2,1,1,1,1,1,1},
	{1,1,1,2,1,1,1,1,1},
	{1,1,1,1,2,1,1,1,1},
	{1,1,1,1,1,2,1,1,1},
	{1,1,1,1,1,1,2,1,1},
	{1,1,1,1,1,1,1,2,1},
	{1,1,1,1,1,1,1,1,2}
};

double matrix2[N][N];

int array1[M];

//Functions required to find the inverse of a matrix
void getCofactor(double A[N][N], double temp[N][N], int p, int q, int n);
int determinant(double A[N][N], int n);
void adjoint(double A[N][N], double adj[N][N]);
bool DOUBLE(double A[N][N], double inverse[N][N]);
template<class T>
void display(T A[N][N]);

//Functions for quicksort
void fill();
void swap(int* a, int *b);
int partition(int arr[], int low, int high);
void INTEGER(int arr[], int low, int high);

//Function to find current second
double sec(void);

int main()
{
	int NINT = 0, NFLOAT = 0;
	double START, VINT, VFLOAT, AverageSpeed;

	START = sec();
	
	
	while (sec() < START + 10) //Runs for 10 Seconds
	{
		//1.03 Seconds
		DOUBLE(matrix1, matrix2);
		NFLOAT++;
	}
	VFLOAT = 60 * NFLOAT / (sec() - START);
	

	
	START = sec();
	while (sec() < START + 10) //Runs for 10 Seconds
	{
		//Approximately 1.02 seconds
		fill();
		INTEGER(array1, 0, M);
		NINT++;
	}
	VINT = 60 * NINT / (sec() - START);
	AverageSpeed = (2 * VFLOAT * VINT) / (VFLOAT + VINT);
	
	cout << "VFloat: " << VFLOAT << endl;
	cout << "VINT: " << VINT << endl;
	cout << "AverageSpeed: " << AverageSpeed << endl;
	system("pause");
	return 0;
}

void getCofactor(double A[N][N], double temp[N][N], int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != p && col != q)
			{
				temp[i][j++] = A[row][col];

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
int determinant(double A[N][N], int n)
{
	int D = 0; // Initialize result 

	//  Base case : if matrix contains single element 
	if (n == 1)
		return A[0][0];

	double temp[N][N]; // To store cofactors 

	int sign = 1;  // To store sign multiplier 

	 // Iterate for each element of first row 
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of A[0][f] 
		getCofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign 
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(double A[N][N], double adj[N][N])
{
	if (N == 1)
	{
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][] 
	double sign = 1, temp[N][N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// Get cofactor of A[i][j] 
			getCofactor(A, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row 
			// and column indexes is even. 
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the 
			// transpose of the cofactor matrix 
			adj[j][i] = (sign)*(determinant(temp, N - 1));
		}
	}
}

// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool DOUBLE(double A[N][N], double inverse[N][N])
{
	// Find determinant of A[][] 
	int det = determinant(A, N);
	if (det == 0)
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint 
	double adj[N][N];
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			inverse[i][j] = adj[i][j] / double(det);

	return true;
}

//Need this function to display elements of array of doubles
template<class T>
void display(T A[N][N])
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}
}

void fill()
{
	srand(time(NULL));
	for (int i = 0; i < M; i++)
	{
		array1[i] = (rand() % 10);
	}
}

void swap(int* a, int *b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

int partition(int arr[], int low, int high)
{
	int pivot = arr[high];    // pivot 
	int i = (low - 1);  // Index of smaller element 

	for (int j = low; j <= high - 1; j++)
	{
		// If current element is smaller than or 
		// equal to pivot 
		if (arr[j] <= pivot)
		{
			i++;    // increment index of smaller element 
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}

void INTEGER(int arr[], int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
		   at right place */
		int pi = partition(arr, low, high);

		// Separately sort elements before 
		// partition and after partition 
		INTEGER(arr, low, pi - 1);
		INTEGER(arr, pi + 1, high);
	}
}

double sec(void)
{
	return double(clock()) / double(CLOCKS_PER_SEC);
}