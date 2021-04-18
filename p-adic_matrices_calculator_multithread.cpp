#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <thread>
#include <ctime>

using namespace std;
typedef vector < vector < int > > Matrix;
/*
* Programme to evaluate the group G_{L, q', q} of matrices modulo q. Given a lattice coset NL+p (where N is not necessarily a power of p), this programme directly evaluates
* the required q and q'.
*
* Multi-threaded algorithm, with number of threads given below. This programme may take a few hours for large inputs (say q' ~ 64) even with 6 threads.
* The algorithm used is the bruteforce algorithm, and is not very efficient.
*
* If 'Y' is entered after the first input, the programme waits for a second input.
* The .txt files it generates should be placed in the same folder as the genus_enumerator.py file.
* Input formatting is as follows:
* >> a1 a2 a3
* >> N
* >> p
* followed by a 'Y' and the next set of inputs, or any other character to stop the programme.
*/
#define NUM_THREADS 2

int mod(int a, int N) {
	if (a >= 0) return (a % N);
	return (N - ((-a) % N));
}

int det(vector<int> col1, vector<int> col2, vector<int> col3) {
	return (
		col1[0] * col2[1] * col3[2]
		- col1[0] * col2[2] * col3[1]
		+ col1[1] * col2[2] * col3[0]
		- col1[1] * col2[0] * col3[2]
		+ col1[2] * col2[0] * col3[1]
		- col1[2] * col2[1] * col3[0]
		);
}

struct Form {
	int d[3];
	Form(int d1 = 1, int d2 = 1, int d3 = 1) {
		this->d[0] = d1;
		this->d[1] = d2;
		this->d[2] = d3;
	}
	int B(vector<int>& x, vector<int>& y) {
		int s = 0;
		for (int i = 0; i <= 2; i++)
			s += this->d[i] * x[i] * y[i];
		return s;
	}
	int Q(vector<int>& x) { return this->B(x, x); }
};

int order(int n, int p) {
	int k = 0;
	while (n % p == 0) {
		k++;
		n /= p;
	}
	return k;
}

int getModulus(Form& a, int N, int p) {
	int max_ord = -1;
	for (int i = 0; i <= 2; i++) {
		int o = order(a.d[i], p);
		if (o >= max_ord) max_ord = o;
	}
	int M = 1;
	for (int index = 0; index < max_ord + order(N, p); index++) M *= p;
	//(int)pow(p, max_ord + order(N, p));
	if (p == 2) M *= 2;
	if (N % 4 == 2 && p == 2) M *= 2;
	return M;
}

struct customCompare {
	bool operator()(const int& lhs, const int& rhs) const { return lhs < rhs; }
	bool operator()(const Matrix& lhs, const Matrix& rhs) const {
		int i = 0, j = 0;
		while (i <= 2 && j <= 2 && lhs[i][j] == rhs[i][j]) {
			j++;
			if (j == 3) {
				j = 0;
				i++;
			}
		}
		if (i > 2 || j > 2 || lhs[i][j] >= rhs[i][j])
			return false;
		return true;
	}
};

void getColumns(Form& a, int N, int p, int M, vector<Matrix>& cols) {
	cout << "M = " << M << endl;
	for (int i = 0; i <= 2; i++) {
		bool optimized = false;
		for (int j = 0; j < i; j++) {
			if (mod(a.d[i] - a.d[j], M) == 0) {
				optimized = true;
				cols[i] = cols[j];
			}
		}
		if (!optimized) {
			vector<int> x(3, 0);
			for (x[0] = 0; x[0] < M; x[0]++) {
				for (x[1] = 0; x[1] < M; x[1]++) {
					for (x[2] = 0; x[2] < M; x[2]++) {
						if (mod(a.Q(x) - a.d[i], M) == 0)
							cols[i].push_back(x);
					}
				}
			}
		}
	}
}

struct Parameters {
	Form a;
	int N;
	int p;
	vector<vector<vector<int>>> Cols;
	int num_thread;
	int i;
	int j;
	int k;
	int M;
	set<Matrix, customCompare> output;
};

int sign(int i, int j, int k) {
	vector<int> col0 = { 0,0,0 }, col1 = { 0,0,0 }, col2 = { 0,0,0 };
	col0[i] = 1;
	col1[j] = 1;
	col2[k] = 1;
	return det(col0, col1, col2);
}

void findMatrices(Parameters& param) {
	Form a = param.a;
	int N = param.N;
	int p = param.p;
	vector<vector<vector<int>>> cols = param.Cols;
	int num_thread = param.num_thread;
	cout << "Created Thread " << num_thread << endl;
	int i = param.i;
	int j = param.j;
	int k = param.k;
	int M = param.M;

	long long int iter = 0;
	long long int tot_iter = (long long int)cols[0].size() * (long long int)cols[1].size() * (long long int)cols[2].size();

	int percentage_complete = 0;
	//cout << "i, j, k = " << i << j << k << endl;
	Matrix tempmatrix(3, vector<int>(3));
	//set <Matrix, customCompare> matrices
	long long int col1Size = cols[i].size();
	//for (vector<vector<int>>::iterator col1 = cols[i].begin() + num_thread; col1 != cols[i].end(); col1 += NUM_THREADS) {
	int s = sign(i, j, k);
	//cout << s;
	for (int index = num_thread; index < col1Size; index += NUM_THREADS) {
		vector<int> col1 = cols[i][index];
		for (vector<vector<int>>::iterator col2 = cols[j].begin(); col2 != cols[j].end(); col2++) {
			if (mod(a.B(col1, *col2), M) == 0) {
				for (vector<vector<int>>::iterator col3 = cols[k].begin(); col3 != cols[k].end(); col3++) {
					iter++;
					if (mod(a.B(*col2, *col3), M) == 0 && mod(a.B(col1, *col3), M) == 0 && mod(s * det(col1, *col2, *col3), M) == 1) {
						tempmatrix[i] = col1;
						tempmatrix[j] = *col2;
						tempmatrix[k] = *col3;
						for (int r = 0; r <= 2; r++) {
							for (int c = 0; c <= 2; c++)
								tempmatrix[r][c] = mod(tempmatrix[r][c], N);
						}
						param.output.insert(tempmatrix);
					}
					if (NUM_THREADS * 100 * iter >= tot_iter) {
						percentage_complete++;
						cout << "Thread " << num_thread << " - " << percentage_complete << "% done||" << flush;
						iter = 0;
					}
				}
			}
			else {
				iter += (long long int)cols[k].size();
			}

			if (NUM_THREADS * 100 * iter >= tot_iter) {
				percentage_complete++;
				cout << "Thread " << num_thread << " - " << percentage_complete << "% done||" << flush;
				iter = 0;
			}
		}
	}
	//param.output = matrices;
}



int main() {
	char cont;
	cin >> cont;

	while (cont == 'Y') {
		int d1, d2, d3, N, p;
		cin >> d1 >> d2 >> d3;
		cin >> N; cin >> p;
		int prod = 1;
		while (N % p == 0) {
			prod *= p;
			N /= p;
		}
		N = prod;
		Form a = { d1, d2, d3 };
		int M = getModulus(a, N, p);
		vector<Matrix> cols(3, vector<vector<int>>(0, vector<int>(3, 0)));
		//getColumns(Form& a, int N, int p, vector<Matrix>& cols)
		//auto start = chrono::high_resolution_clock::now();
		//auto start_t = chrono::high_resolution_clock::to_time_t(start);
		//cout << "Started at " << ctime(&start_t);
		getColumns(a, N, p, M, cols);
		long long int numcol1 = cols[0].size();
		long long int numcol2 = cols[1].size();
		long long int numcol3 = cols[2].size();
		int i, j, k;
		if (numcol1 >= numcol2 && numcol1 >= numcol3) k = 0;
		else if (numcol2 >= numcol1 && numcol2 >= numcol3) k = 1;
		else k = 2;
		i = mod(k - 1, 3);
		j = mod(k + 1, 3);
		cout << i << ' ' << j << ' ' << k << '\n';

		thread allThreads[NUM_THREADS];
		Parameters params[NUM_THREADS];
		for (int t = 0; t < NUM_THREADS; t++) {
			params[t].a = a;
			params[t].i = i;
			params[t].j = j;
			params[t].k = k;
			params[t].N = N;
			params[t].p = p;
			params[t].M = M;
			params[t].num_thread = t;
			params[t].Cols = cols;
			if (t != NUM_THREADS - 1)
				allThreads[t] = thread(findMatrices, ref(params[t]));
			else
				findMatrices(params[NUM_THREADS - 1]);
		}
		for (int t = 0; t < NUM_THREADS - 1; t++)
			allThreads[t].join();

		cout << endl;
		ofstream fout;
		ostringstream filename;
		filename << p << "-adic matrices a=(" << a.d[0] << ',' << a.d[1] << ',' << a.d[2] << ") N=" << N << "(1).txt";
		fout.open(filename.str());

		set<Matrix, customCompare> matrices;
		for (int t = 0; t < NUM_THREADS; t++)
			matrices.insert(params[t].output.begin(), params[t].output.end());
		int num_matrices = 0;

		while (!fout);
		for (set<Matrix, customCompare>::iterator it = matrices.begin(); it != matrices.end(); it++) {
			num_matrices++;
			for (int r = 0; r <= 2; r++) {
				for (int c = 0; c <= 2; c++) {
					fout << (*it)[c][r] << ' ';
				}
				fout << endl;
			}
			fout << "\n";//mod(det((*it)[0], (*it)[1], (*it)[2]), N) << "\n\n";
		}
		fout.close();

		//auto stop = chrono::high_resolution_clock::now();
		//auto stop_t = chrono::high_resolution_clock::to_time_t(stop);

		//auto hours = chrono::duration_cast<chrono::hours>(stop - start);
		//auto minutes = chrono::duration_cast<chrono::minutes>(stop - start);
		//auto seconds = chrono::duration_cast<chrono::seconds>(stop - start);

		cout << "Number of Matrices = " << num_matrices << endl;
		//cout << "Ended at " << ctime(&stop_t) << "Duration: " << hours.count() << "h " << minutes.count() << "m " << seconds.count() << "s\n";
		//cout << "\n\n\n";

		cin >> cont;
	}
	return 0;
}

