#include <iostream>
#include <math.h>
#include <thread>
#include <vector>
#include <fstream>
#include <ctime>

/*
* Code to check an identity of the form r_{a, h, N}(n) = (r_ratio_num / r_ratio_den) r_a(n) + character[sqrt(n/t)] * (t_ratio_num / t_ratio_den) * sqrt(n/t) for all n congruent to m modulo M.
* Here, sqrt(n/t) defaults to 0 if t does not divide n, or if n/t is not a perfect square. 'character' is a Dirichlet character.
* This is implemented as a multi-threaded programme, with the number of threads '#define'd below.
* Implemented in C++ as opposed to python for efficiency reasons
* For any unexplained notation see accompanying paper.
* Input text file requires the following format (all integers, space separated. See example input given below):
* a1 a2 a3
* N
* h1 h2 h3
* m M
* <number of coefficients needed to be checked>
* r_ratio_num  r_ratio_den
* t (Note: if t = 0, then program stops scanning for inputs, and ignores the unary theta series term)
* t_ratio_num  t_ratio_den
* <modulus of character>
* <list of outputs of the character as space separated integers, for the inputs 0, 1, 2, ..., <modulus>-1 >
*
* If more inputs needed to be given in the same run, input 'Y' followed by the next set of inputs. 'N' (or any other character other than 'Y') stops the programme.
*
* E.g: Want to check first for a= (1,1,1), N=4, h = (1,0,0), and want to check for n = 1 mod 8, and then n = 2 mod 8, upto n <= 5000. The expected identity is given in the paper in the Introduction Section.
* Ideal Input format:
* >> 1 1 1
* >> 4
* >> 1 0 0
* >> 1 8
* >> 5000
* >> 1 12
* >> 1
* >> 1 2
* >> 4
* >> 0 1 0 -1
* >>
* >> Y
* >>
* >> 1 1 1
* >> 4
* >> 1 0 0
* >> 2 8
* >> 5000
* >> 0 1
* >> 0
* >> N
*/

using namespace std;
#define NUM_THREADS 3

bool isSquare(int n) {
	float sr = sqrt(n);
	return ((sr - floor(sr)) == 0);
}

bool isBadSquareClass(int n, int t) {
	if (t == 0) { return false; }
	if (n % t != 0) { return false; }
	if (isSquare(n / t)) { return true; }
}

int mod(int a, int N) {
	if (a >= 0) return (a % N);
	int neg_res = (-a) % N;
	if (neg_res == 0) return 0;
	else return (N - neg_res);
}

int countincr(int abs) { return (abs == 0 ? 1 : 2); }

auto findNumSolutions(int n, int a0, int a1, int a2, int N, int h0, int h1, int h2) {
	int count_nocond = 0, count_withcond = 0;
	//int lowerboundx = floor(-sqrt((float)n / a0));
	int upperboundx = ceil(sqrt((float)n / a0));
	for (int abs_x = 0; abs_x <= upperboundx; abs_x++) {
		//int lowerboundy = floor(-sqrt((float)(n - a0 * x * x) / a1));
		int upperboundy = ceil(sqrt((float)(n - a0 * abs_x * abs_x) / a1));
		for (int abs_y = 0; abs_y <= upperboundy; abs_y++) {
			int a2zsq = (n - a0 * abs_x * abs_x - a1 * abs_y * abs_y);
			if (a2zsq % a2 == 0 && isSquare(a2zsq / a2)) {
				int abs_z = (int)sqrt(a2zsq / a2);
				count_nocond += (countincr(abs_x) * countincr(abs_y) * countincr(abs_z));
				if (mod(abs_x, N) == h0 && mod(abs_y, N) == h1) {
					if (abs_z % N == h2) count_withcond++;
					if (abs_z != 0 && mod(-abs_z, N) == h2) count_withcond++;
				}
				if ((abs_x != 0 && mod(-abs_x, N) == h0) && mod(abs_y, N) == h1) {
					if (abs_z % N == h2) count_withcond++;
					if (abs_z != 0 && mod(-abs_z, N) == h2) count_withcond++;
				}
				if (mod(abs_x, N) == h0 && (abs_y != 0 && mod(-abs_y, N) == h1)) {
					if (abs_z % N == h2) count_withcond++;
					if (abs_z != 0 && mod(-abs_z, N) == h2) count_withcond++;
				}
				if ((abs_x != 0 && mod(-abs_x, N) == h0) && (abs_y != 0 && mod(-abs_y, N) == h1)) {
					if (abs_z % N == h2) count_withcond++;
					if (abs_z != 0 && mod(-abs_z, N) == h2) count_withcond++;
				}
			}
		}
	}
	return make_pair(count_nocond, count_withcond);
}

struct Parameters {
	int a1;
	int a2;
	int a3;
	int N;
	int h1;
	int h2;
	int h3;
	int m;
	int M;
	int targetCoeffs;
	int r_ratio_num;
	int r_ratio_den;
	int t = 0;
	int t_ratio_num = 0;
	int t_ratio_den = 0; 
	int char_mod = 1;
	int num_thread;
	int flag = -1;
	vector<int> character;
};

void verifyIdentity(Parameters &Param) {
	cout << "Creating Thread " << Param.num_thread;
	int a1 = Param.a1;
	int a2 = Param.a2;
	int a3 = Param.a3;
	int N = Param.N;
	int h1 = Param.h1;
	int h2 = Param.h2;
	int h3 = Param.h3;
	int m = Param.m; 
	int M = Param.M;
	int targetCoeffs = Param.targetCoeffs; 
	int r_ratio_num = Param.r_ratio_num; 
	int r_ratio_den = Param.r_ratio_den;
	int t = Param.t;
	int t_ratio_num = Param.t_ratio_num; 
	int t_ratio_den = Param.t_ratio_den;
	int char_mod = Param.char_mod;
	vector<int> character = Param.character;

	int tot_iter = (targetCoeffs - m - (M * Param.num_thread) + 1);
	int percent_complete = 0, iter = 0;
	for (int n = m + (M * Param.num_thread); n <= targetCoeffs; n += M * NUM_THREADS) {
		pair<int, int> temp = findNumSolutions(n, a1, a2, a3, N, h1, h2, h3);
		int count_nocond = temp.first, count_withcond = temp.second;
		if (!isBadSquareClass(n, t)) {
			if (count_nocond * r_ratio_num != count_withcond * r_ratio_den) {
				Param.flag = n;
				cout << "| Thread " << Param.num_thread << " failed at n = " << n << " |";
				return;
			}
		}
		else {
			int k = (int)sqrt(n / t);
			int lhs = (count_withcond * r_ratio_den - count_nocond * r_ratio_num) * t_ratio_den;
			int rhs = r_ratio_den * character[mod(k, char_mod)] * t_ratio_num * k;
			if (lhs != rhs) {
				Param.flag = n;
				cout << "| Thread " << Param.num_thread << " failed at n = " << n << " |";
				return;
			}
		}
		iter++;
		if (100 * (M * NUM_THREADS) * iter >= tot_iter) {
			percent_complete++;
			cout << "| Thread " << Param.num_thread << ' ' << percent_complete << "% done |";
			iter = 0;
		}
	}
	Param.flag = -1;
}

int main() {
	ofstream output;
	for (char nextInput = 'Y'; nextInput == 'Y'; cin >> nextInput) {
		int a1, a2, a3; cin >> a1 >> a2 >> a3;
		int N; cin >> N;
		int h1, h2, h3; cin >> h1 >> h2 >> h3;
		int m; cin >> m;
		int M; cin >> M;
		int targetCoeffs; cin >> targetCoeffs;
		int r_ratio_num; cin >> r_ratio_num;
		int r_ratio_den; cin >> r_ratio_den;
		int t; cin >> t;
		int t_ratio_num = 0; 
		int t_ratio_den = 1;
		int char_mod = 1;
		vector<int> character = vector<int>(1,0);
		if (t != 0) {
			cin >> t_ratio_num >> t_ratio_den;
			cin >> char_mod;
			character = vector<int>(char_mod, 0);
			for (int i = 0; i < char_mod; i++)
				cin >> character[i];
		}
		
		//auto start = chrono::high_resolution_clock::now();
		//auto start_t = chrono::high_resolution_clock::to_time_t(start);
		//cout << "Started at " << ctime(&start_t);
		thread allThreads[NUM_THREADS];
		Parameters params[NUM_THREADS];
		for (int thr = 0; thr < NUM_THREADS; thr++) {
			params[thr].a1 = a1;
			params[thr].a2 = a2;
			params[thr].a3 = a3;
			params[thr].N = N;
			params[thr].h1 = h1;
			params[thr].h2 = h2;
			params[thr].h3 = h3;
			params[thr].m = m;
			params[thr].M = M;
			params[thr].targetCoeffs = targetCoeffs;
			params[thr].num_thread = thr;
			params[thr].t = t;
			params[thr].r_ratio_num = r_ratio_num;
			params[thr].r_ratio_den = r_ratio_den;
			params[thr].t_ratio_num = t_ratio_num;
			params[thr].t_ratio_den = t_ratio_den;
			params[thr].char_mod = char_mod;
			params[thr].flag = -1;
			params[thr].character = character;
			if (thr != NUM_THREADS - 1)
				allThreads[thr] = thread(verifyIdentity, ref(params[thr]));
			else
				verifyIdentity(params[NUM_THREADS - 1]);
		}
		for (int thr = 0; thr < NUM_THREADS - 1; thr++)
			allThreads[thr].join();
		cout << '\n';
		output.open("output.txt", ios::out | ios::app);
		output << "a = " << a1 << ' ' << a2 << "; N = " << N << "; h = " << h1 << ' ' << h2 << ' ' << h3 << "; m = " << m << "; M = " << M << '\n';
		for (int thr = 0; thr < NUM_THREADS; thr++) {
			output << "Thread " << thr << " returned: ";
			if (params[thr].flag != -1) output << "Identity fails at n = " << params[thr].flag;
			else output << "Identity Holds";
			output << '\n';
		}
		output << '\n';
		output.close();

		/*auto stop = chrono::high_resolution_clock::now();
		auto stop_t = chrono::high_resolution_clock::to_time_t(stop);

		auto hours = chrono::duration_cast<chrono::hours>(stop - start);
		auto minutes = chrono::duration_cast<chrono::minutes>(stop - start);
		auto seconds = chrono::duration_cast<chrono::seconds>(stop - start);

		cout << "Ended at " << ctime(&stop_t) << "Duration: " << hours.count() << "h " << minutes.count() << "m " << seconds.count() << "s\n";
		cout << "\n\n\n";*/

	}
	return 0;
}