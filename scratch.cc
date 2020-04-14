#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <map>
using namespace std;


size_t rand_int(size_t min, size_t max)
{
	return min + rand() % (max + 1 - min);
}

int main()
{
	// srand (time(NULL));
	// size_t min = 2;
	// size_t max = 5;
	// map<size_t, size_t> counts;
	// for (size_t i = 0; i < 100; i++)
	// {
	// 	size_t r = rand_int(min, max);
	// 	cout << r << endl;
	// 	size_t a = 3;
	// 	size_t b = 5;
	// 	cout << a + b << endl;
	// 	counts[r] += 1;
	// }
	// map<size_t, size_t>::iterator itr;
	// for (itr = counts.begin(); itr != counts.end(); ++itr)
	// {
	// 	cout << itr->first << ": " << itr->second << endl;
	// }

	ifstream input("label");
	vector<int> label;
	int l;
	while (input >> l)
	{
		cout << l << endl;
		label.push_back(l);
	}

	return 0;

}

