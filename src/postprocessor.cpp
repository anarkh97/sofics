#include <iostream>
#include <iomanip>
#include <iterator>
#include <fstream>
#include <sstream>
#include <climits>
#include <cassert>
#include <string>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char* argv[])
{

	if(argc != 5) {
		cout << "Usage: [binary] <topo-file> <mass> <epstrain>"
                     << " <results>\n";
		exit(-1);
	}

	ifstream input1(argv[1], ios::in);
	ifstream input2(argv[2], ios::in);
	ifstream input3(argv[3], ios::in);
	ofstream output(argv[4], ios::out);

	// common string for reading
	string line;

	set<int> node_idx;
	input1.ignore(512, '\n'); // ignore header
	int r;
	for(r=0; r<INT_MAX; ++r) {
		getline(input1, line);
		istringstream iss(line);

		int id, type;
		iss >> id >> type;

		int num_nodes = (type == 3) ? 3 : 4;
		vector<int> data(num_nodes, -1);

		bool done = false;
		for(int i=0; i<num_nodes; ++i) {
			iss >> data[i];
			if(iss.fail()) {
				done = true;
				break;
			}
		}

		if(done) break;
		for(int i=0; i<num_nodes; ++i)
			node_idx.insert(data[i]-1);	
	}
        
	getline(input2, line);

	double mass_value = stod(line);
	output << "    ";
	output << scientific << setprecision(6) 
	       << mass_value << "  " 
	       << "mass" << endl;

	// reading plastic strain
	int num_rows = 0;
	int num_data_rows = node_idx.size();
	input3.ignore(512, '\n'); // ignore header
	getline(input3, line); // number of elements/nodes
	istringstream is(line);
	is >> num_rows;

	assert(num_data_rows <= num_rows);
	
	bool done = false;
	double max_plastic_strain = 0;
	//int max_id = 0;
	while(!done) {
		is.clear(); is.str("");
		// get time stamp
		double time;
		getline(input3, line);
		is.str(line);	
		is >> time;

		if(is.fail()) {
			done = true;
			break;
		}
		// start reading values
		vector<double> data(num_data_rows, -1);
		int idx=0;
		for(int i=0; i<num_rows; ++i) {
			getline(input3, line);
			if(node_idx.find(i) != node_idx.end()) {
				is.clear();
				is.str(line);
				is >> data[idx];
				++idx;
				if(is.fail()) {
					done = true;
					break;
				}
			}
		}

		// find max
		auto max_ret = max_element(data.begin(), 
		  data.end());
	
		double possible_max = *max_ret;

		// compare with previous max
		if(possible_max >= max_plastic_strain) {
			max_plastic_strain = possible_max;

/*
			// Debug: print the node Id for this max
			auto offset = distance(data.begin(), max_ret);
			auto node_iter = node_idx.begin();
			advance(node_iter, offset);
			max_id = *node_iter;
*/
		}
	}

	output << "    ";
	output << scientific << setprecision(6) 
	       << max_plastic_strain << "  " 
	       << "max_plastic_strain" << endl;

	//cout << "Max_ID = " << max_id << endl;

	input1.close();
	input2.close();
	input3.close();
	output.close();

	return 0;

}
