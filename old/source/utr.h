#include <fstream>
#include <string>
#include <iostream>
#include <vector>

struct Bed {
	std::string chrom;
	int start_pos;
	int end_pos;
	std::string name;
	int score;
	std::string strand;
};

class Notopen{}; // Exception throw

class UTR {
public:
	
	UTR(){}

	UTR(const std::string filename) { //constructor open with filename
		//openfile here
		std::ifstream inFile;
		inFile.open(filename.c_str());
		if (!inFile.is_open()) {
			throw Notopen();
		}

		std::string chr,id,strd;
		int start,end,score;

		while(inFile >> chr >> start >> end >> id >> score >> strd) {
			row inrow;
			inrow.chrom = chr;
			inrow.start_pos = start;
			inrow.end_pos = end;
			inrow.name = id;
			inrow.score = score;
			inrow.strand = strd;

			rows.push_back(inrow);
		}
		inFile.close();
	}

	//EFFECTS: prints the utr as it was read
	void print_data() const{
		for (unsigned int i = 0; i < rows.size(); ++i){
			std::cout << rows[i].chrom << " " << rows[i].start_pos << " " << rows[i].end_pos 
			<< " " << rows[i].name << " " << rows[i].score << " "  << rows[i].strand << std::endl;
		}
		return;
	}
	
	//EFECTS: returns how many rows there are
	int size() const{
		return rows.size();
	}

	std::string get_chr(const int &index) const{ //returns chromosome from index
		return rows[index].chrom;
	}

	std::string get_chr(const std::string &nmID) const{ //returns chromosome from NM ID
		for (int i = 0; i < size(); ++i){
			if (get_NM_ID(i) == nmID) return rows[i].chrom;
		}
		return "chrXX";
	}

	int get_start(const int &index) const{
		return rows[index].start_pos;
	}

	int get_end(const int &index) const{
		return rows[index].end_pos;
	}

	std::string get_strand(const int &index) const{
		return rows[index].strand;
	}

	std::string get_strand(const std::string &nmID) const{ //returns strand from NM ID
		for (int i = 0; i < size(); ++i){
			if(get_NM_ID(i) == nmID) return rows[i].strand;
		}
		return "+";
	}

	std::string get_NM_ID(const int &index) const{
		std::string nmID = rows[index].name;
		std::string returnval = nmID.substr(0,3); //NM_
		int i = 3;
		while (nmID[i] != '_'){
			returnval.push_back(nmID[i]);
			++i;
		}
		return returnval;
	}


private:
	struct row {
		std::string chrom;
		int start_pos;
		int end_pos;
		std::string name;
		int score;
		std::string strand;
	};

	std::vector<row> rows; 
};


////////////////////////////////////////////RELATIVE POSITION CLASS////////////////////////////////////


class Rel_Pos {
public:

	Rel_Pos(const std::string filename){
		std::ifstream inFile;
		inFile.open(filename.c_str());
		if (!inFile.is_open()){
			throw Notopen();
		}

		std::string nmID_in, seq_5_in, seq_3_in;
		int rel_5_in, rel_3_in;
		double val_in;

		while(inFile >> nmID_in >> rel_5_in >> rel_3_in >> seq_5_in >> seq_3_in >> val_in){
			row inrow;
			inrow.nmID = nmID_in;
			inrow.rel_start_5 = rel_5_in;
			inrow.rel_start_3 = rel_3_in;
			inrow.seq_5 = seq_5_in;
			inrow.seq_3 = seq_3_in;
			inrow.val = val_in;

			rows.push_back(inrow);
		}
		inFile.close();
		return;
	}

	//EFFECTS: prints to std::cout the relative positions file
	void print_data() const{
		for(unsigned int i = 0; i < rows.size(); ++i){
			std::cout << rows[i].nmID << " " << rows[i].rel_start_5 << " " << rows[i].rel_start_3 << " " << rows[i].seq_5
				<< " " << rows[i].seq_3 << " " << rows[i].val << std::endl;
		}
	}

	//EFFECTS: returns number of rows in data
	int size() const{
		return rows.size();
	}

	std::string get_NM_ID(const int & index) const{
		std::string returnval = rows[index].nmID.substr(0,3); //NM_
		int i = 3;
		while (rows[index].nmID[i] != '_' && rows[index].nmID[i] != '-'){
			returnval.push_back(rows[index].nmID[i]);
			++i;
		}
		return returnval;
	}

	std::string get_full_ID(const int & index) const{
		return rows[index].nmID;
	}

	int get_5_start(const int & index) const{
		return rows[index].rel_start_5;
	}

	int get_3_start(const int & index) const{
		return rows[index].rel_start_3;
	}
	
	//REQUIRES: utr == 3 || utr == 5
	//EFFECTS: returns the start position based on index and utr 
	int get_start(const int & index, const int & utr){
		if (utr == 3) return get_3_start(index);
		else if (utr == 5) return get_5_start(index);
		return 0;
	}

	std::string get_5_seq(const int & index) const{
		return rows[index].seq_5;
	}
	
	std::string get_3_seq(const int & index) const{
		return rows[index].seq_3;
	}
	//REQUIRES: utr == 3 || utr == 5
	//EFFECTS: returns the sequence based on index and utr
	std::string get_seq(const int & index, const int & utr){
		if (utr == 3) return get_3_seq(index);
		else if (utr == 5) return get_5_seq(index);
		return "";
	}

	double get_val(const int & index) const{
		return rows[index].val;
	}

private:
	struct row {
		std::string nmID;
		int rel_start_5;
		int rel_start_3;
		std::string seq_5;
		std::string seq_3;
		double val;
	};

	std::vector<row> rows;
};


