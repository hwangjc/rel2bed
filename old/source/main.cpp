#include <string>
#include <iostream>
#include <stdlib.h>
#include <set>
#include "utr.h"

void compute_position(UTR &utr5, UTR &utr3, Rel_Pos &relpos, char **argv);
void calc_start_pos(std::vector<std::vector<std::pair<int,int> > > &NM_ID_container, UTR &utr5, UTR &utr3, Rel_Pos &relpos, 
										std::set<std::string> &nmIDs, std::ofstream &oututr, int &utr);
void print_out(std::ofstream &oututr, int i, int utr, Bed &outputBed, Rel_Pos &relpos);


int main(int argc, char **argv){
	//////////////////////////////////////////////////////////
	//open 5 utr starting sites, store to class // done, variable name utr5
	//open 3 utr starting sites, store to class // done, variable name utr3
	//open relative_to_start file and store information... // done, variable name relpos
	//compute position... //done
	//output to bed file //done
	//////////////////////////////////////////////////////////

	/////format for input ./main 5utr 3utr relative_pos

	//opens all necessary files
	
	if (argc == 2) {
		if (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
			std::cerr << "Usage:" << std::endl;
			std::cerr << "\t " << argv[0] <<  " 5utr_ref 3utr_ref rel_pos_file 5utr_out 3utr_out" << std::endl;
			return 0;
		}
	}

	if (argc != 6) {
		std::cerr << "Usage error:" << std::endl;
		std::cerr << "\t " << argv[0] <<  " 5utr_ref 3utr_ref rel_pos_file 5utr_out 3utr_out" << std::endl;
		return 1;
	}

	try {
		UTR utr5(argv[1]); //open 5utr file
	} catch (Notopen) {
		std::cerr << argv[1] << " cannot be opened. The file is missing or corrupted." << std::endl;
	}
	UTR utr5(argv[1]);
	//utr5.print_data();

	try {
		UTR utr3(argv[2]); //open 3utr file
	} catch (Notopen) {
		std::cerr << argv[2] << " cannot be opened. The file is missing or corrupted." << std::endl;
	}
	UTR utr3(argv[2]);
	//utr3.print_data();

	try {
		Rel_Pos test(argv[3]); //open rel pos file
	} catch (Notopen) {
		std::cerr << argv[3] << " cannot be opened. The file is missing or corrupted." << std::endl;
	}
	Rel_Pos relpos(argv[3]);
	//relpos.print_data();

	//all necessary files are open
	compute_position(utr5,utr3,relpos,argv);

	std::cout << "...done!" << std::endl;

	return 0;
}


void compute_position(UTR &utr5, UTR &utr3, Rel_Pos &relpos, char **argv){

	/*computations:
			first get NM_IDs that we need to work with // done / var name: nmIDs
			look through 5utr file to get positions for the corresponding NM_ID / get strand - or +
	*/
	std::ofstream out5utr;
	std::ofstream out3utr;
	std::string utrfilename5;
	std::string utrfilename3;

	out5utr.open(argv[4]);
	out3utr.open(argv[5]);

	std::set<std::string> nmIDs; //set of the NM_IDs that we need to do computations on
	std::pair<int,int> start_end; //start end pair for one region of one NMID
	std::vector<std::pair<int,int> > regions; //a container of all the start end regions for one NM_ID
	std::vector<std::vector<std::pair<int,int> > > NM_ID_container; // a container of NM_IDs that hold regions

	//NM_ID sorted alphabetically (ASCII)
	//regions sorted increasingly (hopefully - according to file)

	int utr = 5;

	for (int i = 0; i < relpos.size(); ++i){
		nmIDs.insert(relpos.get_NM_ID(i)); //gets all the NM_IDs we need to do computation on 
	}

	for (auto i:nmIDs){ //creates container of nmid/region/startend
		for (int j = 0; j < utr5.size(); ++j){
			if (i == utr5.get_NM_ID(j)){
				start_end = std::make_pair(utr5.get_start(j),utr5.get_end(j));
				regions.push_back(start_end);
			}
		}
		NM_ID_container.push_back(regions);
	}

	calc_start_pos(NM_ID_container,utr5,utr3,relpos,nmIDs,out5utr,utr); //5UTR calculation
	out5utr.close();

	regions.clear();
	NM_ID_container.clear();

	for (auto i:nmIDs){ //creates container of nmid/region/startend
		for (int j = 0; j < utr3.size(); ++j){
			if (i == utr3.get_NM_ID(j)){
				start_end = std::make_pair(utr3.get_start(j),utr3.get_end(j));
				regions.push_back(start_end);
			}
		}
		NM_ID_container.push_back(regions);
	}

	utr = 3;

	calc_start_pos(NM_ID_container,utr5,utr3,relpos,nmIDs,out3utr,utr); //3UTR calculation
	out3utr.close();

	return;	
}

void calc_start_pos(std::vector<std::vector<std::pair<int,int> > > &NM_ID_container, UTR &utr5, UTR &utr3, Rel_Pos &relpos,
										std::set<std::string> &nmIDs, std::ofstream &oututr, int &utr){
	
	Bed outputBed;
	for (int i = 0; i < relpos.size(); ++i){ //iterate through all the relative positions
		int position = 0;
		int it = 0;
		for (auto j:nmIDs){ //iterate through the set of NM_IDs
			if (j == relpos.get_NM_ID(i)){ //if one of the NM_IDs matches the NM_ID of the relative position, compute absolute position
				
				std::string current_chr;
				int remaining = relpos.get_start(i,utr);
				UTR utr_current;
				if (utr == 5) utr_current = utr5;
				if (utr == 3) utr_current = utr3;
				if (utr_current.get_strand(relpos.get_NM_ID(i)) == "-") { //Negative strand computation
					int k = 1;
					while (remaining > (NM_ID_container[it][NM_ID_container[it].size()-k].second - NM_ID_container[it][NM_ID_container[it].size()-k].first)){
						remaining = remaining - (NM_ID_container[it][NM_ID_container[it].size()-k].second - NM_ID_container[it][NM_ID_container[it].size()-k].first);
						++k;
					}
					position = NM_ID_container[it][NM_ID_container[it].size()-k].second - remaining;
					current_chr = utr_current.get_chr(relpos.get_NM_ID(i));
				
					outputBed.chrom = current_chr;
					outputBed.end_pos = position;
					outputBed.name = relpos.get_full_ID(i);
					outputBed.score = 0;
					outputBed.strand = "-";
					
					int length = relpos.get_seq(i,utr).length();
					while (outputBed.end_pos - length < NM_ID_container[it][NM_ID_container[it].size()-k].first){
						length = NM_ID_container[it][NM_ID_container[it].size()-k].first - (outputBed.end_pos - length);
						outputBed.start_pos = NM_ID_container[it][NM_ID_container[it].size()-k].first;
						print_out(oututr, i, utr, outputBed, relpos);
						++k;
						outputBed.end_pos = NM_ID_container[it][NM_ID_container[it].size()-k].second;
					}
	
					outputBed.start_pos = outputBed.end_pos - length;
					print_out(oututr, i, utr, outputBed, relpos);
				}
				else { //Positive strand computation
					
					int k = 0;
					while (remaining > (NM_ID_container[it][k].second - NM_ID_container[it][k].first)){
						remaining = remaining - (NM_ID_container[it][k].second - NM_ID_container[it][k].first);
						++k;
					}
					position = NM_ID_container[it][k].first + remaining;
					current_chr = utr_current.get_chr(relpos.get_NM_ID(i));

					outputBed.chrom = current_chr;
					outputBed.start_pos = position;
					outputBed.name = relpos.get_full_ID(i);
					outputBed.score = 0;
					outputBed.strand = "+";

					int length = relpos.get_seq(i,utr).length();
					while (outputBed.start_pos + length > NM_ID_container[it][k].second){
						length = (outputBed.start_pos + length) - NM_ID_container[it][k].second;
						outputBed.end_pos = NM_ID_container[it][k].second;
						print_out(oututr, i, utr, outputBed, relpos);
						++k;
						outputBed.start_pos = NM_ID_container[it][k].first;
					}

					outputBed.end_pos = outputBed.start_pos + length;
					print_out(oututr, i, utr, outputBed, relpos);
				}
			}
			++it;
		}
	}
}

void print_out(std::ofstream &oututr, int i, int utr, Bed &outputBed, Rel_Pos &relpos) {

	if (utr == 5) {
		oututr << outputBed.chrom << '\t' << outputBed.start_pos << '\t' << outputBed.end_pos << '\t' << outputBed.name << '\t'
			<< outputBed.score << '\t' << outputBed.strand << '\t' << relpos.get_start(i, utr) << '\t'
			<< relpos.get_seq(i, utr) << std::endl;
	}
	else if (utr == 3) {	
		std::string seq = relpos.get_seq(i, utr);
		std::string rev_seq;
		for (int j = (seq.length() - 1); j >= 0; --j) {
			rev_seq.push_back(seq.at(j));
		}
		oututr << outputBed.chrom << '\t' << outputBed.start_pos << '\t' << outputBed.end_pos << '\t' << outputBed.name << '\t'
			<< outputBed.score << '\t' << outputBed.strand << '\t' << rev_seq << '\t' 
			<< relpos.get_start(i, utr) << std::endl;
	}
}

//////////CURRENT STATUS///////////////////
/*
	04/14/2017
		Finished making basic UTR class; stores data from the file and is able to print it
		Success: able to open 5utr.txt and 3utr.txt
		Next: make basic class for relative_positions and make able to read in data/print data
	
	04/26/2017
		Finished making relative_positions class and is able to read in/print data
		Finished making member functions for utr and relative posiition classes
		Success: able to open relative_positions.txt, able to access data from classes
		Next: implement computations
	04/28/2017
		Finished negative computations for 5utr and 3utr
		Next: implement positive strand, implement more efficient method of computing both 5utr 3utr (functor?)
*/
