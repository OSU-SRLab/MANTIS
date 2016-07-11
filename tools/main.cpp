#include "circularLinkedList.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::uppercase;
using std::cout;
using std::cerr;

using namespace mantis;

// #define IGNORE_COMMENTS

inline bool isBase(char c)
{
	if (c == 'A'|| c == 'T' || c == 'C' || c == 'G' || c == 'N') return true;
	else return false;
}

inline void outputLocus(CircularLinkedList* cll, ofstream& ofs, string chr,
	unsigned long curPos)
{
	unsigned int seqLength = cll->getLastLength();
	int repeatCount = cll->getRepeatCount();
	string kmer = cll->getKmer();
	unsigned long start = curPos - seqLength;
	unsigned long end = curPos; //BED is half-open!
	ofs << chr << "\t" << start << "\t" << end << "\t" << "(" << kmer << ")"
		<< repeatCount << "\t0\t+" << endl;
}

int main (int argc, char** argv)
{
	int minRepeatLength = 1;
	int maxRepeatLength = 5;
	int minBases = 10;
	int minRepeats = 3;
	int opt;
	string fasta;
	string output;
	bool fastaProvided = false;
	bool outputProvided = false;
	extern char* optarg;
	if (argc <= 1 || strncmp(argv[1],"-h",2) == 0)
	{
		cout << "RepeatFinder 1.0" << endl;
		cout << "\t-m: minimum number of bases a repeat region must span to"
			<< " be called a microsatellite. Default: 10" << endl;
		cout << "\t-r: minimum number of repeats for a microsatellite to be"
			<< " called. Default: 3" << endl;
		cout << "\t-l: minimum k-mer length. Default: 1" << endl;
		cout << "\t-L: maximum k-mer length. Default: 5" << endl;
		cout << endl;
		cout << "\t-i: input FASTA file" << endl;
		cout << "\t-o: output microsatellites file" << endl;
		return 1;
	}
	while((opt = getopt(argc, argv,"m:r:l:L:i:o:")) != -1)
	{
		switch(opt)
		{
			case 'm':
				minBases = atoi(optarg);
				break;
			case 'r':
				minRepeats = atoi(optarg);
				break;
			case 'l':
				minRepeatLength = atoi(optarg);
				break;
			case 'L':
				maxRepeatLength = atoi(optarg);
				break;
			case 'i':
				fasta = optarg;
				fastaProvided = true;
				break;
			case 'o':
				output = optarg;
				outputProvided = true;
				break;
			default:
				//cerr << "ERROR: invalid option " << (char) opt << endl;
				return 1;
		}
	}
	if (!fastaProvided) cerr << "ERROR: Must supply input FASTA." << endl;
	if (!outputProvided) cerr << "ERROR: Must supply output path." << endl;
	if (!fastaProvided || !outputProvided) return 1;

	ifstream fastaFile(fasta.c_str());
	ofstream outputFile(output.c_str());
	int numRepeatFinders = maxRepeatLength - minRepeatLength + 1;
	CircularLinkedList** clls = new CircularLinkedList*[numRepeatFinders];
	for (int i = 0; i < numRepeatFinders; i++)
	{
		clls[i] = new CircularLinkedList(i + minRepeatLength, minBases,
			minRepeats, "N");
	}

	string curChr;
	unsigned long curPos = 0;
	char x;
	bool first = true;
	while (fastaFile.get(x))
	{
		if (x != '>')
		{
			//safer
			/* http://www.aosabook.org/en/posa/
				working-with-big-data-in-bioinformatics.html */
			//x &= 0xdf;
			//if (isBase(x))
			
			//potentially faster, but expects a perfectly formed fasta file
			//with LF terminators
			if (x != '\n')
			{
#ifdef IGNORE_COMMENTS
				if (x != ';')
				{
#endif
					x &= 0xdf;
					bool locusFound = false;
					for (int i = 0; i < numRepeatFinders; i++)
					{
						if (!locusFound)
						{
							bool result = clls[i]->next(x);
							if (result)
							{
								outputLocus(clls[i], outputFile,
									curChr, curPos);
								locusFound = true;
							}
						}
						else clls[i]->nextNoResultCheck(x);
					}
					curPos++;
#ifdef IGNORE_COMMENTS
				}
				else
				{
					while (x != '\n')
					{
						fastaFile.get(x);
					}
				}
#endif
			}
		}
		else
		{
			bool gettingChrName = true;
			if (!first)
			{
				for (int i = 0; i < numRepeatFinders; i++)
				{
					bool result = clls[i]->endSequence();
					if (result)
					{
						outputLocus(clls[i], outputFile, curChr, curPos);
						break;
					}
				}
				for (int i = 0; i < numRepeatFinders; i++)
				{
					delete clls[i];
					clls[i] = new CircularLinkedList(i + minRepeatLength,
						minBases, minRepeats, "N");
				}
			}
			else first = false;
			curChr = "";
			curPos = 0;
			while (x != '\n')
			{
				fastaFile.get(x);
				if (gettingChrName)
				{
					if (/* x == '\r' || */ x == '\n' || x == ' ')
					{
						gettingChrName = false;
					}
					else curChr += x;
				}
			}
		}
	}
	for (int i = 0; i < numRepeatFinders; i++)
	{
		bool result = clls[i]->endSequence();
		if (result)
		{
			outputLocus(clls[i], outputFile, curChr, curPos);
			break;
		}
	}
	outputFile.close();
	fastaFile.close();
	for (int i = 0; i < numRepeatFinders; i++)
	{
		delete clls[i];
	}
	delete[] clls;
	return 0;
}
