#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include "marker.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Global variables
vector<char *> popNames;

////////////////////////////////////////////////////////////////////////////////
// Used to store necessary details about the source and length of each segment:
struct Segment { int popNum, endMarker; }; // note: start marker is implict
typedef vector<Segment>*   SegmentList;

////////////////////////////////////////////////////////////////////////////////
// Function decls
unsigned int getRandSeed();
void readDatAndSimulate(char *datFile, char *outFile);
void simulate(int numSampsToSimulate, double *populationProportion, int numPops,
	      vector<SegmentList> &simuOutput,
	      vector<SegmentList> &prevSimulated, int prevGeneration,
	      bool useOnlyAdmixed = false);
bool decideIfRecomb(double geneticDistance);
void recordSegment(SegmentList simuOutput, int popNum, int copyChrom,
		   int endMarker, vector<SegmentList> &prevSimulated,
		   int numSampsToSimulate);
int  choosePop(double *populationProportion, int numPops);
void swapSimuListsAndClear(vector<SegmentList> *&prevSimulated,
			   vector<SegmentList> *&simuOutput,
			   int numSampsToSimulate);
void output(char *outFile, vector<SegmentList> &simulated,
	    int numSamples);
void printUsage(char **argv);



////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  if (argc != 4) {
    printUsage(argv);
  }

  // seed random number generator:
  unsigned int seed;
  seed = getRandSeed();
//  seed = 4040981314u; // for testing
  printf("Using random seed: %u\n\n", seed);
  srand(seed);

  // read in the map:
  printf("Reading marker file... ");
  fflush(stdout);
  Marker::readMarkerFile(/*markerFile=*/ argv[2]);
  printf("done.\n");

  // read the spec (dat file) and simulate:
  readDatAndSimulate(/*datFile=*/ argv[1], /*outFile=*/ argv[3]);

  return 0;
}

unsigned int getRandSeed() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned int ret = tv.tv_sec + tv.tv_usec * 100000;
  return ret;
}


// Reads .dat file from <in> and simulates the chromosomes according to each
// line in the .dat file.
void readDatAndSimulate(char *datFile, char *outFile) {
  const int BUF_SIZE = 80;
  char buf[BUF_SIZE];
  char c;
  // Arrays with indexes corresponding to each simulated sample.  For each
  // sample, we store a list of segments that contain informatino on the source
  // population for the segment and the end marker number (start is implicit)
  //
  // We have two lists since the admixed samples from one generation can mix to
  // form subsequent generations and thus we store the values from the previous
  // simulation for use during the next simulation
  vector<SegmentList> simu1, simu2;
  vector<SegmentList> *simuOutput = &simu1;
  vector<SegmentList> *prevSimulated = &simu2;

  // open dat file:
  FILE *in = fopen(datFile, "r");
  if (!in) {
    printf("Error, could not open dat file %s!\n", datFile);
    exit(1);
  }

  int numSamples;
  int numRead = fscanf(in, "%d", &numSamples);
  if (numRead != 1) {
    printf("Error reading number of samples.\n");
    exit(1);
  }

  // Simulate recombination break points for at least 10k samples (or 20*N)
  // haplotypes so that we're sure to get a diverse population to sample from
  // in each generation
  int numSampsToSimulate = max(10000, 20 * numSamples);

  // Make space in arrays that will store simulation results for each sample
  for(int i = 0; i < numSampsToSimulate; i++) {
    vector<Segment> *a = new vector<Segment>();
    simu1.push_back(a);
    a = new vector<Segment>();
    simu2.push_back(a);
  }


  // compulsory label for the admixture proportion; can ignore
  int ret = fscanf(in, "%*s");
  assert(ret == 0);

  // read the population labels
  while(true) {
    // read until next non-white space character; stop at newline (ends header)
    while (isspace(c = fgetc(in)) && c != '\n');
    if (c == '\n' || c == EOF)
      break;  // end of header
    ungetc(c, in);

    int len = 0;
    while (len < BUF_SIZE - 1 && !isspace(c = fgetc(in)) && c != EOF) {
      buf[len] = c;
      len++;
    }
    ungetc(c, in);
    buf[len] = '\0'; // null terminate the string
    char *newName = new char[len+1];
    strcpy(newName, buf);
    popNames.push_back(newName);
  }

  // The generation number that was last simulated; initially starting from
  // generation 0:
  int prevGeneration = 0;
  int line = 2; // what line number are we on?
  int numPops = popNames.size() + 1; // + 1 for admixed pop
  double *popProportions = new double[numPops];

  // read each line of simulation specs and perform the simulations by calling
  // simulate()
  while (true) {

    // read up to first non-whitespace character; bug out if EOF
    while (isspace(c = fgetc(in)));
    if (c == EOF)
      break;  // done simulating
    ungetc(c, in);

    // read number of generations to simulate
    int simToGeneration;
    numRead = fscanf(in, "%d", &simToGeneration);
    if (numRead != 1) {
      printf("Error on line %d: failed to read the generation number.\n",
	     line);
      exit(1);
    }

    if (simToGeneration <= prevGeneration) {
      printf("Error on line %d: generation %d is not greater than the ",
	     line, simToGeneration);
      printf("previous line.\n");
      exit(1);
    }

    printf("Simulating: Generation: %2d, ", simToGeneration);

    double totalProportions = 0.0f;

    // read in proportions for each of the populations
    for(int i = 0; i < numPops; i++) {
      numRead = fscanf(in, "%lf", &popProportions[i]);
      if (numRead != 1) {
	printf("Error on line %d: failed to read proportion for ", line);
	if (i == 0) {
	  printf("admixed population.\n");
	}
	else {
	  printf("population %s.\n", popNames[i-1]);
	}
	exit(1);
      }
      totalProportions += popProportions[i];

      if (i == 0)
	printf("Admixed: ");
      else
	printf("%s: ", popNames[i-1]);
      printf("%.3f ", popProportions[i]);
    }

    if (totalProportions != 1.0) {
      printf("Error on line %d: population proportions do not sum to 1.0\n",
	     line);
      exit(1);
    }

    // read to end of line: there shouldn't be any more non-whitespace chars
    while (isspace(c = fgetc(in)) && c != '\n');
    if (c != '\n' && c != EOF) {
      printf("Error on line %d: more than %d populations specified",
	  line, numPops);
      exit(1);
    }

    printf("\n");

    int simGenerations = simToGeneration - prevGeneration;
    // simulate the first generation:
    simulate(numSampsToSimulate, popProportions, numPops, *simuOutput,
	     *prevSimulated, prevGeneration);
    for(int i = 1; i < simGenerations; i++) {
      // ready the lists for this new simulation:
      swapSimuListsAndClear(prevSimulated, simuOutput, numSampsToSimulate);
      // simulate the rest of the generations, which will only descend from the
      // admixed samples generated in the previous simulate() call.
      simulate(numSampsToSimulate, popProportions, numPops, *simuOutput,
	       *prevSimulated, prevGeneration+i, /*useOnlyAdmixed=*/ true);
    }

    if (c == EOF)
      break; // done simulating

    ////////////////////////////////////////////////////////////////////////
    // prepare for next line of simulation instructions
    prevGeneration = simToGeneration;
    line++;
    
    // ready the lists for next simulation:
    swapSimuListsAndClear(prevSimulated, simuOutput, numSampsToSimulate);
  }

  fclose(in);

  // done: print the simulation results file:
  output(outFile, *prevSimulated, numSamples);
}

void simulate(int numSampsToSimulate, double *popProportions, int numPops,
	      vector<SegmentList> &simuOutput,
	      vector<SegmentList> &prevSimulated, int prevGeneration,
	      bool useOnlyAdmixed) {
  for(int ind = 0; ind < numSampsToSimulate; ind++) {

    int prevChrom = -1;
    double prevMapPos = 0.0;

    // get populations for the diploid parent of the haploid individual being
    // simulated. Will randomly choose populations for both chromosomes unless
    // we're only sampling admixed haplotypes or if we're in the first
    // generation (see below)
    int copyPop = 0; // initially assume <useOnlyAdmixed> is true...
    if (!useOnlyAdmixed)
      copyPop = choosePop(popProportions, numPops);

    // for when a homolog(s) is from the admixed population in the previous
    // generation, get chromosome number(s) to copy from; require that these
    // are different if copyPop == 0 (both from prev generation
    // admixture)
    int copyChrom[2] = { -1, -1 }; // initially assume copyPop != 0
    if (copyPop == 0) {
      for(int i = 0; i < 2; i++) {
	copyChrom[i] = rand() % numSampsToSimulate;
      }
    }
    while (copyChrom[0] == copyChrom[1]) { // make sure these are different
      copyChrom[1] = rand() % numSampsToSimulate;
    }

    // choose one of the homologs to copy from
    int curHomolog = rand() % 2;

    int numMarkers = Marker::getNumMarkers();
    for(int marker = 0; marker < numMarkers; marker++) {
      const Marker *curMarker = Marker::getMarker(marker);
      int curChrom = curMarker->getChrom();
      double curMapPos = curMarker->getMapPos();
      if (curChrom != prevChrom) {
	// End of prev chrom: record the segment (if we've sampled: marker > 0)
	if (marker > 0)
	  recordSegment(simuOutput[ind], copyPop, copyChrom[curHomolog],
			/*endMarker=*/ marker - 1, prevSimulated,
			numSampsToSimulate);

	// resample copy homolog for new chromosome
	curHomolog = rand() % 2;

	// no previous marker on this chromosome, so no genetic distance -- skip
	prevChrom = curChrom;
	prevMapPos = curMapPos;
	continue;
      }

      double geneticDistance = curMapPos - prevMapPos;
      bool isRecomb = decideIfRecomb(geneticDistance);
      if (isRecomb) {
	// end of the currently running segment; record it:
	// Note: we use endMarker = marker - 1 since the recombination takes
	// place between the current marker and the previous one.
	recordSegment(simuOutput[ind], copyPop, copyChrom[curHomolog],
		      /*endMarker=*/ marker - 1, prevSimulated,
		      numSampsToSimulate);

	// switch the homolog that we're copying from
	curHomolog ^= 1;

	// don't need to resample now that we're only sampling whole chromosomes
	// from a population during the first generation and then only sampling
	// from these chromosomes (which become the admixed population) in
	// subsequent generations
//	curPop = choosePop(popProportions, numPops);
      }
      
      prevMapPos = curMapPos;
      prevChrom = curChrom;
    }
    // record final segment
    recordSegment(simuOutput[ind], copyPop, copyChrom[curHomolog],
		  /*endMarker=*/ numMarkers - 1, prevSimulated,
		  numSampsToSimulate);

    if ((ind+1) % 128 == 0) {
      printf("Simulated individual %d / %d in generation %d\r",
	     ind+1, numSampsToSimulate, prevGeneration + 1);
      fflush(stdout);
    }
  }

  printf("Simulated all %d individuals in generation %d            \n",
	 numSampsToSimulate, prevGeneration + 1);
}

bool decideIfRecomb(double geneticDistance) {
  double randVal = (double) rand() / RAND_MAX;
  double recombProb = 1 - exp(-geneticDistance);
  return randVal < recombProb;
}

void recordSegment(SegmentList outputInd, int popNum, int copyChrom,
		   int endMarker, vector<SegmentList> &prevSimulated,
		   int numSampsToSimulate) {
  if (popNum != 0) {
    Segment s;
    s.popNum = popNum;
    s.endMarker = endMarker;
    outputInd->push_back(s);
  }
  else {// popNum == 0: Sample from one of the previously sampled admixed chroms

    // Determine starting marker:
    int startMarker;
    int numSimuSegments = outputInd->size();
    if (numSimuSegments == 0)
      startMarker = 0;
    else
      startMarker = (*outputInd)[numSimuSegments - 1].endMarker + 1;

    assert(copyChrom >= 0 && copyChrom < numSampsToSimulate);

    vector<Segment>::const_iterator copySeg;
    for(copySeg = prevSimulated[copyChrom]->begin();
	copySeg != prevSimulated[copyChrom]->end(); copySeg++) {
      if (copySeg->endMarker >= startMarker) {
	break;  // found the segment to start copying from
      }
    }

    assert(copySeg != prevSimulated[copyChrom]->end());

    while(copySeg->endMarker < endMarker) {
      // copy from this previously admixed sample
      Segment s = *copySeg;
      outputInd->push_back(s);
      copySeg++;

      assert(copySeg != prevSimulated[copyChrom]->end());
    }

    // Now copy the last segment (so we get all values <= endMarker)
    Segment s = *copySeg;
    s.endMarker = endMarker; // only copy to endMarker
    outputInd->push_back(s);
  }
}

int choosePop(double *popProportions, int numPops) {
  double randVal = (double) rand() / RAND_MAX;

  double totalPrevProportions = 0.0f;
  for(int i = 0; i < numPops; i++) {
    if (popProportions[i] == 0.0f) {
      continue; // no contribution from this population
    }

    if (randVal <= popProportions[i] + totalPrevProportions) {
      // use population i
      return i;
    }
    totalPrevProportions += popProportions[i];
  }

  printf("Error: got to end of choosePop() function without population!\n");
  exit(10);

  return 0;
}

void swapSimuListsAndClear(vector<SegmentList> *&prevSimulated,
			   vector<SegmentList> *&simuOutput,
			   int numSampsToSimulate) {
  vector<SegmentList> *tmp = prevSimulated;
  prevSimulated = simuOutput;
  simuOutput = tmp;
  for(int i = 0; i < numSampsToSimulate; i++) 
    // output list to be updated next round
    (*simuOutput)[i]->clear();
}

void output(char *outFile, vector<SegmentList> &simulated, int numSamples) {
  FILE *out = fopen(outFile, "w");
  if (!out) {
    printf("Error, could not open output file %s!\n", outFile);
    exit(1);
  }

  // print out the population labels:
  vector<char *>::const_iterator i;
  for(i = popNames.begin(); i != popNames.end(); i++) {
    fprintf(out, "%s ", *i);
  }
  fprintf(out, "\n");

  // We've simulated double the number of samples we need (see
  // numSampsToSimulate), so we'll randomly sample from among these which
  // individuals to print; the following records which individuals have been
  // printed
  int numSimulated = simulated.size();
  bool *samplePrinted = new bool[numSimulated];
  for(int i = 0; i < numSimulated; i++) {
    samplePrinted[i] = false;
  }

  // print <numSamples> samples, randomly choosing one of the simulated chroms
  // to print
  for(int i = 0; i < numSamples; i++) {
    int chrom = rand() % numSimulated;
    while(samplePrinted[chrom]) {
      chrom = rand() % numSimulated;
    }

    samplePrinted[chrom] = true; // will now print:

    int numSegments = simulated[chrom]->size();

    for(int j = 0; j < numSegments; j++) {
      Segment s = (*simulated[chrom])[j];
      fprintf(out, "%d:%d ", s.popNum-1, s.endMarker);
    }
    fprintf(out, "\n");
  }

  fclose(out);
}

void printUsage(char **argv) {
  printf("Usage:\n");
  printf("  %s [in.dat] [in.snp] [out.bp]\n\n", argv[0]);
  printf("Note: [in.snp] is expected to be in Eigenstrat format with marker ");
  printf("distances\nspecified in Morgans\n");
  exit(1);
}
