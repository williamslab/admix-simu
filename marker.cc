#include <stdio.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include "marker.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// initialize static members
vector<Marker *> Marker::_allMarkers;
int Marker::_firstMarkerNum[LAST_CHROM + 1] = { 0 };

void Marker::readMarkerFile(char *markerFile) {
  FILE *in;
  char markerName[256];
  int chrom;
  Marker *prevMarker = NULL;
  double mapPos;
  int physPos;
  char alleles[2];

  assert(_allMarkers.size() == 0); // shouldn't have any markers yet
  _allMarkers.reserve(600000);     // allocate space for 600k markers

  int numMarkersCurChrom = 0;

  in = fopen(markerFile, "r");
  // Note: I assume the map positions are in Morgans.
  while (fscanf(in, "%s %d %lf %d %c %c", markerName, &chrom, &mapPos, &physPos,
		&alleles[0], &alleles[1]) == 6) {
    Marker *m = new Marker(markerName, chrom, mapPos, physPos, alleles);

    if (prevMarker != NULL) {
      int prevChrom = prevMarker->_chrom;
      if (chrom == prevChrom) {
	assert(mapPos >= prevMarker->_mapPos && physPos > prevMarker->_physPos);
      }
      else  { // Have a valid prev chrom?  Update marker counts
	assert(chrom > prevChrom);

	numMarkersCurChrom = 0; // reset

	// about to append the first marker for this chrom to this index:
	_firstMarkerNum[chrom] = _allMarkers.size();
      }
    }

    _allMarkers.push_back(m);
    numMarkersCurChrom++;
    prevMarker = m;
  }

  fclose(in);
}

Marker::Marker(char *markerName, int chrom, double mapPos, int physPos,
	       char alleles[2]) {
  _name = new char[strlen(markerName)+1];
  strcpy(_name, markerName);
  _chrom = chrom;
  _mapPos = mapPos;
  _physPos = physPos;
  _alleles[0] = alleles[0];
  _alleles[1] = alleles[1];

  // No longer necessary:
  // These values are calculated as the genotypes are read in, using the method
  // Marker::addMarkerGenotype()
//  _majAlleleCount = _totalAlleles = 0;
}


// No longer necessary:
//void Marker::addMarkerGenotype(int markerNum, int genotype) {
//  if (genotype != 9) { // as long as it's not missing data:
//    // note: a value of 0 is two major alleles, and a genotype of 2 is 0 major
//    // alleles
//    Marker *theMarker = _allMarkers[markerNum];
//    theMarker->_majAlleleCount += 2 - genotype;
//    theMarker->_totalAlleles += 2;
//  }
//}
//
//double Marker::getAlleleFreq(int markerNum) {
//  Marker *theMarker = _allMarkers[markerNum];
//  return (double) theMarker->_majAlleleCount / theMarker->_totalAlleles;
//}
